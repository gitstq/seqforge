"""Motif search (IUPAC aware), restriction sites and primer thermodynamics."""

import math
import re

from .alphabet import normalize, reverse_complement

#: IUPAC nucleotide code -> regex fragment.
IUPAC_REGEX = {
    "A": "A", "C": "C", "G": "G", "T": "T", "U": "U",
    "R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]",
    "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]",
    "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]",
}

#: Common restriction enzymes: name -> recognition sequence.
RESTRICTION_ENZYMES = {
    "EcoRI": "GAATTC",
    "BamHI": "GGATCC",
    "HindIII": "AAGCTT",
    "NotI": "GCGGCCGC",
    "XhoI": "CTCGAG",
    "SalI": "GTCGAC",
    "PstI": "CTGCAG",
    "SmaI": "CCCGGG",
    "KpnI": "GGTACC",
    "SacI": "GAGCTC",
    "NcoI": "CCATGG",
    "NdeI": "CATATG",
    "XbaI": "TCTAGA",
    "SpeI": "ACTAGT",
    "EcoRV": "GATATC",
    "PvuII": "CAGCTG",
    "SphI": "GCATGC",
    "NheI": "GCTAGC",
    "AgeI": "ACCGGT",
    "BglII": "AGATCT",
}


def motif_to_regex(motif):
    """Convert an IUPAC DNA motif into a case-insensitive regex."""
    parts = []
    for position, char in enumerate(motif):
        fragment = IUPAC_REGEX.get(char.upper())
        if fragment is None:
            raise ValueError(
                "invalid IUPAC code %r in motif at position %d"
                % (char, position + 1)
            )
        parts.append(fragment)
    return re.compile("".join(parts), re.IGNORECASE)


def find_motif(seq, motif, both_strands=False):
    """Find all occurrences of *motif* in *seq*.

    Returns a list of ``(start, end, strand)`` tuples (0-based, half-open).
    ``both_strands`` additionally matches the reverse complement.
    """
    pattern = motif_to_regex(motif)
    matches = []
    for match in pattern.finditer(seq):
        matches.append((match.start(), match.end(), "+"))
    if both_strands:
        rc_pattern = motif_to_regex(reverse_complement(motif))
        for match in rc_pattern.finditer(seq):
            matches.append((match.start(), match.end(), "-"))
    matches.sort()
    return matches


def restriction_sites(seq, enzymes=None):
    """Locate recognition sites for (palindromic) restriction enzymes.

    ``enzymes`` is an iterable of names or ``None`` for all built-ins.
    Returns a list of ``(enzyme, start, end)`` tuples sorted by position.
    """
    seq = normalize(seq)
    names = sorted(enzymes) if enzymes else sorted(RESTRICTION_ENZYMES)
    sites = []
    for name in names:
        recognition = RESTRICTION_ENZYMES.get(name)
        if recognition is None:
            raise ValueError(
                "unknown enzyme %r (use restriction_enzymes() to list "
                "built-ins)" % name
            )
        pattern = motif_to_regex(recognition)
        for match in pattern.finditer(seq):
            sites.append((name, match.start(), match.end()))
    sites.sort(key=lambda item: (item[1], item[0]))
    return sites


def restriction_enzymes():
    """Return the built-in enzyme table (name -> recognition sequence)."""
    return dict(RESTRICTION_ENZYMES)


def gc_percent(primer):
    """GC percentage of a primer (A/C/G/T/U only)."""
    primer = normalize(primer)
    if not primer:
        return 0.0
    gc = sum(primer.count(base) for base in "GC")
    return round(100.0 * gc / len(primer), 2)


def tm_wallace(primer):
    """Wallace rule melting temperature: ``2*(A+T) + 4*(G+C)`` in degC."""
    primer = normalize(primer)
    at = primer.count("A") + primer.count("T")
    gc = primer.count("G") + primer.count("C")
    return 2 * at + 4 * gc


def tm_salt_adjusted(primer, na=0.05):
    """Salt-adjusted Tm for primers of ~14-60 nt.

    Uses the common approximation
    ``Tm = 81.5 + 0.41*%%GC - 675/N + 16.6*log10([Na+])`` where ``[Na+]``
    is in mol/L. Documented as an estimate, not a substitute for
    nearest-neighbour models.
    """
    primer = normalize(primer)
    length = len(primer)
    if length == 0:
        return 0.0
    if na <= 0:
        raise ValueError("sodium concentration must be positive")
    return round(
        81.5 + 0.41 * gc_percent(primer) - 675.0 / length
        + 16.6 * math.log10(na), 1,
    )


def primer_stats(primer, na=0.05):
    """Return length / GC / Tm statistics for one primer."""
    primer = normalize(primer)
    return {
        "primer": primer,
        "length": len(primer),
        "gc": gc_percent(primer),
        "tm_wallace": tm_wallace(primer),
        "tm_adjusted": tm_salt_adjusted(primer, na=na) if primer else 0.0,
    }
