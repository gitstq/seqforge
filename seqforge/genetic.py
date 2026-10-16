"""Genetic code tables, translation and transcription.

Amino-acid strings are written in the classic ``TCAG`` codon order
(base1 x base2 x base3, each cycling T, C, A, G), which makes them easy to
verify against the NCBI tables.
"""

from .alphabet import normalize, reverse_complement

#: NCBI genetic code table id -> human readable name.
GENETIC_CODES = {
    1: "Standard",
    2: "Vertebrate Mitochondrial",
    3: "Yeast Mitochondrial",
    5: "Invertebrate Mitochondrial",
    6: "Ciliate Nuclear",
}

_BASES = "TCAG"

_STANDARD_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
# Table 2: ATA->M, AGA/AGG->*
_VERTEBRATE_MT_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG"
# Table 3: CTN->T, ATA->M, TGA->W
_YEAST_MT_AA = "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
# Table 5: AGA/AGG->S, ATA->M, TGA->W
_INVERTEBRATE_MT_AA = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG"
# Table 6: TAA/TAG->Q
_CILIATE_AA = "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"

_AA_STRINGS = {
    1: _STANDARD_AA,
    2: _VERTEBRATE_MT_AA,
    3: _YEAST_MT_AA,
    5: _INVERTEBRATE_MT_AA,
    6: _CILIATE_AA,
}


def _build_table(aa_string):
    table = {}
    index = 0
    for base1 in _BASES:
        for base2 in _BASES:
            for base3 in _BASES:
                table[base1 + base2 + base3] = aa_string[index]
                index += 1
    return table


def get_table(table=1):
    """Return the codon -> amino-acid dictionary for an NCBI table id."""
    aa_string = _AA_STRINGS.get(table)
    if aa_string is None:
        raise ValueError(
            "unsupported genetic code table %r (supported: %s)"
            % (table, ", ".join(str(k) for k in sorted(_AA_STRINGS)))
        )
    return _build_table(aa_string)


def stop_codons(table=1):
    """Return the sorted stop codons of a genetic code table."""
    codons = get_table(table)
    return sorted(codon for codon, aa in codons.items() if aa == "*")


def start_codons(table=1):
    """Return start codons for a table (ATG, plus alternative starts).

    SeqForge deliberately ORFs only on ATG by default; this helper exposes
    the common alternative initiators for callers that need them.
    """
    if table in (2, 3, 5):
        return ["ATG", "ATA"]
    if table in (4,):
        return ["ATG", "GTG", "TTG"]
    return ["ATG"]


def transcribe(dna):
    """Transcribe DNA to RNA (T -> U)."""
    return normalize(dna).replace("T", "U")


def reverse_transcribe(rna):
    """Reverse-transcribe RNA to DNA (U -> T)."""
    return rna.upper().replace("U", "T")


def translate(seq, table=1, to_stop=False, trim_stop=True):
    """Translate a nucleic acid sequence to protein.

    * RNA input (U) is accepted and converted to DNA first.
    * Codons containing ambiguity codes translate to ``X``.
    * ``to_stop`` halts translation at the first stop codon.
    * ``trim_stop`` removes a single trailing ``*`` (ignored when to_stop).
    """
    codons = get_table(table)
    seq = normalize(seq)
    protein = []
    for start in range(0, len(seq) - len(seq) % 3, 3):
        codon = seq[start:start + 3]
        if len(codon) < 3:
            break
        if all(ch in "ACGT" for ch in codon):
            aa = codons.get(codon, "X")
        else:
            aa = "X"
        if aa == "*" and to_stop:
            protein.append("*")
            break
        protein.append(aa)
    result = "".join(protein)
    if trim_stop and not to_stop:
        if result.endswith("*"):
            result = result[:-1]
    return result


def translate_frames(seq, table=1, both_strands=False):
    """Translate a sequence in all 3 (or 6) reading frames.

    Returns a list of ``(frame_label, protein)`` tuples where labels are
    ``+1``, ``+2``, ``+3`` and (optionally) ``-1``, ``-2``, ``-3``.
    """
    seq = normalize(seq)
    frames = [("+1", seq), ("+2", seq[1:]), ("+3", seq[2:])]
    if both_strands:
        rc = reverse_complement(seq)
        frames += [("-1", rc), ("-2", rc[1:]), ("-3", rc[2:])]
    return [(label, translate(s, table=table)) for label, s in frames]
