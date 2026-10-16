"""IUPAC nucleotide alphabets and complement mappings.

All helpers are pure-Python and dependency free. Unknown characters raise
``ValueError`` so that malformed data fails loudly instead of silently
producing wrong biology.
"""

DNA_IUPAC = frozenset("ACGTURYSWKMBDHVN")
RNA_IUPAC = frozenset("ACGURYSWKMBDHVN")
PROTEIN_IUPAC = frozenset("ACDEFGHIKLMNPQRSTVWYBXZJUO*-")

#: Case-insensitive DNA -> DNA complement (IUPAC ambiguity codes included).
DNA_COMPLEMENT = {
    "A": "T", "T": "A", "G": "C", "C": "G", "U": "A", "N": "N",
    "R": "Y", "Y": "R", "S": "S", "W": "W", "K": "M", "M": "K",
    "B": "V", "V": "B", "D": "H", "H": "D",
    "-": "-", ".": ".",
}

#: Case-insensitive RNA -> RNA complement (IUPAC ambiguity codes included).
RNA_COMPLEMENT = {
    "A": "U", "U": "A", "T": "A", "G": "C", "C": "G", "N": "N",
    "R": "Y", "Y": "R", "S": "S", "W": "W", "K": "M", "M": "K",
    "B": "V", "V": "B", "D": "H", "H": "D",
    "-": "-", ".": ".",
}


def reverse_complement(seq, rna=False):
    """Return the reverse complement of a nucleic acid sequence.

    Case is preserved. Ambiguity codes (IUPAC) are complemented correctly.
    Characters outside the IUPAC alphabets raise ``ValueError``.
    """
    table = RNA_COMPLEMENT if rna else DNA_COMPLEMENT
    out = []
    for pos, ch in enumerate(seq):
        comp = table.get(ch)
        if comp is None:
            comp = table.get(ch.upper())
        if comp is None:
            raise ValueError(
                "invalid nucleotide %r at position %d (not IUPAC)" % (ch, pos + 1)
            )
        out.append(comp)
    return "".join(reversed(out))


def complement(seq, rna=False):
    """Return the complement (not reversed) of a nucleic acid sequence."""
    table = RNA_COMPLEMENT if rna else DNA_COMPLEMENT
    out = []
    for pos, ch in enumerate(seq):
        comp = table.get(ch)
        if comp is None:
            comp = table.get(ch.upper())
        if comp is None:
            raise ValueError(
                "invalid nucleotide %r at position %d (not IUPAC)" % (ch, pos + 1)
            )
        out.append(comp)
    return "".join(out)


def is_nucleic(seq):
    """True when every non-gap character is an IUPAC nucleotide code."""
    return all(ch in DNA_IUPAC for ch in seq.upper() if ch not in "-. ")


def is_rna(seq):
    """True when the sequence contains U and no T (heuristic)."""
    upper = seq.upper()
    return "U" in upper and "T" not in upper


def is_protein(seq):
    """True when every non-gap character is a valid amino acid code."""
    return all(ch in PROTEIN_IUPAC for ch in seq.upper() if ch not in "-. ")


def normalize(seq):
    """Upper-case a sequence and convert RNA ``U`` to DNA ``T``."""
    return seq.upper().replace("U", "T")
