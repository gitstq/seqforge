"""Open reading frame (ORF) finding on both strands."""

from .alphabet import normalize, reverse_complement
from .genetic import get_table, translate


def find_orfs(seq, min_len=30, table=1, both_strands=True, starts=None):
    """Find ORFs in a nucleic acid sequence.

    Returns a list of dictionaries with keys ``strand``, ``frame``,
    ``start``, ``end`` (0-based, half-open on the *input* strand),
    ``length`` and ``protein``. Nested ORFs are reported individually so
    callers can filter as needed.
    """
    seq = normalize(seq)
    if starts is None:
        starts = ("ATG",)
    stops = set(codon for codon, aa in get_table(table).items() if aa == "*")
    starts = set(normalize(s) for s in starts)
    results = []
    targets = [("+", seq)]
    if both_strands:
        targets.append(("-", reverse_complement(seq)))
    for strand, target in targets:
        for frame in range(3):
            index = frame
            while index + 3 <= len(target):
                codon = target[index:index + 3]
                if codon in starts:
                    cursor = index
                    stop_at = None
                    while cursor + 3 <= len(target):
                        if target[cursor:cursor + 3] in stops:
                            stop_at = cursor
                            break
                        cursor += 3
                    if stop_at is None:
                        break  # no stop codon in the rest of this frame
                    orf = target[index:stop_at + 3]
                    if len(orf) >= min_len:
                        if strand == "+":
                            start, end = index, stop_at + 3
                        else:
                            length = len(seq)
                            start, end = length - (stop_at + 3), length - index
                        results.append({
                            "strand": strand,
                            "frame": frame + 1,
                            "start": start,
                            "end": end,
                            "length": len(orf),
                            "protein": translate(orf, table=table),
                        })
                    index += 3  # keep scanning for alternative starts
                else:
                    index += 3
    results.sort(key=lambda item: (item["start"], item["end"]))
    return results


def longest_orf(seq, min_len=30, table=1, both_strands=True):
    """Return the longest ORF dictionary (or None when none is found)."""
    orfs = find_orfs(seq, min_len=min_len, table=table,
                     both_strands=both_strands)
    if not orfs:
        return None
    return max(orfs, key=lambda item: item["length"])
