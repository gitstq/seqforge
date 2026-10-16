"""Sequence statistics: lengths, N50/L50, GC content and GC skew."""

import statistics


def length_stats(lengths):
    """Return a summary dictionary for a list of sequence lengths.

    Includes count, total, min, max, mean, median, population standard
    deviation, N50 and L50. Empty input yields zeroed statistics.
    """
    if not lengths:
        return {
            "count": 0, "total": 0, "min": 0, "max": 0, "mean": 0.0,
            "median": 0.0, "stdev": 0.0, "n50": 0, "l50": 0,
        }
    ordered = sorted(lengths, reverse=True)
    total = sum(ordered)
    n50 = 0
    l50 = 0
    cumulative = 0
    for position, length in enumerate(ordered, start=1):
        cumulative += length
        if cumulative * 2 >= total:
            n50 = length
            l50 = position
            break
    return {
        "count": len(ordered),
        "total": total,
        "min": min(ordered),
        "max": max(ordered),
        "mean": round(statistics.mean(ordered), 2),
        "median": round(statistics.median(ordered), 2),
        "stdev": round(statistics.pstdev(ordered), 2),
        "n50": n50,
        "l50": l50,
    }


def _counts(seq):
    upper = seq.upper()
    return {
        "A": upper.count("A"), "C": upper.count("C"),
        "G": upper.count("G"), "T": upper.count("T"),
        "U": upper.count("U"), "N": upper.count("N"),
    }


def gc_fraction(seq):
    """Return ``(gc_count, at_count, called_count)`` for a sequence.

    ``called_count`` counts only unambiguous A/C/G/T/U bases, so callers can
    weight correctly even when the input contains ``N`` or other ambiguity
    codes.
    """
    upper = seq.upper()
    gc = upper.count("G") + upper.count("C")
    at = upper.count("A") + upper.count("T") + upper.count("U")
    return gc, at, gc + at


def gc_content(seq):
    """GC percentage over unambiguous bases.

    RNA input is handled: ``U`` counts as the A/T-family base. Ambiguity
    codes such as ``N`` are excluded from the denominator; a sequence with
    no called bases yields ``0.0``.
    """
    gc, at, called = gc_fraction(seq)
    if called == 0:
        return 0.0
    return round(100.0 * gc / called, 2)


def gc_skew(seq):
    """GC skew ``(G - C) / (G + C)``; 0.0 when G + C is zero."""
    counts = _counts(seq)
    g, c = counts["G"], counts["C"]
    if g + c == 0:
        return 0.0
    return round((g - c) / float(g + c), 4)


def sliding_gc(seq, window=100, step=50):
    """Yield ``(start, gc_percent)`` tuples over sliding windows.

    ``window``/``step`` are in bases; the final partial window is skipped.
    A non-positive ``step`` falls back to the window size.
    """
    if window <= 0:
        raise ValueError("window size must be positive")
    if step <= 0:
        step = window
    seq = seq.upper()
    for start in range(0, len(seq) - window + 1, step):
        chunk = seq[start:start + window]
        yield start, gc_content(chunk)


def sequence_summary(records):
    """Summarise a list of FASTA/FASTQ records.

    Returns ``(summary, per_sequence)`` where ``summary`` aggregates the
    whole file and ``per_sequence`` is a list of per-record dictionaries.
    """
    per_sequence = []
    lengths = []
    gc_bases = 0
    called_bases = 0
    for record in records:
        length = len(record.seq)
        lengths.append(length)
        gc, at, called = gc_fraction(record.seq)
        gc_bases += gc
        called_bases += called
        per_sequence.append({
            "id": record.id,
            "length": length,
            "gc": round(100.0 * gc / called, 2) if called else 0.0,
        })
    summary = length_stats(lengths)
    summary["gc"] = (
        round(100.0 * gc_bases / called_bases, 2) if called_bases else 0.0
    )
    return summary, per_sequence
