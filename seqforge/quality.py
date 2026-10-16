"""FASTQ quality handling: Phred offsets, trimming and filtering."""

from .alphabet import normalize


def detect_phred_offset(qual):
    """Guess the Phred offset (33 or 64) of a quality string.

    Heuristic: any character below ``@`` (64) can only be Phred+33;
    characters above ``J`` (74) can only be Phred+64; otherwise Phred+33
    (the modern default) is assumed.
    """
    if not qual:
        return 33
    lowest = min(ord(ch) for ch in qual)
    highest = max(ord(ch) for ch in qual)
    if lowest < 64:
        return 33
    if highest > 74:
        return 64
    return 33


def phred_scores(qual, offset=None):
    """Convert a quality string to integer Phred scores."""
    if offset is None:
        offset = detect_phred_offset(qual)
    return [ord(ch) - offset for ch in qual]


def mean_quality(qual, offset=None):
    """Mean Phred score of a quality string (0.0 when empty)."""
    scores = phred_scores(qual, offset)
    if not scores:
        return 0.0
    return round(sum(scores) / float(len(scores)), 2)


def per_base_quality(records):
    """Mean quality per base position, up to the shortest record.

    Returns a list of floats (position 0 first).
    """
    if not records:
        return []
    offset = None
    min_len = min(len(r.qual) for r in records)
    if min_len == 0:
        return []
    columns = []
    for position in range(min_len):
        column = []
        for record in records:
            if offset is None:
                offset = detect_phred_offset(record.qual)
            column.append(ord(record.qual[position]) - offset)
        columns.append(round(sum(column) / float(len(column)), 2))
    return columns


def trim_sliding(scores, threshold, window):
    """Return ``(start, end)`` of the best region under sliding-window QC.

    Every window of ``window`` bases fully inside the region must have a
    mean score >= ``threshold`` (BWA-style quality trimming). Returns
    ``(0, 0)`` when no base qualifies.
    """
    if not scores:
        return 0, 0
    if window <= 0:
        window = 1
    if len(scores) <= window:
        mean = sum(scores) / float(len(scores))
        return (0, len(scores)) if mean >= threshold else (0, 0)
    window_means = []
    for start in range(0, len(scores) - window + 1):
        chunk = scores[start:start + window]
        window_means.append(sum(chunk) / float(window))
    best_start, best_end = 0, 0
    run_start = None
    for index, ok in enumerate(window_means):
        if ok and run_start is None:
            run_start = index
        elif not ok and run_start is not None:
            span_start, span_end = run_start, index - 1 + window
            if span_end - span_start > best_end - best_start:
                best_start, best_end = span_start, span_end
            run_start = None
    if run_start is not None:
        span_start = run_start
        span_end = len(window_means) - 1 + window
        if span_end - span_start > best_end - best_start:
            best_start, best_end = span_start, span_end
    return best_start, best_end


def trim_record(record, threshold=20, window=4):
    """Trim a FASTQ record with the sliding-window algorithm."""
    scores = phred_scores(record.qual)
    start, end = trim_sliding(scores, threshold, window)
    record.seq = record.seq[start:end]
    record.qual = record.qual[start:end]
    return record


def n_fraction(seq):
    """Fraction of ``N`` bases (0.0 for empty sequences)."""
    if not seq:
        return 0.0
    return normalize(seq).count("N") / float(len(seq))


def passes_filters(record, min_len=0, max_len=None, min_mean_q=0.0,
                   max_n_frac=1.0):
    """True when a record passes length / quality / N filters."""
    length = len(record.seq)
    if length < min_len:
        return False
    if max_len is not None and length > max_len:
        return False
    if min_mean_q > 0 and mean_quality(record.qual) < min_mean_q:
        return False
    if n_fraction(record.seq) > max_n_frac:
        return False
    return True


def filter_fastq(records, min_len=0, max_len=None, min_mean_q=0.0,
                 max_n_frac=1.0):
    """Yield records that pass all filters."""
    for record in records:
        if passes_filters(record, min_len, max_len, min_mean_q, max_n_frac):
            yield record
