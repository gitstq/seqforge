"""Format conversion: FASTA/FASTQ <-> CSV/TSV/NDJSON."""

import csv
import io
import json

from .stats import gc_content
from .seqio import write_fasta, write_fastq


def fasta_to_fastq(records, default_qual="I"):
    """Convert FASTA records to FASTQ records.

    ``default_qual`` is the placeholder quality character (``I`` = Phred33
    Q40 by default) since FASTA carries no quality information.
    """
    from .seqio import FastqRecord

    converted = []
    for record in records:
        qual = default_qual * len(record.seq)
        converted.append(
            FastqRecord(record.id, record.description, record.seq, qual)
        )
    return converted


def fastq_to_fasta(records):
    """Convert FASTQ records to FASTA records (quality is dropped)."""
    from .seqio import FastaRecord

    return [
        FastaRecord(record.id, record.description, record.seq)
        for record in records
    ]


def record_rows(records):
    """Build plain rows (id, description, length, gc, sequence)."""
    rows = []
    for record in records:
        rows.append([
            record.id,
            record.description,
            len(record.seq),
            gc_content(record.seq),
            record.seq,
        ])
    return rows


ROW_HEADERS = ["id", "description", "length", "gc", "sequence"]


def to_csv(records, path=None):
    """Write records as CSV (or return the CSV text when path is None)."""
    if path is None:
        buffer = io.StringIO()
        writer = csv.writer(buffer, lineterminator="\n")
        writer.writerow(ROW_HEADERS)
        writer.writerows(record_rows(records))
        return buffer.getvalue()
    with io.open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(ROW_HEADERS)
        writer.writerows(record_rows(records))
    return None


def to_tsv(records, path=None):
    """Write records as TSV (or return the TSV text when path is None)."""
    if path is None:
        buffer = io.StringIO()
        writer = csv.writer(buffer, delimiter="\t", lineterminator="\n")
        writer.writerow(ROW_HEADERS)
        writer.writerows(record_rows(records))
        return buffer.getvalue()
    with io.open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(ROW_HEADERS)
        writer.writerows(record_rows(records))
    return None


def to_ndjson(records, path=None):
    """Write records as newline-delimited JSON (or return the text)."""
    lines = []
    for record in records:
        entry = {
            "id": record.id,
            "description": record.description,
            "length": len(record.seq),
            "gc": gc_content(record.seq),
            "sequence": record.seq,
        }
        if hasattr(record, "qual") and record.qual:
            entry["quality"] = record.qual
        lines.append(json.dumps(entry, ensure_ascii=False))
    text = "\n".join(lines) + ("\n" if lines else "")
    if path is None:
        return text
    with io.open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write(text)
    return None


def convert(records, target, path=None, source_format="fasta"):
    """Convert records to ``target`` format.

    ``target`` is one of ``fasta``, ``fastq``, ``csv``, ``tsv``,
    ``ndjson``. Returns output text when *path* is None.
    """
    target = target.lower()
    if target == "fasta":
        if source_format == "fastq":
            records = fastq_to_fasta(records)
        return write_fasta(records, path)
    if target == "fastq":
        if source_format == "fasta":
            records = fasta_to_fastq(records)
        return write_fastq(records, path)
    if target == "csv":
        return to_csv(records, path)
    if target == "tsv":
        return to_tsv(records, path)
    if target in ("ndjson", "jsonl"):
        return to_ndjson(records, path)
    raise ValueError(
        "unsupported target format %r (fasta, fastq, csv, tsv, ndjson)"
        % target
    )
