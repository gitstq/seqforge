"""FASTA / FASTQ reading, writing, detection and validation.

Only the Python standard library is used. Files are streamed with
generators where practical so that multi-gigabyte genomes stay affordable.
``.gz`` inputs are transparently decompressed.
"""

import gzip
import io
import os

from .alphabet import DNA_IUPAC, PROTEIN_IUPAC


class SeqForgeError(Exception):
    """Base class for all SeqForge errors."""


class SequenceFormatError(SeqForgeError):
    """Raised when a file cannot be parsed as FASTA/FASTQ."""


class FastaRecord(object):
    """A single FASTA record (id, description, sequence)."""

    __slots__ = ("id", "description", "seq")

    def __init__(self, id, description="", seq=""):
        self.id = id
        self.description = description or ""
        self.seq = seq or ""

    def header(self):
        if self.description:
            return "%s %s" % (self.id, self.description)
        return self.id

    def __len__(self):
        return len(self.seq)

    def __repr__(self):
        return "FastaRecord(id=%r, len=%d)" % (self.id, len(self.seq))


class FastqRecord(object):
    """A single FASTQ record (id, description, sequence, quality)."""

    __slots__ = ("id", "description", "seq", "qual")

    def __init__(self, id, description="", seq="", qual=""):
        self.id = id
        self.description = description or ""
        self.seq = seq or ""
        self.qual = qual or ""

    def header(self):
        if self.description:
            return "%s %s" % (self.id, self.description)
        return self.id

    def __len__(self):
        return len(self.seq)

    def __repr__(self):
        return "FastqRecord(id=%r, len=%d)" % (self.id, len(self.seq))


def detect_encoding(path, sample_size=65536):
    """Guess the text encoding of *path*.

    Order matters: ``utf-8`` is tried first, then the Chinese ``gb18030``
    superset (common for files produced on Windows), then ``latin-1`` which
    never fails and keeps bytes reversible.
    """
    with open(path, "rb") as handle:
        raw = handle.read(sample_size)
    for encoding in ("utf-8-sig", "utf-8", "gb18030", "latin-1"):
        try:
            raw.decode(encoding)
            return encoding
        except (UnicodeDecodeError, LookupError):
            continue
    return "latin-1"


def open_text(path, encoding=None):
    """Open *path* for text reading, transparently handling ``.gz``.

    Returns a context manager (use with ``with``).
    """
    if encoding is None:
        encoding = detect_encoding(path)
    if str(path).lower().endswith(".gz"):
        return io.TextIOWrapper(
            gzip.open(path, "rb"), encoding=encoding, newline=None
        )
    return io.open(path, "r", encoding=encoding, newline=None)


def _split_header(line):
    """Split a FASTA/FASTQ header into (id, description)."""
    header = line[1:].strip()
    if not header:
        return "", ""
    parts = header.split(None, 1)
    if len(parts) == 1:
        return parts[0], ""
    return parts[0], parts[1].strip()


def iter_fasta_text(text):
    """Yield :class:`FastaRecord` objects from a FASTA *text* blob."""
    record_id = None
    description = ""
    chunks = []
    started = False
    for lineno, raw in enumerate(text.splitlines(), start=1):
        line = raw.strip()
        if not line:
            continue
        if line.startswith(";"):  # legacy comment line
            continue
        if line.startswith(">"):
            if started:
                yield FastaRecord(record_id, description, "".join(chunks))
            record_id, description = _split_header(line)
            chunks = []
            started = True
            if not record_id:
                raise SequenceFormatError(
                    "line %d: FASTA header is missing an identifier" % lineno
                )
        elif started:
            chunks.append("".join(line.split()))
        else:
            raise SequenceFormatError(
                "line %d: FASTA data before the first '>' header" % lineno
            )
    if started:
        yield FastaRecord(record_id, description, "".join(chunks))


def read_fasta_text(text):
    """Return a list of :class:`FastaRecord` from a FASTA *text* blob."""
    return list(iter_fasta_text(text))


def iter_fasta(path, encoding=None):
    """Yield :class:`FastaRecord` objects from a FASTA file."""
    with open_text(path, encoding) as handle:
        text = handle.read()
    for record in iter_fasta_text(text):
        yield record


def read_fasta(path, encoding=None):
    """Return a list of :class:`FastaRecord` from a FASTA file."""
    return list(iter_fasta(path, encoding))


def iter_fastq_text(text):
    """Yield :class:`FastqRecord` objects from a FASTQ *text* blob.

    Only the canonical 4-line layout is accepted; wrapped/multi-line
    FASTQ is rejected with a clear error instead of being mis-parsed.
    """
    lines = text.splitlines()
    index = 0
    total = len(lines)
    while index < total:
        # Skip blank padding between records.
        while index < total and not lines[index].strip():
            index += 1
        if index >= total:
            break
        lineno = index + 1
        header = lines[index]
        if not header.startswith("@"):
            raise SequenceFormatError(
                "line %d: expected '@' FASTQ header, got %r" % (lineno, header[:40])
            )
        if index + 3 >= total:
            raise SequenceFormatError(
                "line %d: truncated FASTQ record (need 4 lines)" % lineno
            )
        seq = lines[index + 1].strip()
        plus = lines[index + 2].strip()
        qual = lines[index + 3].strip()
        if not plus.startswith("+"):
            raise SequenceFormatError(
                "line %d: expected '+' separator, got %r"
                % (index + 3, plus[:40])
            )
        if len(seq) != len(qual):
            raise SequenceFormatError(
                "line %d: sequence length %d != quality length %d"
                % (lineno, len(seq), len(qual))
            )
        record_id, description = _split_header(header)
        if not record_id:
            raise SequenceFormatError(
                "line %d: FASTQ header is missing an identifier" % lineno
            )
        yield FastqRecord(record_id, description, seq, qual)
        index += 4


def read_fastq_text(text):
    """Return a list of :class:`FastqRecord` from a FASTQ *text* blob."""
    return list(iter_fastq_text(text))


def iter_fastq(path, encoding=None):
    """Yield :class:`FastqRecord` objects from a FASTQ file."""
    with open_text(path, encoding) as handle:
        text = handle.read()
    for record in iter_fastq_text(text):
        yield record


def read_fastq(path, encoding=None):
    """Return a list of :class:`FastqRecord` from a FASTQ file."""
    return list(iter_fastq(path, encoding))


def detect_format(path, encoding=None):
    """Return ``"fasta"``, ``"fastq"`` or ``"unknown"`` for *path*."""
    with open_text(path, encoding) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                return "fasta"
            if line.startswith("@"):
                return "fastq"
            return "unknown"
    return "unknown"


def read_any(path, encoding=None):
    """Read FASTA or FASTQ, returning ``(format, records)``."""
    fmt = detect_format(path, encoding)
    if fmt == "fasta":
        return fmt, read_fasta(path, encoding)
    if fmt == "fastq":
        return fmt, read_fastq(path, encoding)
    raise SequenceFormatError(
        "%s: unrecognized sequence format (first non-empty line must "
        "start with '>' for FASTA or '@' for FASTQ)" % path
    )


def format_record(record, wrap=60):
    """Render a FASTA/FASTQ record as text (FASTA layout)."""
    lines = [">" + record.header()]
    seq = record.seq
    if wrap and wrap > 0:
        for start in range(0, len(seq), wrap):
            lines.append(seq[start:start + wrap])
    else:
        lines.append(seq)
    return "\n".join(lines)


def write_fasta(records, path=None, wrap=60, encoding="utf-8"):
    """Write records in FASTA format to *path* (or return text when None)."""
    if path is None:
        return "\n".join(format_record(r, wrap) for r in records) + "\n"
    with io.open(path, "w", encoding=encoding, newline="\n") as handle:
        for record in records:
            handle.write(format_record(record, wrap))
            handle.write("\n")


def format_fastq_record(record):
    """Render a FASTQ record as its canonical 4-line text block."""
    return "@%s\n%s\n+\n%s\n" % (record.header(), record.seq, record.qual)


def write_fastq(records, path=None, encoding="utf-8"):
    """Write records in FASTQ format to *path* (or return text when None)."""
    if path is None:
        return "".join(format_fastq_record(r) for r in records)
    with io.open(path, "w", encoding=encoding, newline="\n") as handle:
        for record in records:
            handle.write(format_fastq_record(record))


def _invalid_chars(seq, alphabet):
    return sorted({ch for ch in seq.upper() if ch not in alphabet and ch not in "-. "})


def validate_file(path, encoding=None):
    """Validate a FASTA/FASTQ file and return a report dictionary.

    The report contains counts plus a list of issues; it never raises for
    content problems (only for unreadable files).
    """
    fmt = detect_format(path, encoding)
    if fmt == "unknown":
        raise SequenceFormatError(
            "%s: unrecognized sequence format" % path
        )
    issues = []
    seen = {}
    total_records = 0
    total_bases = 0
    if fmt == "fasta":
        records = iter_fasta(path, encoding)
        for record in records:
            total_records += 1
            total_bases += len(record.seq)
            seen[record.id] = seen.get(record.id, 0) + 1
            if not record.seq:
                issues.append(
                    {"id": record.id, "severity": "error",
                     "message": "empty sequence"}
                )
            bad = _invalid_chars(record.seq, DNA_IUPAC)
            if bad:
                issues.append(
                    {"id": record.id, "severity": "error",
                     "message": "invalid characters: %s" % ", ".join(bad)}
                )
    else:
        records = iter_fastq(path, encoding)
        for record in records:
            total_records += 1
            total_bases += len(record.seq)
            seen[record.id] = seen.get(record.id, 0) + 1
            if not record.seq:
                issues.append(
                    {"id": record.id, "severity": "error",
                     "message": "empty sequence"}
                )
            bad = _invalid_chars(record.seq, DNA_IUPAC)
            if bad:
                issues.append(
                    {"id": record.id, "severity": "error",
                     "message": "invalid characters: %s" % ", ".join(bad)}
                )
    duplicates = sorted([rid for rid, n in seen.items() if n > 1])
    for rid in duplicates:
        issues.append(
            {"id": rid, "severity": "warning",
             "message": "duplicate identifier appears %d times" % seen[rid]}
        )
    if total_records == 0:
        issues.append({"id": "-", "severity": "error",
                       "message": "file contains no records"})
    return {
        "file": os.path.basename(str(path)),
        "format": fmt,
        "encoding": encoding or detect_encoding(path),
        "records": total_records,
        "total_bases": total_bases,
        "duplicate_ids": duplicates,
        "issues": issues,
        "ok": not [i for i in issues if i["severity"] == "error"],
    }


# Protein alphabet is re-exported for callers that validate translations.
__all__ = [
    "SeqForgeError", "SequenceFormatError", "FastaRecord", "FastqRecord",
    "detect_encoding", "open_text", "iter_fasta", "read_fasta",
    "iter_fastq", "read_fastq", "detect_format", "read_any",
    "format_record", "write_fasta", "format_fastq_record", "write_fastq",
    "validate_file", "PROTEIN_IUPAC",
]
