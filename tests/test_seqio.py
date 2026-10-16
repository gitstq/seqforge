import gzip
import io
import os
import shutil
import tempfile
import unittest

from seqforge.seqio import (
    FastaRecord,
    FastqRecord,
    SequenceFormatError,
    detect_format,
    read_any,
    read_fasta,
    read_fasta_text,
    read_fastq,
    read_fastq_text,
    validate_file,
    write_fasta,
    write_fastq,
)

FASTA_TEXT = (
    ">s1 first sequence\n"
    "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG\n"
    ">s2\n"
    "GGGCCT\n"
)

FASTQ_TEXT = (
    "@r1 read one\n"
    "ACGTACGT\n"
    "+\n"
    "IIIIIIII\n"
    "@r2\n"
    "TTTTAAAA\n"
    "+r2\n"
    "########\n"
)


class FastaParsingTests(unittest.TestCase):
    def test_parse_text(self):
        records = read_fasta_text(FASTA_TEXT)
        self.assertEqual(len(records), 2)
        self.assertEqual(records[0].id, "s1")
        self.assertEqual(records[0].description, "first sequence")
        self.assertEqual(records[1].id, "s2")
        self.assertEqual(records[1].description, "")
        self.assertEqual(len(records[0]), 39)

    def test_multi_line_sequence(self):
        records = read_fasta_text(">s\nATGC\nATGC\n")
        self.assertEqual(records[0].seq, "ATGCATGC")

    def test_blank_lines_tolerated(self):
        records = read_fasta_text("\n>s\n\nATGC\n\n")
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].seq, "ATGC")

    def test_data_before_header_raises(self):
        with self.assertRaises(SequenceFormatError):
            read_fasta_text("ATGC\n>s\nATGC\n")

    def test_roundtrip(self):
        records = read_fasta_text(FASTA_TEXT)
        text = write_fasta(records, None, wrap=10)
        again = read_fasta_text(text)
        self.assertEqual([r.seq for r in records], [r.seq for r in again])


class FastqParsingTests(unittest.TestCase):
    def test_parse_text(self):
        records = read_fastq_text(FASTQ_TEXT)
        self.assertEqual(len(records), 2)
        self.assertEqual(records[0].id, "r1")
        self.assertEqual(records[0].description, "read one")
        self.assertEqual(records[0].seq, "ACGTACGT")
        self.assertEqual(records[0].qual, "IIIIIIII")
        self.assertEqual(records[1].qual, "########")

    def test_length_mismatch_raises(self):
        with self.assertRaises(SequenceFormatError):
            read_fastq_text("@r\nACGT\n+\nII\n")

    def test_missing_plus_raises(self):
        with self.assertRaises(SequenceFormatError):
            read_fastq_text("@r\nACGT\n-\nIIII\n")

    def test_truncated_raises(self):
        with self.assertRaises(SequenceFormatError):
            read_fastq_text("@r\nACGT\n+\n")

    def test_roundtrip(self):
        records = read_fastq_text(FASTQ_TEXT)
        text = write_fastq(records, None)
        again = read_fastq_text(text)
        self.assertEqual([r.qual for r in records], [r.qual for r in again])


class FileTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _write(self, name, text, encoding="utf-8"):
        path = os.path.join(self.tmp, name)
        with io.open(path, "w", encoding=encoding, newline="\n") as handle:
            handle.write(text)
        return path

    def test_detect_format(self):
        fasta = self._write("a.fasta", FASTA_TEXT)
        fastq = self._write("a.fastq", FASTQ_TEXT)
        self.assertEqual(detect_format(fasta), "fasta")
        self.assertEqual(detect_format(fastq), "fastq")

    def test_read_fasta_file(self):
        path = self._write("a.fasta", FASTA_TEXT)
        records = read_fasta(path)
        self.assertEqual(len(records), 2)

    def test_read_fastq_file(self):
        path = self._write("a.fastq", FASTQ_TEXT)
        records = read_fastq(path)
        self.assertEqual(len(records), 2)

    def test_read_any(self):
        fasta = self._write("a.fasta", FASTA_TEXT)
        fastq = self._write("a.fastq", FASTQ_TEXT)
        self.assertEqual(read_any(fasta)[0], "fasta")
        self.assertEqual(read_any(fastq)[0], "fastq")

    def test_gzip_roundtrip(self):
        path = os.path.join(self.tmp, "a.fasta.gz")
        with gzip.open(path, "wb") as handle:
            handle.write(FASTA_TEXT.encode("utf-8"))
        records = read_fasta(path)
        self.assertEqual(len(records), 2)

    def test_validate_ok(self):
        path = self._write("a.fasta", FASTA_TEXT)
        report = validate_file(path)
        self.assertTrue(report["ok"])
        self.assertEqual(report["records"], 2)
        self.assertEqual(report["issues"], [])

    def test_validate_reports_issues(self):
        text = ">s1\n\n>s1\nATGC\n>s2\nATGZ\n"
        path = self._write("bad.fasta", text)
        report = validate_file(path)
        self.assertFalse(report["ok"])
        self.assertIn("s1", report["duplicate_ids"])
        messages = " ".join(i["message"] for i in report["issues"])
        self.assertIn("empty sequence", messages)
        self.assertIn("invalid characters", messages)

    def test_validate_unknown_format(self):
        path = self._write("weird.txt", "hello world\n")
        with self.assertRaises(SequenceFormatError):
            validate_file(path)

    def test_missing_file_raises(self):
        with self.assertRaises(IOError):
            read_fasta(os.path.join(self.tmp, "nope.fasta"))


if __name__ == "__main__":
    unittest.main()
