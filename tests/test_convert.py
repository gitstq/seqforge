import csv
import io
import json
import unittest

from seqforge.convert import (
    convert,
    fasta_to_fastq,
    fastq_to_fasta,
    to_csv,
    to_ndjson,
    to_tsv,
)
from seqforge.seqio import read_fasta_text, read_fastq_text


class ConversionTests(unittest.TestCase):
    def setUp(self):
        self.fasta = read_fasta_text(">s1 d1\nATGGCC\n>s2\nGGCCGG\n")
        self.fastq = read_fastq_text("@r1\nACGT\n+\nIIII\n")

    def test_fasta_to_fastq(self):
        records = fasta_to_fastq(self.fasta)
        self.assertEqual(len(records), 2)
        self.assertEqual(records[0].qual, "I" * 6)
        self.assertEqual(records[0].id, "s1")

    def test_fastq_to_fasta(self):
        records = fastq_to_fasta(self.fastq)
        self.assertEqual(records[0].seq, "ACGT")
        self.assertFalse(hasattr(records[0], "qual"))

    def test_csv_text(self):
        text = to_csv(self.fasta)
        rows = list(csv.reader(io.StringIO(text)))
        self.assertEqual(rows[0], ["id", "description", "length", "gc",
                                   "sequence"])
        self.assertEqual(rows[1][0], "s1")
        self.assertEqual(rows[1][1], "d1")
        self.assertEqual(rows[1][2], "6")
        self.assertEqual(rows[1][3], "66.67")

    def test_tsv_text(self):
        text = to_tsv(self.fasta)
        rows = list(csv.reader(io.StringIO(text), delimiter="\t"))
        self.assertEqual(rows[1][0], "s1")

    def test_ndjson_text(self):
        text = to_ndjson(self.fasta)
        lines = [line for line in text.splitlines() if line]
        self.assertEqual(len(lines), 2)
        first = json.loads(lines[0])
        self.assertEqual(first["id"], "s1")
        self.assertEqual(first["length"], 6)

    def test_ndjson_keeps_quality(self):
        text = to_ndjson(self.fastq)
        first = json.loads(text.splitlines()[0])
        self.assertEqual(first["quality"], "IIII")

    def test_convert_fasta_to_fastq(self):
        text = convert(self.fasta, "fastq", None, source_format="fasta")
        self.assertTrue(text.startswith("@s1 d1"))

    def test_convert_fastq_to_fasta(self):
        text = convert(self.fastq, "fasta", None, source_format="fastq")
        self.assertTrue(text.startswith(">r1"))

    def test_convert_csv(self):
        text = convert(self.fasta, "csv", None, source_format="fasta")
        self.assertIn("s1", text)

    def test_convert_bad_target(self):
        with self.assertRaises(ValueError):
            convert(self.fasta, "xml", None, source_format="fasta")


if __name__ == "__main__":
    unittest.main()
