import unittest

from seqforge.seqio import FastaRecord
from seqforge.stats import (
    gc_content,
    gc_fraction,
    gc_skew,
    length_stats,
    sequence_summary,
    sliding_gc,
)


class LengthStatsTests(unittest.TestCase):
    def test_basic(self):
        stats = length_stats([4, 3, 2, 1])
        self.assertEqual(stats["count"], 4)
        self.assertEqual(stats["total"], 10)
        self.assertEqual(stats["min"], 1)
        self.assertEqual(stats["max"], 4)
        self.assertEqual(stats["mean"], 2.5)
        self.assertEqual(stats["median"], 2.5)
        self.assertEqual(stats["n50"], 3)
        self.assertEqual(stats["l50"], 2)

    def test_empty(self):
        stats = length_stats([])
        self.assertEqual(stats["count"], 0)
        self.assertEqual(stats["n50"], 0)

    def test_single(self):
        stats = length_stats([7])
        self.assertEqual(stats["n50"], 7)
        self.assertEqual(stats["l50"], 1)


class GCTests(unittest.TestCase):
    def test_gc_content(self):
        self.assertEqual(gc_content("ATGC"), 50.0)
        self.assertEqual(gc_content("GGCC"), 100.0)
        self.assertEqual(gc_content("ATAT"), 0.0)

    def test_gc_content_rna(self):
        self.assertEqual(gc_content("AUGC"), 50.0)

    def test_gc_content_empty(self):
        self.assertEqual(gc_content(""), 0.0)
        self.assertEqual(gc_content("NNNN"), 0.0)

    def test_gc_content_ignores_ambiguous_bases(self):
        # N is excluded from the denominator: 2 GC of 3 called bases
        self.assertEqual(gc_content("GCNAN"), 66.67)

    def test_gc_fraction(self):
        # "GCNAT": 2 G/C, 2 A/T, N excluded from the called count
        self.assertEqual(gc_fraction("GCNAT"), (2, 2, 4))

    def test_gc_skew(self):
        self.assertEqual(gc_skew("GGCC"), 0.0)
        self.assertEqual(gc_skew("GGG"), 1.0)
        self.assertEqual(gc_skew("CCC"), -1.0)
        self.assertEqual(gc_skew("ATAT"), 0.0)

    def test_sliding_gc(self):
        windows = list(sliding_gc("ATATATGC", window=4, step=4))
        self.assertEqual(windows, [(0, 0.0), (4, 50.0)])

    def test_sliding_gc_skips_partial_window(self):
        windows = list(sliding_gc("ATGCAA", window=4, step=4))
        self.assertEqual(windows, [(0, 50.0)])

    def test_sliding_gc_rejects_bad_window(self):
        with self.assertRaises(ValueError):
            list(sliding_gc("ATGC", window=0))


class SummaryTests(unittest.TestCase):
    def test_sequence_summary(self):
        records = [
            FastaRecord("s1", "", "ATGC"),
            FastaRecord("s2", "", "GGCCGG"),
        ]
        summary, per_sequence = sequence_summary(records)
        self.assertEqual(summary["count"], 2)
        self.assertEqual(summary["total"], 10)
        self.assertEqual(summary["max"], 6)
        self.assertEqual(summary["gc"], 80.0)
        self.assertEqual(len(per_sequence), 2)
        self.assertEqual(per_sequence[0]["gc"], 50.0)

    def test_summary_empty(self):
        summary, per_sequence = sequence_summary([])
        self.assertEqual(summary["count"], 0)
        self.assertEqual(per_sequence, [])


if __name__ == "__main__":
    unittest.main()
