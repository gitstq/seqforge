import unittest

from seqforge.quality import (
    detect_phred_offset,
    filter_fastq,
    mean_quality,
    n_fraction,
    passes_filters,
    per_base_quality,
    phred_scores,
    trim_record,
    trim_sliding,
)
from seqforge.seqio import FastqRecord


def make_record(seq="ACGTACGT", qual="IIIIIIII", id="r1"):
    return FastqRecord(id, "", seq, qual)


class PhredTests(unittest.TestCase):
    def test_offset_33(self):
        self.assertEqual(detect_phred_offset("IIII"), 33)
        self.assertEqual(detect_phred_offset("!"), 33)

    def test_offset_64(self):
        self.assertEqual(detect_phred_offset("hhhh"), 64)

    def test_offset_empty(self):
        self.assertEqual(detect_phred_offset(""), 33)

    def test_phred_scores(self):
        self.assertEqual(phred_scores("I"), [40])
        self.assertEqual(phred_scores("#"), [2])

    def test_mean_quality(self):
        self.assertEqual(mean_quality("II"), 40.0)
        self.assertEqual(mean_quality(""), 0.0)


class TrimTests(unittest.TestCase):
    def test_trim_sliding_keeps_good_region(self):
        scores = [40, 40, 2, 2, 40, 40]
        self.assertEqual(trim_sliding(scores, 20, 2), (0, 3))

    def test_trim_sliding_all_good(self):
        scores = [30, 30, 30]
        self.assertEqual(trim_sliding(scores, 20, 2), (0, 3))

    def test_trim_sliding_all_bad(self):
        self.assertEqual(trim_sliding([5, 5, 5], 20, 2), (0, 0))

    def test_trim_sliding_empty(self):
        self.assertEqual(trim_sliding([], 20, 2), (0, 0))

    def test_trim_sliding_short_read(self):
        self.assertEqual(trim_sliding([40], 20, 4), (0, 1))
        self.assertEqual(trim_sliding([5], 20, 4), (0, 0))

    def test_trim_record(self):
        record = make_record(seq="ACGTACGT", qual="II##IIII")
        trim_record(record, threshold=25, window=4)
        # window means: 21, 21, 21, 30.5, 40 -> good windows start at index 3
        # and 4, so the kept span is [3, 8)
        self.assertEqual(record.seq, "TACGT")
        # qual[3:8] keeps the '#' at index 3, so seq and qual stay aligned
        self.assertEqual(record.qual, "#IIII")
        self.assertEqual(len(record.seq), len(record.qual))


class FilterTests(unittest.TestCase):
    def test_n_fraction(self):
        self.assertEqual(n_fraction("ACGT"), 0.0)
        self.assertEqual(n_fraction("ACGN"), 0.25)
        self.assertEqual(n_fraction(""), 0.0)

    def test_passes_filters(self):
        record = make_record(seq="ACGT", qual="IIII")
        self.assertTrue(passes_filters(record, min_len=4))
        self.assertFalse(passes_filters(record, min_len=5))
        self.assertFalse(passes_filters(record, max_len=3))
        self.assertFalse(passes_filters(record, min_mean_q=41))
        # ACGT has no N at all, so any N fraction limit is satisfied
        self.assertTrue(passes_filters(record, max_n_frac=0.5))

    def test_passes_filters_rejects_ambiguous_reads(self):
        record = make_record(seq="ACGN", qual="IIII", id="r2")
        self.assertFalse(passes_filters(record, max_n_frac=0.1))
        self.assertTrue(passes_filters(record, max_n_frac=0.5))

    def test_filter_fastq(self):
        records = [
            make_record(seq="ACGT", qual="IIII", id="a"),
            make_record(seq="A", qual="I", id="b"),
        ]
        kept = list(filter_fastq(records, min_len=4))
        self.assertEqual([r.id for r in kept], ["a"])


class PerBaseTests(unittest.TestCase):
    def test_per_base_quality(self):
        records = [
            make_record(seq="AC", qual="II", id="a"),
            make_record(seq="GT", qual="55", id="b"),
        ]
        # '5' = 53 -> Q20; 'I' = 73 -> Q40
        self.assertEqual(per_base_quality(records), [30.0, 30.0])

    def test_per_base_quality_empty(self):
        self.assertEqual(per_base_quality([]), [])

    def test_per_base_quality_truncates_to_shortest(self):
        records = [
            make_record(seq="ACGT", qual="IIII", id="a"),
            make_record(seq="AC", qual="II", id="b"),
        ]
        self.assertEqual(per_base_quality(records), [40.0, 40.0])


if __name__ == "__main__":
    unittest.main()
