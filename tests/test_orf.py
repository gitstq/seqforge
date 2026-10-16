import unittest

from seqforge.orf import find_orfs, longest_orf


class OrfTests(unittest.TestCase):
    def test_simple_forward_orf(self):
        orfs = find_orfs("ATGAAATAA", min_len=6)
        self.assertEqual(len(orfs), 1)
        orf = orfs[0]
        self.assertEqual(orf["strand"], "+")
        self.assertEqual(orf["frame"], 1)
        self.assertEqual(orf["start"], 0)
        self.assertEqual(orf["end"], 9)
        self.assertEqual(orf["length"], 9)
        self.assertEqual(orf["protein"], "MK")

    def test_reverse_strand_orf(self):
        # revcomp("ATGAAATAA") = "TTATTTCAT"
        orfs = find_orfs("TTATTTCAT", min_len=6)
        self.assertEqual(len(orfs), 1)
        orf = orfs[0]
        self.assertEqual(orf["strand"], "-")
        self.assertEqual(orf["start"], 0)
        self.assertEqual(orf["end"], 9)
        self.assertEqual(orf["protein"], "MK")

    def test_min_len_filter(self):
        self.assertEqual(find_orfs("ATGAAATAA", min_len=30), [])

    def test_forward_only(self):
        orfs = find_orfs("TTATTTCAT", min_len=6, both_strands=False)
        self.assertEqual(orfs, [])

    def test_alternative_starts(self):
        orfs = find_orfs("ATGATGAAATAA", min_len=6)
        starts = [(o["start"], o["end"], o["length"]) for o in orfs]
        self.assertIn((0, 12, 12), starts)
        self.assertIn((3, 12, 9), starts)

    def test_no_orf_without_stop(self):
        self.assertEqual(find_orfs("ATGAAACCC", min_len=6), [])

    def test_orf_requires_start_codon(self):
        self.assertEqual(find_orfs("TAAAAAATAA", min_len=6), [])

    def test_orf_skipped_when_frame_misaligned(self):
        # ATG starts at index 1: frame 2
        orfs = find_orfs("AATGAAATAA", min_len=6)
        self.assertEqual(len(orfs), 1)
        self.assertEqual(orfs[0]["frame"], 2)
        self.assertEqual(orfs[0]["start"], 1)

    def test_longest_orf(self):
        orf = longest_orf("ATGATGAAATAA", min_len=6)
        self.assertIsNotNone(orf)
        self.assertEqual(orf["length"], 12)

    def test_longest_orf_none(self):
        self.assertIsNone(longest_orf("CCCCCC", min_len=6))


if __name__ == "__main__":
    unittest.main()
