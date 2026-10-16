import unittest

from seqforge.motif import (
    find_motif,
    gc_percent,
    motif_to_regex,
    primer_stats,
    restriction_enzymes,
    restriction_sites,
    tm_salt_adjusted,
    tm_wallace,
)


class MotifTests(unittest.TestCase):
    def test_exact_matches(self):
        self.assertEqual(
            find_motif("ACGTACGT", "ACG"),
            [(0, 3, "+"), (4, 7, "+")],
        )

    def test_both_strands(self):
        matches = find_motif("ACGTACGT", "ACG", both_strands=True)
        self.assertEqual(matches, [(0, 3, "+"), (1, 4, "-"),
                                   (4, 7, "+"), (5, 8, "-")])

    def test_iupac_motif(self):
        # "ATGCGAT": positions 3-5 spell "CGA", which matches CGW and CGN
        self.assertEqual(find_motif("ATGCGAT", "CGW"), [(3, 6, "+")])
        self.assertEqual(find_motif("ATGCGAT", "CGN"), [(3, 6, "+")])

    def test_case_insensitive(self):
        self.assertEqual(find_motif("acgtacgt", "ACG"),
                         [(0, 3, "+"), (4, 7, "+")])

    def test_no_match(self):
        self.assertEqual(find_motif("ATGC", "TTTT"), [])

    def test_invalid_motif_raises(self):
        with self.assertRaises(ValueError):
            find_motif("ATGC", "AXG")

    def test_motif_to_regex_flags(self):
        pattern = motif_to_regex("acgt")
        self.assertTrue(pattern.search("ACGT"))


class EnzymeTests(unittest.TestCase):
    def test_builtin_table(self):
        table = restriction_enzymes()
        self.assertEqual(table["EcoRI"], "GAATTC")
        self.assertIn("BamHI", table)

    def test_find_sites(self):
        sites = restriction_sites("GAATTCGAATTC", ["EcoRI"])
        self.assertEqual(sites, [("EcoRI", 0, 6), ("EcoRI", 6, 12)])

    def test_find_sites_all_enzymes(self):
        sites = restriction_sites("GGATCC", ["BamHI"])
        self.assertEqual(len(sites), 1)

    def test_unknown_enzyme_raises(self):
        with self.assertRaises(ValueError):
            restriction_sites("GAATTC", ["NotAnEnzyme"])

    def test_lowercase_sequence(self):
        sites = restriction_sites("gaattc", ["EcoRI"])
        self.assertEqual(len(sites), 1)


class PrimerTests(unittest.TestCase):
    def test_gc_percent(self):
        self.assertEqual(gc_percent("ATGC"), 50.0)
        self.assertEqual(gc_percent(""), 0.0)

    def test_tm_wallace(self):
        self.assertEqual(tm_wallace("ATGCGT"), 18)

    def test_tm_salt_adjusted(self):
        self.assertAlmostEqual(
            tm_salt_adjusted("ATGCGTATGCGTATGCGT", na=0.05), 42.9, places=1)

    def test_tm_salt_adjusted_requires_positive_na(self):
        with self.assertRaises(ValueError):
            tm_salt_adjusted("ATGC", na=0)

    def test_primer_stats(self):
        stats = primer_stats("ATGCGT", na=0.05)
        self.assertEqual(stats["length"], 6)
        self.assertEqual(stats["gc"], 50.0)
        self.assertEqual(stats["tm_wallace"], 18)
        self.assertIn("tm_adjusted", stats)


if __name__ == "__main__":
    unittest.main()
