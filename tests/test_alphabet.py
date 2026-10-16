import unittest

from seqforge.alphabet import (
    complement,
    is_nucleic,
    is_protein,
    is_rna,
    normalize,
    reverse_complement,
)


class ReverseComplementTests(unittest.TestCase):
    def test_simple_dna(self):
        self.assertEqual(reverse_complement("ATGC"), "GCAT")

    def test_case_preserved(self):
        self.assertEqual(reverse_complement("atGC"), "GCat")

    def test_rna_complement(self):
        self.assertEqual(reverse_complement("AUGC", rna=True), "GCAU")

    def test_iupac_ambiguity(self):
        # complement-then-reverse: RYKM -> YRMK -> KMRY
        self.assertEqual(reverse_complement("RYKM"), "KMRY")
        # BDHVN -> VHDBN -> NBDHV
        self.assertEqual(reverse_complement("BDHVN"), "NBDHV")

    def test_gaps_preserved(self):
        self.assertEqual(reverse_complement("AT-GC"), "GC-AT")

    def test_invalid_character_raises(self):
        with self.assertRaises(ValueError):
            reverse_complement("ATGX")

    def test_complement_without_reversing(self):
        self.assertEqual(complement("ATGC"), "TACG")
        self.assertEqual(complement("AUGC", rna=True), "UACG")


class AlphabetPredicateTests(unittest.TestCase):
    def test_is_nucleic(self):
        self.assertTrue(is_nucleic("ATGCRYSWKMBDHVN"))
        self.assertFalse(is_nucleic("ATGZ"))

    def test_is_rna(self):
        self.assertTrue(is_rna("AUGC"))
        self.assertFalse(is_rna("ATGC"))

    def test_is_protein(self):
        self.assertTrue(is_protein("MAKVL*"))
        self.assertFalse(is_protein("MAK1"))

    def test_normalize(self):
        self.assertEqual(normalize("aUgc"), "ATGC")


if __name__ == "__main__":
    unittest.main()
