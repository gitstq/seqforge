import unittest

from seqforge.genetic import (
    GENETIC_CODES,
    get_table,
    reverse_transcribe,
    stop_codons,
    transcribe,
    translate,
    translate_frames,
)


class TranslateTests(unittest.TestCase):
    def test_simple(self):
        self.assertEqual(translate("ATGGCC"), "MA")

    def test_stop_trimmed_by_default(self):
        self.assertEqual(translate("ATGTAA"), "M")

    def test_stop_kept_when_requested(self):
        self.assertEqual(translate("ATGTAA", trim_stop=False), "M*")

    def test_to_stop_halts(self):
        self.assertEqual(translate("ATGTAATAA"), "M*")
        self.assertEqual(translate("ATGTAATAA", to_stop=True), "M*")

    def test_rna_input(self):
        self.assertEqual(translate("AUGGCC"), "MA")

    def test_lowercase_input(self):
        self.assertEqual(translate("atggcc"), "MA")

    def test_ambiguous_codon_is_x(self):
        self.assertEqual(translate("ATNGGN"), "MX")

    def test_partial_codon_ignored(self):
        self.assertEqual(translate("ATGGC"), "M")

    def test_vertebrate_mito_table(self):
        self.assertEqual(translate("ATA", table=2), "M")
        self.assertEqual(translate("AGA", table=2), "")

    def test_yeast_mito_table(self):
        self.assertEqual(translate("CTG", table=3), "T")
        self.assertEqual(translate("TGA", table=3), "W")

    def test_ciliate_table(self):
        self.assertEqual(translate("TAA", table=6), "Q")

    def test_unsupported_table_raises(self):
        with self.assertRaises(ValueError):
            get_table(99)

    def test_stop_codons_standard(self):
        self.assertEqual(stop_codons(1), ["TAA", "TAG", "TGA"])


class FrameTests(unittest.TestCase):
    def test_transcribe(self):
        self.assertEqual(transcribe("ATGCT"), "AUGCU")

    def test_reverse_transcribe(self):
        self.assertEqual(reverse_transcribe("AUGCU"), "ATGCT")

    def test_three_frames(self):
        frames = translate_frames("AATGGCCTAA")
        labels = [label for label, _ in frames]
        self.assertEqual(labels, ["+1", "+2", "+3"])

    def test_six_frames(self):
        frames = translate_frames("AATGGCCTAA", both_strands=True)
        self.assertEqual(len(frames), 6)
        self.assertEqual(frames[3][0], "-1")

    def test_genetic_codes_registry(self):
        self.assertIn(1, GENETIC_CODES)
        self.assertEqual(GENETIC_CODES[1], "Standard")


if __name__ == "__main__":
    unittest.main()
