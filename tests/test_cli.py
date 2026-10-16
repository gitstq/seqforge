import io
import json
import os
import shutil
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout

from seqforge.cli import main

FASTA = ">s1 demo sequence\nATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG\n>s2\nGGCCGG\n"
FASTQ = "@r1\nACGTACGT\n+\nIIIIIIII\n@r2\nTTTTAAAA\n+r2\n########\n"


class CliTestCase(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.fasta_path = self._write("demo.fasta", FASTA)
        self.fastq_path = self._write("demo.fastq", FASTQ)

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _write(self, name, text):
        path = os.path.join(self.tmp, name)
        with io.open(path, "w", encoding="utf-8", newline="\n") as handle:
            handle.write(text)
        return path

    def _run(self, argv):
        """Run the CLI, returning (exit_code, stdout, stderr)."""
        out, err = io.StringIO(), io.StringIO()
        with redirect_stdout(out), redirect_stderr(err):
            code = main(argv)
        return code, out.getvalue(), err.getvalue()

    # -- stats -------------------------------------------------------------

    def test_stats_table(self):
        code, out, _ = self._run(["stats", self.fasta_path])
        self.assertEqual(code, 0)
        self.assertIn("records", out)
        self.assertIn("N50", out)

    def test_stats_json(self):
        code, out, _ = self._run(["stats", self.fasta_path, "--json"])
        self.assertEqual(code, 0)
        data = json.loads(out)
        self.assertEqual(data["summary"]["records"], 2)
        self.assertEqual(data["format"], "fasta")

    def test_stats_per_seq(self):
        code, out, _ = self._run(["stats", self.fasta_path, "--per-seq"])
        self.assertEqual(code, 0)
        self.assertIn("s1", out)

    # -- validate ----------------------------------------------------------

    def test_validate_ok(self):
        code, out, _ = self._run(["validate", self.fasta_path])
        self.assertEqual(code, 0)
        self.assertIn("OK", out)

    def test_validate_errors_exit_code(self):
        path = self._write("bad.fasta", ">s1\n\n")
        code, _, _ = self._run(["validate", path])
        self.assertEqual(code, 2)

    def test_validate_json(self):
        code, out, _ = self._run(["validate", self.fasta_path, "--json"])
        self.assertEqual(code, 0)
        self.assertTrue(json.loads(out)["ok"])

    # -- gc -----------------------------------------------------------------

    def test_gc(self):
        code, out, _ = self._run(["gc", self.fasta_path])
        self.assertEqual(code, 0)
        self.assertIn("gc (%)", out)

    def test_gc_window(self):
        code, out, _ = self._run(["gc", self.fasta_path, "--window", "10"])
        self.assertEqual(code, 0)
        self.assertIn("Sliding GC", out)

    # -- revcomp / transcribe ----------------------------------------------

    def test_revcomp(self):
        code, out, _ = self._run(["revcomp", self.fasta_path])
        self.assertEqual(code, 0)
        self.assertIn(">s1", out)
        # revcomp of the 39 bp demo sequence starts with "CTATC..."
        self.assertIn("CTATC", out)

    def test_transcribe(self):
        code, out, _ = self._run(["transcribe", self.fasta_path])
        self.assertEqual(code, 0)
        self.assertIn("AUGGCC", out)

    # -- translate -----------------------------------------------------------

    def test_translate_single_frame(self):
        code, out, _ = self._run(["translate", self.fasta_path])
        self.assertEqual(code, 0)
        self.assertIn(">s1", out)

    def test_translate_six_frames(self):
        code, out, _ = self._run(["translate", self.fasta_path, "--frames", "6"])
        self.assertEqual(code, 0)
        self.assertIn("frame -1", out)

    def test_translate_table_choice(self):
        code, out, _ = self._run(["translate", self.fasta_path, "--table", "2"])
        self.assertEqual(code, 0)

    # -- orf -----------------------------------------------------------------

    def test_orf(self):
        code, out, _ = self._run(["orf", self.fasta_path, "--min-len", "30"])
        self.assertEqual(code, 0)
        self.assertIn("strand", out)

    def test_orf_json(self):
        code, out, _ = self._run(["orf", self.fasta_path, "--json"])
        self.assertEqual(code, 0)
        data = json.loads(out)
        self.assertIn("orfs", data)

    def test_orf_none_found(self):
        path = self._write("short.fasta", ">x\nCCCCCCCCCCCCCC\n")
        code, out, _ = self._run(["orf", path])
        self.assertEqual(code, 0)
        self.assertIn("No ORFs", out)

    # -- find ------------------------------------------------------------------

    def test_find_motif(self):
        code, out, _ = self._run(["find", self.fasta_path, "TATAAA"])
        # TATAAA is a TATA box, not present in the demo sequence
        self.assertEqual(code, 0)
        self.assertIn("not found", out)

    def test_find_motif_hit(self):
        code, out, _ = self._run(["find", self.fasta_path, "GGCC"])
        self.assertEqual(code, 0)
        self.assertIn("strand", out)

    def test_find_invalid_motif(self):
        code, _, err = self._run(["find", self.fasta_path, "AXG"])
        self.assertEqual(code, 1)
        self.assertIn("error", err)

    # -- enzymes -----------------------------------------------------------------

    def test_enzymes_list(self):
        code, out, _ = self._run(["enzymes", "--list"])
        self.assertEqual(code, 0)
        self.assertIn("EcoRI", out)

    def test_enzymes_scan(self):
        code, out, _ = self._run(["enzymes", self.fasta_path,
                                  "--enzyme", "EcoRI"])
        self.assertEqual(code, 0)

    # -- tm ----------------------------------------------------------------------

    def test_tm(self):
        code, out, _ = self._run(["tm", "ATGCGTATGCGTATGCGT"])
        self.assertEqual(code, 0)
        self.assertIn("Tm Wallace", out)

    # -- fqstat / fqtrim -------------------------------------------------------------

    def test_fqstat(self):
        code, out, _ = self._run(["fqstat", self.fastq_path])
        self.assertEqual(code, 0)
        self.assertIn("phred offset", out)

    def test_fqstat_rejects_fasta(self):
        code, _, err = self._run(["fqstat", self.fasta_path])
        self.assertEqual(code, 1)
        self.assertIn("error", err)

    def test_fqtrim_output_file(self):
        out_path = os.path.join(self.tmp, "trimmed.fastq")
        code, _, err = self._run(["fqtrim", self.fastq_path,
                                  "-o", out_path, "--quality", "10"])
        self.assertEqual(code, 0)
        self.assertTrue(os.path.exists(out_path))
        with io.open(out_path, encoding="utf-8") as handle:
            self.assertTrue(handle.read().startswith("@r1"))
        self.assertIn("fqtrim:", err)

    # -- convert ---------------------------------------------------------------------

    def test_convert_fasta_to_csv(self):
        code, out, _ = self._run(["convert", self.fasta_path, "--to", "csv"])
        self.assertEqual(code, 0)
        self.assertIn("id,description,length,gc,sequence", out)

    def test_convert_fasta_to_fastq(self):
        code, out, _ = self._run(["convert", self.fasta_path, "--to", "fastq"])
        self.assertEqual(code, 0)
        self.assertTrue(out.startswith("@s1 demo sequence"))

    def test_convert_fastq_to_fasta(self):
        code, out, _ = self._run(["convert", self.fastq_path, "--to", "fasta"])
        self.assertEqual(code, 0)
        self.assertTrue(out.startswith(">r1"))

    def test_convert_ndjson(self):
        code, out, _ = self._run(["convert", self.fasta_path, "--to", "ndjson"])
        self.assertEqual(code, 0)
        first = json.loads(out.splitlines()[0])
        self.assertEqual(first["id"], "s1")

    # -- sample / dedup / info ----------------------------------------------------------

    def test_sample_reproducible(self):
        code, out, _ = self._run(["sample", self.fasta_path, "-n", "1",
                                  "--seed", "7"])
        self.assertEqual(code, 0)
        self.assertEqual(out.count(">"), 1)

    def test_dedup(self):
        path = self._write("dup.fasta", ">a\nATGC\n>b\nATGC\n>c\nGGCC\n")
        code, out, _ = self._run(["dedup", path])
        self.assertEqual(code, 0)
        self.assertEqual(out.count(">"), 2)

    def test_info(self):
        code, out, _ = self._run(["info", self.fasta_path])
        self.assertEqual(code, 0)
        self.assertIn("format=fasta", out)

    # -- error handling ------------------------------------------------------------------

    def test_no_command_prints_help(self):
        code, out, _ = self._run([])
        self.assertEqual(code, 1)
        self.assertIn("usage", out.lower())

    def test_missing_file(self):
        missing = os.path.join(self.tmp, "nope.fasta")
        code, _, err = self._run(["stats", missing])
        self.assertEqual(code, 1)
        self.assertIn("error", err)

    def test_output_to_file(self):
        out_path = os.path.join(self.tmp, "rc.fasta")
        code, _, _ = self._run(["revcomp", self.fasta_path, "-o", out_path])
        self.assertEqual(code, 0)
        self.assertTrue(os.path.exists(out_path))


if __name__ == "__main__":
    unittest.main()
