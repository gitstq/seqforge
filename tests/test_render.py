import json
import unittest

from seqforge.render import format_bases, render_table, to_json


class RenderTableTests(unittest.TestCase):
    def test_basic_structure(self):
        table = render_table(["name", "value"], [["a", "1"]])
        lines = table.splitlines()
        self.assertEqual(lines[0], "+-------+-------+")
        self.assertIn("| name | value |", lines[1])
        self.assertIn("| a    | 1     |", lines[3])
        self.assertEqual(lines[2], lines[4])

    def test_right_align(self):
        table = render_table(["n"], [["12345"]], aligns=["right"])
        self.assertIn("| 12345 |", table)

    def test_center_align(self):
        table = render_table(["n"], [["x"]], aligns=["center"])
        self.assertIn("| x |", table)

    def test_truncation(self):
        table = render_table(["s"], [["A" * 60]], max_width=10)
        self.assertIn("...", table)

    def test_cjk_width(self):
        table = render_table(["序列"], [["中文内容"]])
        # CJK cells count as two terminal columns; the header is centred.
        self.assertIn("|   序列   |", table.splitlines()[1])
        self.assertIn("| 中文内容 |", table.splitlines()[3])

    def test_none_cell(self):
        table = render_table(["a", "b"], [[None, "x"]])
        self.assertIn("|   |", table)

    def test_empty_rows(self):
        table = render_table(["a"], [])
        self.assertIn("| a |", table)


class JsonTests(unittest.TestCase):
    def test_unicode_preserved(self):
        text = to_json({"key": "中文"})
        self.assertIn("中文", text)
        self.assertEqual(json.loads(text), {"key": "中文"})

    def test_indent(self):
        text = to_json({"a": 1})
        self.assertIn("\n", text)


class FormatBasesTests(unittest.TestCase):
    def test_thousands_separator(self):
        self.assertEqual(format_bases(1234567), "1,234,567")
        self.assertEqual(format_bases(0), "0")


if __name__ == "__main__":
    unittest.main()
