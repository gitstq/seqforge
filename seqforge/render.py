"""Rendering helpers: ASCII tables and JSON output."""

import json

ALIGN_LEFT = "left"
ALIGN_RIGHT = "right"
ALIGN_CENTER = "center"


def render_table(headers, rows, aligns=None, max_width=48):
    """Render a bordered ASCII table.

    ``aligns`` is an optional list of ``left``/``right``/``center`` per
    column (default: all left). Cells wider than ``max_width`` are
    truncated with an ellipsis. Pure ASCII output keeps logs portable.
    """
    headers = [str(h) for h in headers]
    rows = [[("" if cell is None else str(cell)) for cell in row]
            for row in rows]
    if aligns is None:
        aligns = [ALIGN_LEFT] * len(headers)

    def truncate(value):
        value = value.replace("\n", " ").replace("\t", " ")
        if max_width and len(value) > max_width:
            return value[:max_width - 3] + "..."
        return value

    headers = [truncate(h) for h in headers]
    rows = [[truncate(cell) for cell in row] for row in rows]

    def display_width(value):
        # CJK characters occupy two terminal cells.
        width = 0
        for char in value:
            width += 2 if ord(char) > 0x2E7F else 1
        return width

    def pad(value, width, align):
        padding = width - display_width(value)
        if padding <= 0:
            return value
        if align == ALIGN_RIGHT:
            return " " * padding + value
        if align == ALIGN_CENTER:
            left = padding // 2
            return " " * left + value + " " * (padding - left)
        return value + " " * padding

    columns = len(headers)
    widths = [display_width(h) for h in headers]
    for row in rows:
        for index in range(columns):
            if index < len(row):
                widths[index] = max(widths[index], display_width(row[index]))

    def border(left, mid, right):
        return left + mid.join("-" * (w + 2) for w in widths) + right

    def align_of(index):
        if index < len(aligns):
            return aligns[index]
        return ALIGN_LEFT

    lines = [border("+", "+", "+")]
    lines.append(
        "| " + " | ".join(
            pad(headers[i], widths[i], ALIGN_CENTER)
            for i in range(columns)
        ) + " |"
    )
    lines.append(border("+", "+", "+"))
    for row in rows:
        cells = []
        for index in range(columns):
            value = row[index] if index < len(row) else ""
            cells.append(pad(value, widths[index], align_of(index)))
        lines.append("| " + " | ".join(cells) + " |")
    lines.append(border("+", "+", "+"))
    return "\n".join(lines)


def to_json(obj):
    """Pretty JSON with non-ASCII characters preserved."""
    return json.dumps(obj, indent=2, ensure_ascii=False)


def format_bases(count):
    """Format a base count with a thousands separator."""
    return "{:,}".format(count)
