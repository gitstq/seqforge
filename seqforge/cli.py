"""SeqForge command line interface."""

import argparse
import io
import random
import sys

from . import __version__
from .alphabet import normalize, reverse_complement
from .convert import convert as convert_records
from .genetic import GENETIC_CODES, transcribe, translate, translate_frames
from .motif import (
    find_motif,
    primer_stats,
    restriction_enzymes,
    restriction_sites,
)
from .orf import find_orfs
from .quality import (
    detect_phred_offset,
    mean_quality,
    passes_filters,
    per_base_quality,
    trim_record,
)
from .render import format_bases, render_table, to_json
from .seqio import (
    FastaRecord,
    SeqForgeError,
    detect_format,
    read_any,
    validate_file,
    write_fasta,
    write_fastq,
)
from .stats import gc_content, gc_skew, sequence_summary, sliding_gc

#: Exit code returned by ``validate`` when errors are found.
EXIT_VALIDATION_ERROR = 2


def _write(text, path):
    """Write *text* to *path* or stdout when path is None."""
    if path is None:
        sys.stdout.write(text)
        return
    with io.open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write(text)


def _as_fasta_records(records, sequences):
    """Pair records with new sequences, keeping ids/descriptions."""
    return [
        FastaRecord(record.id, record.description, sequence)
        for record, sequence in zip(records, sequences)
    ]


# ---------------------------------------------------------------------------
# command handlers
# ---------------------------------------------------------------------------

def cmd_stats(args):
    fmt, records = read_any(args.file, args.encoding)
    summary, per_sequence = sequence_summary(records)
    if args.json:
        payload = {"file": args.file, "format": fmt, "summary": summary}
        if args.per_seq:
            payload["sequences"] = per_sequence
        sys.stdout.write(to_json(payload) + "\n")
        return 0
    rows = [
        ["format", fmt],
        ["records", summary["count"]],
        ["total bases", format_bases(summary["total"])],
        ["min length", summary["min"]],
        ["max length", summary["max"]],
        ["mean length", summary["mean"]],
        ["median length", summary["median"]],
        ["length stdev", summary["stdev"]],
        ["N50", summary["n50"]],
        ["L50", summary["l50"]],
        ["GC content (%)", summary["gc"]],
    ]
    sys.stdout.write(
        render_table(["metric", "value"], rows, aligns=["left", "right"]) + "\n"
    )
    if args.per_seq:
        sys.stdout.write(render_table(
            ["id", "length", "gc (%)"],
            [[r["id"], r["length"], r["gc"]] for r in per_sequence],
            aligns=["left", "right", "right"],
        ) + "\n")
    return 0


def cmd_validate(args):
    report = validate_file(args.file, args.encoding)
    if args.json:
        sys.stdout.write(to_json(report) + "\n")
    else:
        rows = [
            ["file", report["file"]],
            ["format", report["format"]],
            ["encoding", report["encoding"]],
            ["records", report["records"]],
            ["total bases", format_bases(report["total_bases"])],
            ["duplicate ids", len(report["duplicate_ids"])],
            ["issues", len(report["issues"])],
            ["status", "OK" if report["ok"] else "FAILED"],
        ]
        sys.stdout.write(
            render_table(["field", "value"], rows, aligns=["left", "right"]) + "\n"
        )
        if report["issues"]:
            sys.stdout.write(render_table(
                ["id", "severity", "message"],
                [[i["id"], i["severity"], i["message"]] for i in report["issues"]],
            ) + "\n")
    return 0 if report["ok"] else EXIT_VALIDATION_ERROR


def cmd_gc(args):
    fmt, records = read_any(args.file, args.encoding)
    rows = []
    for record in records:
        rows.append([record.id, len(record.seq), gc_content(record.seq),
                     gc_skew(record.seq)])
    if args.json:
        payload = {
            "file": args.file,
            "sequences": [
                {"id": r[0], "length": r[1], "gc": r[2], "gc_skew": r[3]}
                for r in rows
            ],
        }
        if args.window:
            payload["windows"] = {}
            for record in records:
                windows = list(sliding_gc(record.seq, args.window,
                                          args.step or args.window))
                shown = windows[:args.max_windows]
                payload["windows"][record.id] = [
                    {"start": start + 1, "gc": gc} for start, gc in shown
                ]
                if len(windows) > len(shown):
                    payload["windows"][record.id].append(
                        {"start": "...", "gc": "%d more windows omitted"
                         % (len(windows) - len(shown))})
        sys.stdout.write(to_json(payload) + "\n")
        return 0
    sys.stdout.write(render_table(
        ["id", "length", "gc (%)", "gc skew"],
        rows, aligns=["left", "right", "right", "right"],
    ) + "\n")
    if args.window:
        for record in records:
            windows = list(sliding_gc(record.seq, args.window,
                                      args.step or args.window))
            shown = windows[:args.max_windows]
            sys.stdout.write("Sliding GC for %s (window=%d, step=%d):\n"
                             % (record.id, args.window,
                                args.step or args.window))
            sys.stdout.write(render_table(
                ["start (1-based)", "gc (%)"],
                [[start + 1, gc] for start, gc in shown],
                aligns=["right", "right"],
            ) + "\n")
            if len(windows) > len(shown):
                sys.stdout.write("... %d more windows omitted\n"
                                 % (len(windows) - len(shown)))
    return 0


def cmd_revcomp(args):
    fmt, records = read_any(args.file, args.encoding)
    sequences = [reverse_complement(r.seq, rna=args.rna) for r in records]
    _write(write_fasta(_as_fasta_records(records, sequences), None), args.output)
    return 0


def cmd_transcribe(args):
    fmt, records = read_any(args.file, args.encoding)
    sequences = [transcribe(r.seq) for r in records]
    _write(write_fasta(_as_fasta_records(records, sequences), None), args.output)
    return 0


def cmd_translate(args):
    fmt, records = read_any(args.file, args.encoding)
    output_records = []
    if args.frames == 1:
        for record in records:
            protein = translate(record.seq, table=args.table)
            output_records.append(
                FastaRecord(record.id, record.description, protein))
    else:
        both = args.frames == 6
        for record in records:
            for label, protein in translate_frames(
                    record.seq, table=args.table, both_strands=both):
                output_records.append(
                    FastaRecord("%s frame %s" % (record.id, label),
                                record.description, protein))
    _write(write_fasta(output_records, None), args.output)
    return 0


def cmd_orf(args):
    fmt, records = read_any(args.file, args.encoding)
    orfs = []
    for record in records:
        for orf in find_orfs(record.seq, min_len=args.min_len,
                             table=args.table,
                             both_strands=not args.forward_only):
            orf["id"] = record.id
            orfs.append(orf)
    if args.json:
        sys.stdout.write(to_json({
            "file": args.file,
            "min_len": args.min_len,
            "table": args.table,
            "count": len(orfs),
            "orfs": orfs,
        }) + "\n")
        return 0
    if not orfs:
        sys.stdout.write("No ORFs found (min length %d bp).\n" % args.min_len)
        return 0
    sys.stdout.write(render_table(
        ["id", "strand", "frame", "start", "end", "length", "protein"],
        [[o["id"], o["strand"], o["frame"], o["start"] + 1, o["end"],
          o["length"], o["protein"][:40]] for o in orfs],
        aligns=["left", "center", "right", "right", "right", "right", "left"],
    ) + "\n")
    return 0


def cmd_find(args):
    fmt, records = read_any(args.file, args.encoding)
    matches = []
    for record in records:
        for start, end, strand in find_motif(
                record.seq, args.motif, both_strands=args.both_strands):
            matches.append({"id": record.id, "start": start + 1,
                            "end": end, "strand": strand})
    if args.json:
        sys.stdout.write(to_json({
            "file": args.file, "motif": args.motif,
            "count": len(matches), "matches": matches,
        }) + "\n")
        return 0
    if not matches:
        sys.stdout.write("Motif %s not found.\n" % args.motif)
        return 0
    sys.stdout.write(render_table(
        ["id", "start", "end", "strand"],
        [[m["id"], m["start"], m["end"], m["strand"]] for m in matches],
        aligns=["left", "right", "right", "center"],
    ) + "\n")
    return 0


def cmd_enzymes(args):
    if args.list:
        table = restriction_enzymes()
        if args.json:
            sys.stdout.write(to_json(table) + "\n")
            return 0
        sys.stdout.write(render_table(
            ["enzyme", "recognition"],
            [[name, seq] for name, seq in sorted(table.items())],
            aligns=["left", "left"],
        ) + "\n")
        return 0
    if not args.file:
        raise SeqForgeError("enzymes: input file required (or use --list)")
    fmt, records = read_any(args.file, args.encoding)
    names = None
    if args.enzyme:
        names = [n.strip() for n in args.enzyme.split(",") if n.strip()]
    sites = []
    for record in records:
        for name, start, end in restriction_sites(record.seq, names):
            sites.append({"id": record.id, "enzyme": name,
                          "start": start + 1, "end": end})
    if args.json:
        sys.stdout.write(to_json({
            "file": args.file, "count": len(sites), "sites": sites,
        }) + "\n")
        return 0
    if not sites:
        sys.stdout.write("No recognition sites found.\n")
        return 0
    sys.stdout.write(render_table(
        ["id", "enzyme", "start", "end"],
        [[s["id"], s["enzyme"], s["start"], s["end"]] for s in sites],
        aligns=["left", "left", "right", "right"],
    ) + "\n")
    return 0


def cmd_tm(args):
    stats = [primer_stats(primer, na=args.na) for primer in args.primers]
    if args.json:
        sys.stdout.write(to_json(stats) + "\n")
        return 0
    sys.stdout.write(render_table(
        ["primer", "length", "gc (%)", "Tm Wallace (C)", "Tm adjusted (C)"],
        [[s["primer"], s["length"], s["gc"], s["tm_wallace"], s["tm_adjusted"]]
         for s in stats],
        aligns=["left", "right", "right", "right", "right"],
    ) + "\n")
    return 0


def cmd_fqstat(args):
    fmt, records = read_any(args.file, args.encoding)
    if fmt != "fastq":
        raise SeqForgeError("fqstat requires a FASTQ file (detected: %s)" % fmt)
    if not records:
        raise SeqForgeError("FASTQ file contains no records")
    offset = detect_phred_offset(records[0].qual)
    mean_qs = [mean_quality(r.qual) for r in records]
    per_base = per_base_quality(records)
    summary, _ = sequence_summary(records)
    shown = per_base[:args.max_positions]
    payload = {
        "file": args.file,
        "records": len(records),
        "phred_offset": offset,
        "mean_quality": round(sum(mean_qs) / len(mean_qs), 2),
        "min_mean_quality": min(mean_qs),
        "length": {k: summary[k] for k in
                   ("count", "min", "max", "mean", "median")},
        "per_base_quality": shown,
    }
    if args.json:
        sys.stdout.write(to_json(payload) + "\n")
        return 0
    rows = [
        ["records", len(records)],
        ["phred offset", offset],
        ["mean quality (Q)", payload["mean_quality"]],
        ["worst mean quality (Q)", payload["min_mean_quality"]],
        ["min length", summary["min"]],
        ["max length", summary["max"]],
        ["mean length", summary["mean"]],
        ["median length", summary["median"]],
    ]
    sys.stdout.write(
        render_table(["metric", "value"], rows, aligns=["left", "right"]) + "\n"
    )
    if shown:
        sys.stdout.write("Per-base mean quality (first %d positions):\n"
                         % len(shown))
        sys.stdout.write(render_table(
            ["position", "mean quality (Q)"],
            [[i + 1, q] for i, q in enumerate(shown)],
            aligns=["right", "right"],
        ) + "\n")
    return 0


def cmd_fqtrim(args):
    fmt, records = read_any(args.file, args.encoding)
    if fmt != "fastq":
        raise SeqForgeError("fqtrim requires a FASTQ file (detected: %s)" % fmt)
    total_in = len(records)
    bases_in = sum(len(r.seq) for r in records)
    kept = []
    for record in records:
        trim_record(record, threshold=args.quality, window=args.window)
        if passes_filters(record, min_len=args.min_len, max_len=args.max_len,
                          min_mean_q=args.min_mean_q, max_n_frac=args.max_n):
            kept.append(record)
    bases_out = sum(len(r.seq) for r in kept)
    _write(write_fastq(kept, None), args.output)
    sys.stderr.write(
        "fqtrim: kept %d/%d records, %d/%d bases (quality>=%d, window=%d, "
        "min_len=%d, max_len=%s, min_mean_q=%s, max_n=%s)\n"
        % (len(kept), total_in, bases_out, bases_in, args.quality,
           args.window, args.min_len,
           args.max_len if args.max_len is not None else "inf",
           args.min_mean_q, args.max_n)
    )
    return 0


def cmd_convert(args):
    fmt, records = read_any(args.file, args.encoding)
    text = convert_records(records, args.to, None, source_format=fmt)
    _write(text, args.output)
    return 0


def cmd_sample(args):
    fmt, records = read_any(args.file, args.encoding)
    rng = random.Random(args.seed)
    if args.number >= len(records):
        sampled = list(records)
    else:
        sampled = rng.sample(records, args.number)
    if fmt == "fastq":
        _write(write_fastq(sampled, None), args.output)
    else:
        _write(write_fasta(sampled, None), args.output)
    sys.stderr.write("sample: %d of %d records (seed=%s)\n"
                     % (len(sampled), len(records), args.seed))
    return 0


def cmd_dedup(args):
    fmt, records = read_any(args.file, args.encoding)
    seen = set()
    kept = []
    duplicates = 0
    for record in records:
        key = normalize(record.seq)
        if key in seen:
            duplicates += 1
            continue
        seen.add(key)
        kept.append(record)
    if fmt == "fastq":
        _write(write_fastq(kept, None), args.output)
    else:
        _write(write_fasta(kept, None), args.output)
    sys.stderr.write("dedup: %d unique, %d duplicates removed\n"
                     % (len(kept), duplicates))
    return 0


def cmd_info(args):
    fmt = detect_format(args.file, args.encoding)
    if fmt == "unknown":
        raise SeqForgeError("info: unrecognized sequence format")
    _, records = read_any(args.file, args.encoding)
    preview = records[:3]
    rows = [[r.id, r.description[:40], len(r.seq), gc_content(r.seq)]
            for r in preview]
    sys.stdout.write(render_table(
        ["id", "description", "length", "gc (%)"], rows,
        aligns=["left", "left", "right", "right"],
    ) + "\n")
    sys.stdout.write("format=%s records=%d (showing first %d)\n"
                     % (fmt, len(records), len(preview)))
    return 0


# ---------------------------------------------------------------------------
# parser
# ---------------------------------------------------------------------------

def build_parser():
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument(
        "--encoding", default=None,
        help="force input encoding (default: auto-detect utf-8/gb18030)")
    common.add_argument(
        "-o", "--output", default=None,
        help="write output to FILE instead of stdout")

    parser = argparse.ArgumentParser(
        prog="seqforge",
        description="SeqForge - zero-dependency bioinformatics sequence "
                    "toolkit (FASTA/FASTQ stats, QC, translation, ORF, "
                    "motif search).")
    parser.add_argument("--version", action="version",
                        version="%(prog)s " + __version__)
    sub = parser.add_subparsers(dest="command", metavar="<command>")

    p = sub.add_parser("stats", parents=[common],
                       help="length/GC/N50 summary statistics")
    p.add_argument("file")
    p.add_argument("--per-seq", action="store_true",
                   help="also print per-sequence statistics")
    p.add_argument("--json", action="store_true", help="JSON output")
    p.set_defaults(func=cmd_stats)

    p = sub.add_parser("validate", parents=[common],
                       help="validate FASTA/FASTQ structure")
    p.add_argument("file")
    p.add_argument("--json", action="store_true", help="JSON output")
    p.set_defaults(func=cmd_validate)

    p = sub.add_parser("gc", parents=[common],
                       help="per-sequence GC content and GC skew")
    p.add_argument("file")
    p.add_argument("--window", type=int, default=None,
                   help="sliding-window size in bases")
    p.add_argument("--step", type=int, default=None,
                   help="sliding-window step (default: window size)")
    p.add_argument("--max-windows", type=int, default=20,
                   help="maximum windows printed per sequence")
    p.add_argument("--json", action="store_true", help="JSON output")
    p.set_defaults(func=cmd_gc)

    p = sub.add_parser("revcomp", parents=[common],
                       help="reverse complement sequences")
    p.add_argument("file")
    p.add_argument("--rna", action="store_true",
                   help="treat sequences as RNA (complement A<->U)")
    p.set_defaults(func=cmd_revcomp)

    p = sub.add_parser("transcribe", parents=[common],
                       help="transcribe DNA to RNA")
    p.add_argument("file")
    p.set_defaults(func=cmd_transcribe)

    p = sub.add_parser("translate", parents=[common],
                       help="translate DNA/RNA to protein")
    p.add_argument("file")
    p.add_argument("--table", type=int, default=1, choices=sorted(GENETIC_CODES),
                   help="NCBI genetic code table (default: 1 Standard)")
    p.add_argument("--frames", type=int, default=1, choices=[1, 3, 6],
                   help="translate 1, 3 or 6 reading frames (default: 1)")
    p.set_defaults(func=cmd_translate)

    p = sub.add_parser("orf", parents=[common],
                       help="find open reading frames (both strands)")
    p.add_argument("file")
    p.add_argument("--min-len", type=int, default=30,
                   help="minimum ORF length in bp (default: 30)")
    p.add_argument("--table", type=int, default=1, choices=sorted(GENETIC_CODES),
                   help="NCBI genetic code table (default: 1 Standard)")
    p.add_argument("--forward-only", action="store_true",
                   help="scan the forward strand only")
    p.add_argument("--json", action="store_true", help="JSON output")
    p.set_defaults(func=cmd_orf)

    p = sub.add_parser("find", parents=[common],
                       help="find an IUPAC motif in sequences")
    p.add_argument("file")
    p.add_argument("motif", help="IUPAC DNA motif, e.g. GASTY or TATAAA")
    p.add_argument("--both-strands", action="store_true",
                   help="also match the reverse complement")
    p.add_argument("--json", action="store_true", help="JSON output")
    p.set_defaults(func=cmd_find)

    p = sub.add_parser("enzymes", parents=[common],
                       help="restriction enzyme sites (or --list)")
    p.add_argument("file", nargs="?", default=None)
    p.add_argument("--list", action="store_true",
                   help="list built-in enzymes and exit")
    p.add_argument("--enzyme", default=None,
                   help="comma-separated enzyme names (default: all)")
    p.add_argument("--json", action="store_true", help="JSON output")
    p.set_defaults(func=cmd_enzymes)

    p = sub.add_parser("tm", parents=[common],
                       help="primer melting temperature")
    p.add_argument("primers", nargs="+", help="one or more primer sequences")
    p.add_argument("--na", type=float, default=0.05,
                   help="sodium concentration in mol/L (default: 0.05)")
    p.add_argument("--json", action="store_true", help="JSON output")
    p.set_defaults(func=cmd_tm)

    p = sub.add_parser("fqstat", parents=[common],
                       help="FASTQ quality statistics")
    p.add_argument("file")
    p.add_argument("--max-positions", type=int, default=50,
                   help="maximum per-base positions printed")
    p.add_argument("--json", action="store_true", help="JSON output")
    p.set_defaults(func=cmd_fqstat)

    p = sub.add_parser("fqtrim", parents=[common],
                       help="trim and filter FASTQ reads")
    p.add_argument("file")
    p.add_argument("--quality", type=int, default=20,
                   help="minimum mean quality per window (default: 20)")
    p.add_argument("--window", type=int, default=4,
                   help="sliding window size (default: 4)")
    p.add_argument("--min-len", type=int, default=0,
                   help="discard reads shorter than this")
    p.add_argument("--max-len", type=int, default=None,
                   help="discard reads longer than this")
    p.add_argument("--min-mean-q", type=float, default=0.0,
                   help="discard reads with mean quality below this")
    p.add_argument("--max-n", type=float, default=1.0,
                   help="maximum fraction of N bases (default: 1.0)")
    p.set_defaults(func=cmd_fqtrim)

    p = sub.add_parser("convert", parents=[common],
                       help="convert between FASTA/FASTQ/CSV/TSV/NDJSON")
    p.add_argument("file")
    p.add_argument("--to", required=True,
                   choices=["fasta", "fastq", "csv", "tsv", "ndjson"],
                   help="target format")
    p.set_defaults(func=cmd_convert)

    p = sub.add_parser("sample", parents=[common],
                       help="randomly sample records (reproducible)")
    p.add_argument("file")
    p.add_argument("-n", "--number", type=int, required=True,
                   help="number of records to sample")
    p.add_argument("--seed", type=int, default=42,
                   help="random seed (default: 42)")
    p.set_defaults(func=cmd_sample)

    p = sub.add_parser("dedup", parents=[common],
                       help="remove duplicate sequences (case-insensitive)")
    p.add_argument("file")
    p.set_defaults(func=cmd_dedup)

    p = sub.add_parser("info", parents=[common],
                       help="detect format and preview records")
    p.add_argument("file")
    p.set_defaults(func=cmd_info)

    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if not getattr(args, "func", None):
        parser.print_help()
        return 1
    try:
        return args.func(args)
    except SeqForgeError as exc:
        sys.stderr.write("error: %s\n" % exc)
        return 1
    except ValueError as exc:
        sys.stderr.write("error: %s\n" % exc)
        return 1
    except IOError as exc:
        sys.stderr.write("error: %s\n" % exc)
        return 1


if __name__ == "__main__":
    sys.exit(main())
