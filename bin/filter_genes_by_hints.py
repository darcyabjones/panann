#!/usr/bin/env python3

import json
import sys
import argparse
from collections import namedtuple, defaultdict

from gffpal.gff import GFF


def cli(prog, args):
    parser = argparse.ArgumentParser(
        prog=prog,
        description=""" .
        """
    )

    parser.add_argument(
        "infile",
        type=argparse.FileType('r'),
        help="Input gff files. Use '-' for stdin.",
    )

    parser.add_argument(
        "-c", "--min-cov",
        default=0.5,
        type=float,
        help=("The minimum amount a hint must cover a "
              "gene to be considered as support."),
    )

    parser.add_argument(
        "-t", "--type",
        default="CDS",
        type=str,
        help="The type to filter based on coverage of.",
    )

    parser.add_argument(
        "-l", "--group-level",
        default="mRNA",
        type=str,
        help="The level to group genes by.",
    )

    parser.add_argument(
        "-e", "--exclude",
        default=[],
        nargs="+",
        type=str,
        help="Also exclude matches that are only supported by these sources.",
    )

    parser.add_argument(
        "-s", "--stats",
        default=None,
        type=argparse.FileType('w'),
        help="Write some stats about what was excluded to a file.",
    )

    parser.add_argument(
        "-o", "--outfile",
        default=sys.stdout,
        type=argparse.FileType('w'),
        help="Output gff3 file path. Default stdout.",
    )

    parser.add_argument(
        "-f", "--filtered",
        default=None,
        type=argparse.FileType('w'),
        help="Write filtered records to this file.",
    )

    return parser.parse_args(args)


hint = namedtuple("hint", ["source", "name", "start", "end", "score", "naln"])


def parse_hint(string):
    sstring = string.strip().split()
    assert len(sstring) == 6, string

    source = sstring[0]
    name = sstring[1]
    start = int(sstring[2])
    end = int(sstring[3])

    if sstring[4] == ".":
        score = None
    else:
        score = float(sstring[4])

    naln = int(sstring[5])

    return hint(source, name, start, end, score, naln)


def main():
    args = cli(sys.argv[0], sys.argv[1:])
    gff = GFF.parse(args.infile)

    kept = []
    not_kept = []

    for mrna in gff.select_type(args.group_level):
        length = 0
        aligned = defaultdict(int)

        for child in gff.traverse_children([mrna]):
            if child.type != args.type:
                continue

            length += child.length()

            best_hints = dict()
            hints = child.attributes.custom.get("hint", None)
            if hints is None:
                continue

            for h in hints.split(","):
                this_hint = parse_hint(h)

                if not ((this_hint.source in best_hints) and
                        (best_hints[this_hint.source].naln > this_hint.naln)):
                    best_hints[this_hint.source] = this_hint

            for s, h in best_hints.items():
                aligned[s] += h.naln

        supported = [
            k
            for k, v
            in aligned.items()
            if ((v / length) > args.min_cov)
        ]

        lineage = []
        lineage.extend(gff.traverse_parents([mrna], sort=True))
        # lineage.append(mrna)
        lineage.extend(gff.traverse_children([mrna], sort=True))

        if len(supported) == 0:
            sup = True
            not_kept.extend(lineage)
        elif all(s in args.exclude for s in supported):
            sup = True
            not_kept.extend(lineage)
        else:
            sup = False
            kept.extend(lineage)

        if args.stats:
            line = {k: (v / length)
                    for k, v in aligned.items()}
            line["id"] = mrna.attributes.id
            line["length"] = length
            line["filtered"] = sup
            print(json.dumps(line), file=args.stats)

    seen = set()
    for k in kept:
        if k in seen:
            continue
        else:
            print(k, file=args.outfile)
            seen.add(k)

    if args.filtered is not None:
        seen = set()
        for k in not_kept:
            if k in seen:
                continue
            else:
                print(k, file=args.filtered)
                seen.add(k)
    return


if __name__ == "__main__":
    main()
