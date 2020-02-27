#!/usr/bin/env python3

import json
import sys
import argparse
from collections import namedtuple, defaultdict

from gffpal.gff import GFF
from gffpal.attributes import GFF3Attributes


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


def all_children_failed(parent):
    for child in parent.traverse_children():
        if child == parent:
            continue
        # As soon as we hit one that didn't fail, do early exit.
        elif ((child.attributes is not None) and
                ("is_unreliable" not in child.attributes.custom)):
            return False

    return True


def deal_with_kids(children, is_unreliable, type_, length, aligned):
    for child in children:

        if child.attributes is None:
            child.attributes = GFF3Attributes()

        if is_unreliable:
            child.attributes.custom["is_unreliable"] = "true"

        if child.type != type_:
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

    return length, aligned


def find_coverages(aligned, length):
    return {
        k: (v / length)
        for k, v
        in aligned.items()
    }


def main():
    args = cli(sys.argv[0], sys.argv[1:])
    gff = GFF.parse(args.infile)

    for mrna in gff.select_type(args.group_level):

        if mrna.attributes is None:
            mrna.attributes = GFF3Attributes()

        failed_antifam = "antifam_match" in mrna.attributes.custom
        if failed_antifam:
            mrna.attributes.custom["is_unreliable"] = "true"

        length, aligned = deal_with_kids(
            gff.traverse_children([mrna]),
            failed_antifam,
            args.type,
            0,
            defaultdict(int)
        )

        coverages = find_coverages(aligned, length)
        supported = [k for k, v in coverages.items() if v > args.min_cov]

        lineage = list(gff.traverse_children([mrna], sort=True))

        if ((len(supported) == 0) or
                (all(s in args.exclude for s in supported))):
            sup = True
            for f in lineage:
                if f.attributes is None:
                    f.attributes = GFF3Attributes.empty()

                f.attributes.custom["is_unreliable"] = "true"
        else:
            sup = False

        if args.stats:
            line = coverages
            line["id"] = mrna.attributes.id
            line["length"] = length
            line["filtered"] = sup
            print(json.dumps(line), file=args.stats)

    seen = set()
    for k in gff:
        if (len(k.parents) == 0) and all_children_failed(k):
            if k.attributes is None:
                k.attributes = GFF3Attributes.empty()
            k.attributes.custom["is_unreliable"] = "true"

        if k in seen:
            continue
        else:
            print(k, file=args.outfile)
            seen.add(k)

    return


if __name__ == "__main__":
    main()
