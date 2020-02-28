#!/usr/bin/env python3

import json
import sys
import argparse
from collections import namedtuple, defaultdict

from intervaltree import Interval, IntervalTree

from typing import List
from typing import Optional
from typing import TextIO
from typing import Iterator
from typing import Mapping, DefaultDict, Dict
from typing import Tuple

from gffpal.gff import GFF, GFF3Record
from gffpal.attributes import GFF3Attributes


def cli(prog: str, args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog=prog,
        description=""" .
        """
    )

    parser.add_argument(
        "infile",
        type=argparse.FileType('r'),
        help="The input gff file. Use '-' for stdin.",
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
        "-f", "--filtered",
        default=None,
        type=argparse.FileType('w'),
        help="Write the removed genes results to this file.",
    )

    parser.add_argument(
        "-r", "--threshold",
        default=0.3,
        type=float,
        help="The threshold to use when considering if a locus if new.",
    )

    parser.add_argument(
        "-o", "--outfile",
        default=sys.stdout,
        type=argparse.FileType('w'),
        help="Output gff3 file path. Default stdout.",
    )

    return parser.parse_args(args)


Hint = namedtuple("Hint", ["source", "name", "start", "end", "score", "naln"])


def parse_hint(string: str) -> Hint:
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

    return Hint(source, name, start, end, score, naln)


def all_children_unreliable(parent: GFF3Record) -> bool:
    for child in parent.traverse_children():
        if child == parent:
            continue
        # As soon as we hit one that didn't fail, do early exit.
        elif ((child.attributes is not None) and
                ("is_unreliable" not in child.attributes.custom)):
            return False

    return True


def all_children_excluded(parent: GFF3Record) -> bool:
    for child in parent.traverse_children():
        if child == parent:
            continue
        # As soon as we hit one that didn't fail, do early exit.
        elif ((child.attributes is not None) and
                ("should_exclude" not in child.attributes.custom)):
            return False

    return True


def deal_with_kids(
    children: Iterator[GFF3Record],
    is_unreliable: bool,
    type_: str,
    length: int,
    aligned: DefaultDict[str, int],
) -> Tuple[int, DefaultDict[str, int]]:
    for child in children:

        if child.attributes is None:
            child.attributes = GFF3Attributes()

        if is_unreliable:
            child.attributes.custom["is_unreliable"] = "true"

        if child.type != type_:
            continue

        length += child.length()

        best_hints: Dict[str, Hint] = dict()
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


def find_coverages(
    aligned: Mapping[str, int],
    length: int
) -> Dict[str, float]:
    return {
        k: (v / length)
        for k, v
        in aligned.items()
    }


def gff_to_itree(records: Iterator[GFF3Record]) -> Dict[str, IntervalTree]:
    intervals: DefaultDict[str, List[Interval]] = defaultdict(list)

    for mrna in records:
        intervals[mrna.seqid].append(Interval(mrna.start, mrna.end, mrna))

    itree: Dict[str, IntervalTree] = {
        k: IntervalTree(v)
        for k, v
        in intervals.items()
    }
    return itree


def nintersection(lstart, lend, rstart, rend) -> int:
    assert lstart < lend
    assert rstart < rend
    return min([lend, rend]) - max([lstart, rstart])


def is_novel_locus(
    record: GFF3Record,
    itree: Dict[str, IntervalTree],
    threshold: float
) -> bool:
    overlaps = itree.overlaps(record.start, record.end)

    for overlap in overlaps:
        if overlap.data == record:
            continue

        noverlap = nintersection(
            record.start,
            record.end,
            overlap.begin,
            overlap.end
        )

        if (noverlap / len(overlap)) > threshold:
            return False

    return True


def write_results(
    gff: GFF,
    outfile: TextIO,
    filtered_file: Optional[TextIO]
) -> None:
    seen = set()
    for k in gff:
        if k in seen:
            continue

        if (len(k.parents) == 0) and all_children_unreliable(k):
            if k.attributes is None:
                k.attributes = GFF3Attributes.empty()
            k.attributes.custom["is_unreliable"] = "true"

        if (len(k.parents) == 0) and all_children_excluded(k):
            if k.attributes is None:
                k.attributes = GFF3Attributes.empty()
            k.attributes.custom["should_exclude"] = "true"

        if ((k.attributes is not None)
                and ("should_exclude" in k.attributes.custom)):
            if filtered_file is not None:
                print(k, file=filtered_file)
        else:
            print(k, file=outfile)

        seen.add(k)

    return


def main():
    args = cli(sys.argv[0], sys.argv[1:])
    gff = GFF.parse(args.infile)

    itree = gff_to_itree(gff.select_type(args.group_level))

    for mrna in gff.select_type(args.group_level):

        if mrna.attributes is None:
            mrna.attributes = GFF3Attributes()

        failed_antifam = "antifam_match" in mrna.attributes.custom
        if failed_antifam:
            mrna.attributes.custom["is_unreliable"] = "true"
            mrna.attributes.custom["should_exclude"] = "true"

        length, aligned = deal_with_kids(
            gff.traverse_children([mrna]),
            failed_antifam,
            args.type,
            0,
            defaultdict(int)
        )

        coverages = find_coverages(aligned, length)
        supported = [k for k, v in coverages.items() if v > args.min_cov]
        not_supported = ((len(supported) == 0) or
                         (all(s in args.exclude for s in supported)))

        is_novel = is_novel_locus(mrna, itree, args.threshold)

        if not_supported:
            lineage = list(gff.traverse_children([mrna], sort=True))
            for f in lineage:
                if f.attributes is None:
                    f.attributes = GFF3Attributes.empty()

                f.attributes.custom["is_unreliable"] = "true"
                if not is_novel:
                    f.attributes.custom["should_exclude"] = "true"

        if args.stats:
            line = coverages
            line["id"] = mrna.attributes.id
            line["length"] = length
            line["is_supported"] = not_supported
            line["is_novel_locus"] = is_novel
            line["antifam_match"] = failed_antifam
            line["excluded"] = (failed_antifam
                                or (not_supported and not is_novel))

            print(json.dumps(line), file=args.stats)

    write_results(gff, args.outfile, args.filtered_file)
    return


if __name__ == "__main__":
    main()
