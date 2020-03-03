#!/usr/bin/env python3

import json
import sys
import argparse
from copy import deepcopy
from collections import namedtuple, defaultdict

from intervaltree import Interval, IntervalTree

from typing import List, Set
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
    type_: str,
    length: int,
    aligned: DefaultDict[str, int],
) -> Tuple[int, DefaultDict[str, int]]:
    for child in children:

        if child.attributes is None:
            child.attributes = GFF3Attributes()

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


def write_gff(
    gff: GFF,
    outfile: TextIO,
) -> None:
    print("#gff-version 3", file=outfile)
    for feature in gff.traverse_children(sort=True, breadth=True):
        print(feature, file=outfile)
    return


def split_gffs(gff: GFF, group_level: str) -> Tuple[GFF, GFF]:
    keep: Set[GFF3Record] = set()
    drop: Set[GFF3Record] = set()

    for mrna in gff.select_type(group_level):
        if mrna.attributes is None:
            should_keep = True
        elif "should_exclude" in mrna.attributes.custom:
            should_keep = False
        else:
            should_keep = True

        if should_keep:
            keep.update(gff.traverse_children([mrna]))
            keep.update(gff.traverse_parents([mrna]))
        else:
            drop.update(gff.traverse_children([mrna]))
            drop.update(gff.traverse_parents([mrna]))

    kept = prune_gff(keep)
    dropped = prune_gff(drop)
    return kept, dropped


def prune_gff(records: Set[GFF3Record]) -> GFF:

    new_records: Dict[GFF3Record, GFF3Record] = dict()

    # Create a mapping from old to new objects
    # to preserve hashing/lookup capability
    for record in records:
        new_record = deepcopy(record)
        new_record.children = []
        new_record.parents = []

        new_records[record] = new_record

    for record in records:
        new_record = new_records[record]

        for parent in record.parents:
            # Don't add parents that shouldn't be in this set
            if parent not in records:
                continue

            new_parent = new_records[parent]
            new_record.add_parent(new_parent)

        for child in record.children:
            # Don't add children that shouldn't be in this set
            if child not in records:
                continue

            new_child = new_records[child]
            new_record.add_child(new_child)

        # Update the record parent IDS to reflect the new split set.
        if new_record.attributes is not None:
            new_record.attributes.parent = []
            for parent in new_record.parents:
                # This should always be true, as the ID is necessary
                assert parent.attributes is not None
                assert parent.attributes.id is not None
                new_record.attributes.parent.append(parent.attributes.id)
        else:
            # This necessarily should be true, since attributes define parent
            # child relationships.
            assert len(new_record.children) == 0
            assert len(new_record.parents) == 0

    return GFF(new_records.values())


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
            mrna.attributes.custom["is_unreliable"] = "true"
            if not is_novel:
                mrna.attributes.custom["should_exclude"] = "true"

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

    kept, dropped = split_gffs(gff, args.group_level)

    write_gff(kept, args.outfile)

    if args.filtered is not None:
        write_gff(dropped, args.filtered)
    return


if __name__ == "__main__":
    main()
