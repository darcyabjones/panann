#!/usr/bin/env python3

import re
import sys

def main():
    if len(sys.argv) != 4:
        print("Usage: {} in.gff out.gtf both_utrs.txt".format(sys.argv[0]), file=sys.stderr)
        sys.exit(1)
    else:
        infile = sys.argv[1]
        outfile = sys.argv[2]
        both_file = sys.argv[3]

    id_regex = re.compile(r"ID=[^;]*(asmbl_\d+\.p\d+)")
    parent_regex = re.compile(r"Parent=[^;]*(asmbl_\d+\.p\d+)")
    split_regex = re.compile(r"\s+")

    with open(infile) as inhandle:

        complete = set()
        utr5 = set()
        utr3 = set()
        for line in inhandle:
            if line.startswith("#"):
                continue
            elif line.strip() == "":
                continue

            id_ = id_regex.search(line)
            if id_ is not None:
                id_ = id_.groups()[0]

            parent = parent_regex.search(line)
            if parent is not None:
                parent =parent.groups()[0]

            if "complete" in line:
                if id_ is not None:
                    complete.add(id_)
                if parent is not None:
                    complete.add(parent)

            row_type = split_regex.split(line)
            if row_type[2] == "five_prime_UTR":
                if parent is not None:
                    utr5.add(parent)
            elif row_type[2] == "three_prime_UTR":
                if parent is not None:
                    utr3.add(parent)

    both_utrs = set.intersection(complete, utr5, utr3)
    with open(both_file, "w") as outhandle:
        print("\n".join(both_utrs), file=outhandle)

    with open(infile) as inhandle, open(outfile, "w") as outhandle:

        for line in inhandle:
            if not (("\tCDS\t" in line)
                    or ("\tthree_prime_UTR\t" in line)
                    or ("\tfive_prime_UTR\t" in line)):
                continue

            # Because using utrs, we don't need exons.
            # or ("\texon\t" in line)
            
            id_ = id_regex.search(line).groups()[0]
            parent = parent_regex.search(line)
            if parent is not None:
                parent = parent.groups()[0]

            if (id_ not in complete) and (parent is not None) and (parent not in complete):
                continue

            sline = line.split("\t")
            sline[8] = 'transcript_id "{}";'.format(id_)
            print("\t".join(sline), file=outhandle)

    return

if __name__ == "__main__":
    main()
