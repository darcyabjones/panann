#!/usr/bin/env python3

import sys
import argparse

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.CheckSum import seguid


class IdConverter(object):

    """
    Generates an id based on an integer.
    Essentially it's a small hashing algorithm that allows us to
    give unique ids that aren't too long.
    """

    def __init__(
            self,
            state: int = 0,
            prefix: str = "",
            length: int = 4,
            alphabet: str = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ',
            ) -> None:
        """
        """

        self.state = state
        self.prefix = prefix
        self.length = length
        self.alphabet = alphabet

        from baseconv import BaseConverter
        self.converter = BaseConverter(alphabet)
        return

    def encode(self, number: int) -> str:
        """
        """
        template = "{pre}{{p:{first}>{length}}}".format(
            pre=self.prefix,
            first=self.converter.digits[0],
            length=self.length,
            )
        return template.format(p=self.converter.encode(number))

    def decode(self, pattern: str) -> int:
        """
        """
        pattern = pattern[len(self.prefix):]
        return int(self.converter.decode(pattern))

    def __next__(self):
        string = self.encode(self.state)
        self.state += 1
        return string

    def __iter__(self):
        return self


def cli(prog, args):
    parser = argparse.ArgumentParser(
        prog=prog,
        description="""
        """
    )

    parser.add_argument(
        "-i", "--infiles",
        required=True,
        nargs="+",
        type=argparse.FileType('r'),
        help="Input fasta files. Use '-' for stdin.",
    )

    parser.add_argument(
        "-o", "--outfile",
        type=argparse.FileType('w'),
        default=sys.stdout,
        help="Output fasta file path. Default stdout.",
    )

    parser.add_argument(
        "-m", "--map",
        type=argparse.FileType('w'),
        default=None,
        help="Write mapping from original file and id to new id.",
    )

    parser.add_argument(
        "-p", "--prefix",
        default="PA_",
        help="Prefix to add to new seqids.",
    )

    return parser.parse_args(args)


def main():
    args = cli(prog=sys.argv[0], args=sys.argv[1:])

    id_counter = 1
    id_conv = IdConverter(prefix=args.prefix, length=6)
    seen = dict()

    for infile in args.infiles:
        for record in SeqIO.parse(infile, format="fasta"):
            seq = record.seq.rstrip("*").upper()
            chk = seguid(str(seq))
            if chk in seen:
                new_id = seen[chk]
            else:
                new_id = id_conv.encode(id_counter)
                id_counter += 1
                new_record = SeqRecord(
                    id=new_id,
                    name=new_id,
                    description=new_id,
                    seq=seq
                )
                SeqIO.write(new_record, args.outfile, format="fasta")

            if args.map is not None:
                print(
                    "{}\t{}\t{}".format(infile.name, record.id, new_id),
                    file=args.map
                )
    return


if __name__ == "__main__":
    main()
