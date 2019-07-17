#!/usr/bin/env python3

import sys
import argparse
from typing import NamedTuple
from typing import Optional
from typing import TypeVar
from typing import ClassVar
from typing import List
from typing import Dict

from enum import Enum

FILE_FORMATS = ["gene", "match"]
SOURCES = ["M", "E", "P", "RM", "W", "XNT", "C", "D", "T", "R", "PB"]

GFF_TYPE_MAP = {
    "transcription_start_site": "transcription_start_site",
    "TSS": "transcription_start_site",
    "SO:0000315": "transcription_start_site",
    "transcription_end_site": "transcription_end_site",
    "SO:0000616": "transcription_end_site",
    "exon": "exon",
    "SO:0000147": "exon",
    "interior_coding_exon": "exon",
    "SO:0000004": "exon",
    "coding_exon": "exon",
    "SO:0000195": "exon",
    "five_prime_coding_exon": "exon",
    "SO:0000200": "exon",
    "three_prime_coding_exon": "exon",
    "SO:0000202": "exon",
    "interior_exon": "exon",
    "SO:0000201": "exon",
    "UTR": "UTR",
    "SO:0000203": "UTR",
    "noncoding_exon": "UTR",
    "SO:0000198": "UTR",
    "three_prime_UTR": "three_prime_UTR",
    "SO:0000205": "three_prime_UTR",
    "five_prime_UTR": "five_prime_UTR",
    "SO:0000204": "five_prime_UTR",
    "three_prime_noncoding_exon": "three_prime_UTR",
    "SO:0000444": "three_prime_UTR",
    "five_prime_noncoding_exon": "five_prime_UTR",
    "SO:0000445": "five_prime_UTR",
    "intron": "intron",
    "SO:0000188": "intron",
    "five_prime_intron": "intron",
    "SO:0000190": "intron",
    "interior_intron": "intron",
    "SO:0000191": "intron",
    "three_prime_intron": "intron",
    "SO:0000192": "intron",
    "UTR_intron": "intron",
    "SO:0000446": "intron",
    "SO:0000447": "intron",
    "five_prime_UTR_intron": "intron",
    "SO:0000448": "intron",
    "three_prime_UTR_intron": "intron",
    "CDS": "CDS",
    "coding_sequence": "CDS",
    "SO:0000316": "CDS",
    "CDS_fragment": "CDSpart",
    "SO:0001384": "CDSpart",
    "CDS_supported_by_peptide_spectrum_match": "CDSpart",
    "SO:0002071": "CDSpart",
    "start_codon": "start_codon",
    "SO:0000318": "start_codon",
    "non_canonical_start_codon": "start_codon",
    "SO:0000680": "start_codon",
    "stop_codon": "stop_codon",
    "SO:0000319": "stop_codon",
    "five_prime_cis_splice_site": "five_prime_cis_splice_site",
    "SO:0000163": "five_prime_cis_splice_site",
    "donor splice site": "five_prime_cis_splice_site",
    "five prime splice site": "five_prime_cis_splice_site",
    "splice donor site": "five_prime_cis_splice_site",
    "canonical_five_prime_splice_site": "five_prime_cis_splice_site",
    "SO:0000677": "five_prime_cis_splice_site",
    "non_canonical_five_prime_splice_site": "five_prime_cis_splice_site",
    "SO:0000679": "five_prime_cis_splice_site",
    "three_prime_cis_splice_site": "three_prime_cis_splice_site",
    "SO:0000164": "three_prime_cis_splice_site",
    "acceptor splice site": "three_prime_cis_splice_site",
    "splice acceptor site": "three_prime_cis_splice_site",
    "three prime splice site": "three_prime_cis_splice_site",
    "canonical_three_prime_splice_site": "three_prime_cis_splice_site",
    "SO:0000676": "three_prime_cis_splice_site",
    "non_canonical_three_prime_splice_site": "three_prime_cis_splice_site",
    "SO:0000678": "three_prime_cis_splice_site",
    "match": "CDS",
    "SO:0000343": "CDS",
    "protein_match": "CDS",
    "SO:0000349": "CDS",
    "nucleotide_to_protein_match": "CDS",
    "translated_nucleotide_match": "CDS",
    "SO:0000181": "CDS",
    "nucleotide_match": "exon",
    "SO:0000347": "exon",
    "cDNA_match": "exon",
    "SO:0000689": "exon",
}

HINT_TYPE = [
    "start",
    "stop",
    "tss",
    "tts",
    "ass",
    "dss",
    "exonpart",
    "exon",
    "intronpart",
    "intron",
    "CDSpart",
    "CDS",
    "UTRpart",
    "UTR",
    "irpart",
    "nonexonpart",
    "genicpart",
]


StrandT = TypeVar("StrandT", bound=Strand)

class Strand(Enum):
    PLUS = 0
    MINUS = 1
    UNSTRANDED = 2
    UNKNOWN = 3

    into_str_map: ClassVar[List[str]] = ["+", "-", ".", "?"]
    from_str_map: ClassVar[Dict[str, StrandT]] = {
        "+": cls.PLUS,
        "-": cls.MINUS,
        ".": cls.UNSTRANDED,
        "?": cls.UNKNOWN,
    }

    def __str__(self):
        return into_str_map[self.value]

    @classmethod
    def from_str(cls: StrandT, string: str) -> StrandT:
        try:
            return from_str_map[string]
        except KeyError:
            valid = list(from_str_map.keys())
            raise KeyError(f"Invalid option. Must be one of {valid}")


class Phase(Enum):
    FIRST = 0
    SECOND = 1
    THIRD = 2
    NOT_CDS = 3

    into_str_map: ClassVar[List[str]] = ["0", "1", "2", "3"]
    from_str_map: ClassVar[Dict[str, StrandT]] = {
        "0": cls.FIRST,
        "1": cls.SECOND,
        "2": cls.THIRD,
        "0": cls.NOT_CDS,
    }

    def __str__(self):
        return into_str_map[self.value]

    @classmethod
    def from_str(cls: StrandT, string: str) -> StrandT:
        try:
            return from_str_map[string]
        except KeyError:
            valid = list(from_str_map.keys())
            raise KeyError(f"Invalid option. Must be one of {valid}")

class GFFRecord(NamedTuple):

    seqid: str
    source: str
    type: str
    start: int
    end: int
    score: Optional[float] = None
    strand: Strand = Strand.UNSTRANDED
    phase: Phase = Phase.NOT_CDS



def cli(prog, args):
    parser = argparse.ArgumentParser(
        prog=prog,
        description="""
        """
    )

    parser.add_argument(
        "infile",
        default=sys.stdin,
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
        "-f", "--format",
        default="gene",
        choices=FILE_FORMATS,
        help="The format of the input gff.",
    )

    parser.add_argument(
        "-s", "--source",
        default="M",
        help=f"The type of hint to create. Usually one of {SOURCES}.",
    )

    parser.add_argument(
        "-c", "--cds",
        default="CDSpart",
        choices=HINT_TYPE,
        help="The type to map CDS features to."
    )

    parser.add_argument(
        "-i", "--intron",
        default="intron",
        choices=HINT_TYPE,
        help="The type to map intron features to."
    )

    parser.add_argument(
        "-e", "--exon",
        default="exonpart",
        choices=HINT_TYPE,
        help="The type to map exon features to."
    )

    parser.add_argument(
        "-5", "--utr5",
        default="UTRpart",
        choices=HINT_TYPE,
        help="The type to map five_prime_UTR features to."
    )

    parser.add_argument(
        "-3", "--utr3",
        default="UTRpart",
        choices=HINT_TYPE,
        help="The type to map three_prime_UTR features to."
    )

    parser.add_argument(
        "-u", "--utr",
        default="UTRpart",
        choices=HINT_TYPE,
        help="The type to map UTR features to."
    )

    return parser.parse_args(args)


def main():
    args = cli(prog=sys.argv[0], args=sys.argv[1:])

    gff_to_hints = {
        "CDS": args.cds,
        "CDSpart": "CDSpart",
        "intron": args.intron,
        "intronpart": "intronpart",
        "exon": args.exon,
        "exonpart": "exonpart",
        "five_prime_UTR": args.utr5,
        "three_prime_UTR": args.utr3,
        "UTR": args.utr,
        "UTRpart": "UTRpart",
        "transcription_start_site": "tss",
        "transcription_end_site": "tts",
        "start_codon": "start",
        "stop_codon": "stop",
        "five_prime_cis_splice_site": "dss",
        "three_prime_cis_splice_site": "ass",
        "irpart": "irpart",
        "nonexonpart": "nonexonpart",
        "genicpart": "genicpart",
    }

    return

if __name__ == "__main__":
    main()
