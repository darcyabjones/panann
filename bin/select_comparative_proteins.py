#!/usr/bin/env python3

import sys
import argparse

from Bio import SeqIO
import pandas as pd


def cli(prog, args):
    parser = argparse.ArgumentParser(
        prog=prog,
        description=""" """
    )

    parser.add_argument(
        "--clusters",
        required=True,
        type=argparse.FileType('r'),
        help="Input cluster tsv file.",
    )

    parser.add_argument(
        "--assignments",
        required=True,
        type=argparse.FileType('r'),
        help="Input assignments file.",
    )

    parser.add_argument(
        "--cdsparts",
        required=True,
        type=argparse.FileType('r'),
        help="Input cdsparts file.",
    )

    parser.add_argument(
        "--proteins",
        required=True,
        type=argparse.FileType('r'),
        help="Input proteins file.",
    )

    return parser.parse_args(args)


def main():
    args = cli(sys.argv[0], sys.argv[1:])

    clusters = pd.read_csv(
        args.clusters,
        sep="\t",
        names=["cluster", "member"]
    )["cluster"].unique()

    assignments = pd.read_csv(
        args.assignments,
        sep="\t"
    ).set_index("transcript", drop=False)

    filtered_assignments = assignments.loc[clusters, ]
    filtered_assignments.reset_index(drop=True, inplace=True)
    filtered_assignments.to_csv("assignment.tsv", index=False, sep="\t")

    proteins = SeqIO.parse(args.proteins, format="fasta")
    proteins_to_get = set(clusters)
    filtered_proteins = [p for p in proteins if p.id in proteins_to_get]
    SeqIO.write(filtered_proteins, "proteins.fasta", format="fasta")

    cdsparts = SeqIO.parse(args.cdsparts, format="fasta")

    cds_to_get = set()
    cds_parts_col = filtered_assignments["cds-parts"].str.split(", ")

    for gene, parts in zip(filtered_assignments["#geneID"], cds_parts_col):
        names = [f"{gene}_{p}" for p in parts]
        cds_to_get.update(names)

    filtered_cdss = [c for c in cdsparts if c.id in cds_to_get]
    SeqIO.write(filtered_cdss, "cdsparts.fasta", format="fasta")
    return


if __name__ == "__main__":
    main()
