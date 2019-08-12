#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/


def helpMessage() {
    log.info"""
    # panann

    A pipeline to predict genes in a multiple fungal genomes.

    ## Usage

    ```bash


    ## Exit codes

    - 0: All ok.
    - 1: Incomplete parameter inputs.


    ## Output
    aligned_reads - Fastq reads aligned to each genome as bams.
    """.stripIndent()
}

if (params.help) {
    helpMessage()
    exit 0
}

params.genome = false
params.annotations = false
params.add_utr = false
params.nocrf = false
