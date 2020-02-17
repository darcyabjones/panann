#!/usr/bin/env nextflow

include symmetric_difference from './utils'

// Handy little utility that can replace some if statements in closures.
is_null = { f -> (f == null || f == '') }


// Valid values for the analysis column in the config table.
valid_table_analyses = [
    "genome",
    "known",
    "fastq_forward",
    "fastq_reverse",
    "gmap",
    "spaln_transcripts",
    "spaln_proteins",
    "exonerate",
    "cram",
    "stringtie",
    "augustus_intron_hints",
    "genemark",
    "pasa",
    "codingquarry",
    "codingquarrypm",
    "gemoma_intron_hints",
    "gemoma_forward_coverage",
    "gemoma_reverse_coverage",
    "gemoma",
    "augustus",
    "gemoma_comparative",
    "evm",
    "augustus_gapfiller",
    "augustus_hints",
    "evm_transcript_hints",
    "evm_protein_hints",
    "final"
]


/**
 * Raises a generic error if we hit a branch that should be controlled for.
 * This means we fail early rather than some time downstream.
 */
def param_unexpected_error() {
    log.error "We encountered an error while validating input arguments that " +
        "should be possible. Please raise an issue on github or contact the " +
        "authors."
    exit 1
}


// This probably needs to be moved to each pipeline script.
def validate_params(params) {
    def should_run = [:]

    def error = false

    if ( false ) {
        log.error "message"
        error = true
    }

    if (error) {
        exit 1
    }

    return should_run
}


/**
 * Convenience function to convert a filepath to a value channel.
 * Expects the filepath to not be empty.
 *
 * @param filepath The path to the file to get a channel for.
 * @return Value channel containing a single file object.
 */
def get_file(filepath) {
    if ( filepath ) {
        handle = Channel.value( file(filepath, checkIfExists: true) )
    } else {
        // The expectation is that you would check that filepath is not false
        // before you use this.
        param_unexpected_error()
    }

    return handle
}


/**
 * Convenience function to raise a user error if a required table value was
 * empty.
 *
 * @param inval The value in the table to test if empty.
 * @param field The column name in the table that shouldn't be empty.
 * @param line The full line as a string or hashmap to print so that people
 *             can find the problematic line.
 * @return The input value unaltered. This makes it convenient to use at
 *         the same time as assignment or while passing to other functions.
 */
def table_not_null(inval, field, line) {
    if (is_null(inval)) {
        log.error "A line supplied to `--table` had no ${field} field specified."
        log.error "This value is required."
        log.error "The line was: ${line}"

        exit 1
    }

    return inval
}


/**
 * Convert a cell value filename to a file object, raising a user error if
 * the cell is empty.
 *
 * @param filepath The actual filepath/cell value to load.
 * @param line The whole line so that people can find the issue.
 * @param File The filepath as a file object.
 */
def table_get_file(filepath, line) {
    if ( filepath ) {
        handle = file(filepath, checkIfExists: true)
    } else {
        log.error 'A line supplied to `--table` had no file field.'
        log.error 'The line was: ${it}'

        exit 1
    }

    return handle
}


/**
 * Raise a user error if the value isn't one of a list of options.
 * Used to control user input.
 *
 * @param inval The value that should be one of several options.
 * @param options A list of options that inval can be one of.
 * @param field The column name in the table that shouldn't be empty.
 * @param line The full line as a string or hashmap to print so that people
 *             can find the problematic line.
 * @return The input value unaltered. This makes it convenient to use at
 *         the same time as assignment or while passing to other functions.
 */
def table_one_of(inval, options, field, line) {
    if ( is_null ) {
        return inval
    }

    if ( inval in options ) {
        return inval
    } else {
        log.error "A line supplied to `--table` used an invalid choice for the ${field} field."
        log.error "The line was: ${line}"

        exit 1
    }
}


/**
 * Parse the config table into separate channels.
 *
 * @param genomes A channel containing the input genome fastas.
 *                Should have the structure: tuple(val(name), path(fasta))
 * @param table A channel containing the table tsv file.
 * @return HashMap of channels keyed by analysis type.
 */
def handle_table(genomes, table) {

    // First we validate that the fields exist and fork
    // the channel into many based on analysis type.
    branched = table
        .splitCsv(header: true, sep: '\t', strip: true)
        .map { it ->
            def line = [:]

            // Analysis column always is reqiured and has a restricted set
            // of options.
            line.analysis = table_one_of(
                table_not_null(it.analysis, "analysis", it),
                valid_table_analyses,
                "analysis",
                it
            )

            // The file is always reqiured.
            line.file = table_get_file(
                table_not_null(it.file, "file", it),
                it
            )

            // Name is required for all analyses except fastq files.
            if ( line.analysis in ["fastq_forward", "fastq_reverse"] ) {
                line.name = it.name
            } else {
                line.name = table_not_null(it.name, "name", it)
            }

            // The strand must be provided for fastq and cram input.
            if ( line.analysis in ["fastq_forward", "fastq_reverse", "cram"]) {
                line.strand = table_one_of(
                    table_not_null(it.strand, "strand", it),
                    ["fr", "rf"],
                    "strand",
                    it
                )
            }

            // The read_group must be provided for fastq and cram input.
            if ( line.analysis in ["fastq_forward", "fastq_reverse", "cram"]) {
                line.read_group = table_not_null(it.read_group, "read_group", it)
            }

            return line
        }
        .branch {
            genome: it.analysis == "genome"
            known: it.analysis == "known"
            fastq_forward: it.analysis == "fastq_forward"
            fastq_reverse: it.analysis == "fastq_reverse"
            gmap: it.analysis == "gmap"
            spaln_transcripts: it.analysis == "spaln_transcripts"
            spaln_proteins: it.analysis == "spaln_proteins"
            exonerate: it.analysis == "exonerate"
            cram: it.analysis == "cram"
            stringtie: it.analysis == "stringtie"
            augustus_intron_hints: it.analysis == "augustus_intron_hints"
            genemark: it.analysis == "genemark"
            pasa: it.analysis == "pasa"
            codingquarry: it.analysis == "codingquarry"
            codingquarrypm: it.analysis == "codingquarrypm"
            gemoma_intron_hints: it.analysis == "gemoma_intron_hints"
            gemoma_forward_coverage: it.analysis == "gemoma_forward_coverage"
            gemoma_reverse_coverage: it.analysis == "gemoma_reverse_coverage"
            gemoma: it.analysis == "gemoma"
            augustus: it.analysis == "augustus"
            gemoma_comparative: it.analysis == "gemoma_comparative"
            evm: it.analysis == "evm"
            augustus_gapfiller: it.analysis == "augustus_gapfiller"
            augustus_hints: it.analysis == "augustus_hints"
            evm_transcript_hints: it.analysis == "evm_transcript_hints"
            evm_protein_hints: it.analysis == "evm_protein_hints"
            final: it.analysis == "final"
        }


    // At this point we know that files exist and no fields that should be
    // filled are null.
    table_genomes = genomes.mix(branched.genome.map { [it.name, it.file] })

    // These are the GFFs
    known = branched.known.map { [it.name, it.file] }

    // Make a temporary channel to get strands for read_groups
    // and make sure that all strands for the read_group are the same.
    fastq_strands = branched.fastq_forward
        .mix(branched.fastq_reverse)
        .map { [it.read_group, it.strand] }
        .unique()
        .groupTuple(by: 0)
        .map { rg, strands ->
            if (strands.size() > 1) {
                log.error "The fastq read group ${rg} in the table had more " +
                          "than one strandedness specified."
                log.error "Please make sure they are all the same or split " +
                          "them into multiple read groups."
                exit 1
            }

            return [rg, strands[0]]
        }

    // Join the strands with the actual fastq channel.
    fastq = branched.fastq_forward
        .map { [it.read_group, it.file] }
        .groupTuple(by: 0)
        .join(
            branched.fastq_reverse
                .map { [it.read_group, it.file] }
                .groupTuple(by: 0),
            by: 0
        )
        .combine(fastq_strands, by: 0)

    cram = branched.cram.map { [it.name, it.read_group, it.file, it.strand] }

    // We'll return a hashmap of the channels to avoid crazy tuples.
    // It'll look like the output of .branch()
    out_channels = ["genome": table_genomes, "fastq": fastq, "cram": cram]

    // The other analysis types are all the same so we can just exclude
    // the columns: `analysis`, `read_group`, `strand`
    valid_table_analyses.each { analysis ->
        if ( !(analysis in ["genome", "fastq", "cram"]) ) {
            out_channels[analysis] = branched
                .getProperty(analysis)
                .map { [it.name, it.file] }
        }
    }

    return out_channels
}


/**
 * This checks if there are any values in either channel that are not present
 * in either channel.
 */
def assert_same_names(ch1, ch2, ch1_name, ch2_name, msg) {
    ch1.collect().toList().combine(ch2.collect().toList())
        .map { c1, c2 ->
            symdiff = symmetric_difference(c1, c2)
            if ( symdiff.size() != 0 ) {
                log.error "Some names for the ${ch1_name} and ${ch2_name} " +
                          "don't match up."
                log.error "They are: ${symdiff}."

                if (msg) {
                    log.error msg
                }

                exit 1
            }
        }
}


/**
 * Say you have a bunch of genomes, and you want to make sure you have
 * all of the alignments for the genomes, but you don't care if there are a
 * few extra alignment files that don't match.
 */
def assert_two_covers_one(ch1, ch2, ch1_name, ch2_name, msg) {

    ch1.collect().toList().combine(ch2.collect().toList())
        .map { c1, c2 ->
            diff = c1 - c2
            if ( diff.size() != 0 ) {
                log.error "Some names in the ${ch1_name} are not represented " +
                          "in ${ch2_name}."
                log.error "They are: ${diff}."

                if ( msg ) {
                    log.error msg
                }

                exit 1
            }
        }
}
