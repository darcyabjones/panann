#!/usr/bin/env Rscript

VERSION="0.0.1"


suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("ape"))
suppressPackageStartupMessages(library("phangorn"))
suppressPackageStartupMessages(library("magrittr"))

option_list <- list(
    make_option(
        c("-i", "--infile"),
        type="character",
        action="store",
        help="The distance matrix to use (required)."
    ),
    make_option(
        c("-o", "--outfile"),
        type="character",
        action="store",
        help="The newick file to create."
    ),
    make_option(
        c("--nj"),
        type="logical",
        action="store_true",
        default=FALSE,
        help="Use neighbor joining instead of UPGMA."
    ),
    make_option(
        "--version",
        type="logical",
        action="store_true",
        default=FALSE,
        help="Print version and exit.",
    )
)

parser <- OptionParser(
    usage = "%prog --infile my.tsv --outfile my.nwk",
    option_list = option_list
)

args <- parse_args(parser)

log_stderr <- function(...) {
  cat(sprintf(...), sep='', file=stderr())
}

quit_with_err <- function(...) {
  log_stderr(...)
  quit(save = "no", status = 1, runLast = FALSE)
}

validate_file <- function(path) {
  if (is.null(path)) {
    quit_with_err("Please provide required file")
  }
}

stripper <- function(str) {
  str <- gsub("^X", "", str, perl = TRUE)
  str
}

main <- function(args) {
  if (args$version) {
    cat(VERSION, file=stdout())
    quit(save = "no", status = 0, runLast = FALSE)
  }

  validate_file(args$infile)
  validate_file(args$outfile)

  table <- read.table(
    args$infile,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  row.names(table) <- table[, "sample"]
  table = subset(table, select = -sample)
  mat <- as.dist(as.matrix(table), upper = FALSE, diag = FALSE)

  if (args$nj) {
    tree <- NJ(mat)
  } else {
    tree <- upgma(mat)
  }

  write.tree(tree, args$outfile)
}

main(args)
