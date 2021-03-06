% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/motif_all.R
\name{motif_all}
\alias{motif_all}
\title{Find ALL motifs in a list of sequences using motif-all algorithm}
\usage{
motif_all(seqs, min.seqs, genes = NULL, ncores = 1, verbose = FALSE,
  center = "STY")
}
\arguments{
\item{seqs}{a character vector of sequences}

\item{min.seqs}{minimum number of motif}

\item{genes}{the gene name from where a sequences is discovered, a character vector
has the same length as \code{seqs}}

\item{ncores}{the number of cores to be used, passed to \code{mclapply}.}

\item{verbose}{logical, whether print detailed information}

\item{center}{the amino acid centered at the sequences, sequences with other center AA would be 
removed from the list. To disable this function, set \code{center = NULL}.}
}
\value{
a list consists of:
  $mw - motif wise count, a names integer vector. The names are the sequence motifs and
  the integers indicate the frequency of each motifs in the input sequences.
  $gw - gene wise count, a list of two elements: 1) unique gene and 2) motif genes, i.e.
  genes include a specific type of motif
}
\description{
Find ALL motifs in a list of sequences using motif-all algorithm
}
\references{
He, Zengyou, Can Yang, Guangyu Guo, Ning Li, and Weichuan Yu. 2011.
  "Motif-All: Discovering All Phosphorylation Motifs." BMC Bioinformatics
  12 Suppl 1 (February):S22.
}
\author{
Chen Meng
}
