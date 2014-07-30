\name{gseseq}
\alias{gseseq}

\title{
    Gene Set Enrichment for Sequence Data
}
\description{

    This function implements a naive gene set enrichment method for
    sequence data. It is a toy used during a Bioc2014 Developer Day
    exercise.

}
\usage{
gseseq(genes, isDE, weights, sets)
}
\arguments{

  \item{genes}{\code{character(n)} vector of gene symbols.}

  \item{isDE}{\code{logical(n)} vector indicating statistical
      significance of gene differential expression.}

  \item{weights}{\code{numeric(n)} vector of non-negative values (e.g.,
      read depth) used to weigth contribution of each gene to the set(s)
      it belongs to.}

  \item{sets}{\code{list()} of gene sets, with each element a
      \code{character()} vector gene symbols in the set.}

}

\value{

    A \code{data.frame} listing, for each gene set, the (weighted)
    probability of over- or under-representation of differentially
    expressed genes in the set, relative to the universe of genes.

}

\author{
    Martin Morgan \url{mtmorgan@fhcrc.org}
}

\seealso{
    The Bioconductor goseq package
}

%% \examples{
%% data(gseseq)     # genes, isDE< weights, sets
%% gse <- gseseq(genes, isDE, weights, sets)
%% }

\keyword{manip}
