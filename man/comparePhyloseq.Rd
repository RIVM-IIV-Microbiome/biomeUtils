% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/comparePhyloseq.R
\name{comparePhyloseq}
\alias{comparePhyloseq}
\title{Compare Phyloseq Objects}
\usage{
comparePhyloseq(x)
}
\arguments{
\item{x}{A named list of \code{phyloseq} objects to compare}
}
\value{
A tibble
}
\description{
Compare Phyloseq Objects
}
\details{
Provided a named list of different \code{phyloseq} objects
         `comparePhyloseq` returns a data.frame with basic values in
         each of the \code{phyloseq} objects.
}
\examples{
library(biomeUtils)
data("FuentesIliGutData")
ps1 <- subset_samples(FuentesIliGutData, ILI == "C")
ps1 <- prune_taxa(taxa_sums(ps1) > 0, ps1)

ps2 <- subset_samples(FuentesIliGutData, ILI == "L1")
ps2 <- prune_taxa(taxa_sums(ps2) > 0, ps2)

ps3 <- subset_samples(FuentesIliGutData, ILI == "L2")
ps3 <- prune_taxa(taxa_sums(ps3) > 0, ps3)

ps.list <- c("C" = ps1, "L1" = ps2, "L2" = ps3)

comparePhyloseq(ps.list)
}
\references{
Shetty SA (2021). Utilities for microbiome analytics.
\url{https://github.com/RIVM-IIV-Microbiome/biomeUtils}
}
\author{
Sudarshan A. Shetty
}
