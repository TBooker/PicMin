% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PicMin_functions.R
\name{GenerateNullData}
\alias{GenerateNullData}
\title{Generate data under the null hypothesis}
\usage{
GenerateNullData(adaptation_screen, a, b, n, genes)
}
\arguments{
\item{adaptation_screen}{The threshold used to determine adaptation}

\item{a}{the 'a' parameter of a beta distribution of p-values for the false null}

\item{b}{the 'b' parameter of a beta distribution of p-values for the false null}

\item{n}{the number of species in the test}

\item{genes}{the number of genes in the genome use to calculate empirical p-values}
}
\description{
Generate data under the null hypothesis
}
