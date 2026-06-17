\name{build_a_matrix}
\alias{build_a_matrix}
\title{Compute the additive relationship matrix (A) from a pedigree}
\usage{
build_a_matrix(ped_df)
}
\arguments{
\item{ped_df}{data.frame with columns id, sire, dam (as returned by build_ped)}
}
\value{
Named n-by-n matrix A
}
\description{
A[i,j] = 2 * kinship(i,j). Diagonal = 1 + F_i (inbreeding coefficient).
Founders are assumed non-inbred and unrelated to each other.
}
\details{
Uses pedigreemm::getA() when the package is installed; otherwise falls back to
the Henderson (1976) tabular method (O(n^2), suitable for n < ~5000).
}
\examples{
ped <- build_ped(n_gen = 1, seed = 1)
A <- build_a_matrix(ped)
dim(A)
}
