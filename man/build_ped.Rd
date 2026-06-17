\name{build_ped}
\alias{build_ped}
\title{Simulate a multi-generation pedigree}
\usage{
build_ped(
  n_male_founders = 5,
  n_female_founders = 5,
  n_gen = 3,
  n_males_per_gen = 10,
  n_females_per_gen = 10,
  litter_size = 2,
  seed = NULL
)
}
\arguments{
\item{n_male_founders}{Males in generation 0 (no parents)}

\item{n_female_founders}{Females in generation 0}

\item{n_gen}{Number of offspring generations to simulate}

\item{n_males_per_gen}{Target male offspring per generation}

\item{n_females_per_gen}{Target female offspring per generation}

\item{litter_size}{Offspring per dam per generation}

\item{seed}{Optional RNG seed for reproducibility}
}
\value{
A \code{data.frame} with columns \code{id}, \code{sire}, \code{dam}, \code{sex},
and \code{gen}. Rows are ordered parents before offspring, which is required
for additive relationship matrix construction.
}
\description{
Founders in generation 0 have no parents. Each subsequent generation uses
random mating: dams are sampled without replacement where possible (each dam
mated at most once per generation); sires are sampled with replacement.
}
\examples{
ped <- build_ped(n_gen = 1, seed = 1)
head(ped)
}
