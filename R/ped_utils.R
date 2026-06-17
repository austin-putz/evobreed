#' Simulate a multi-generation pedigree
#'
#' Founders in generation 0 have no parents. Each subsequent generation uses
#' random mating: dams are sampled without replacement where possible (each
#' dam mated at most once per generation); sires are sampled with replacement.
#'
#' @param n_male_founders Males in generation 0 (no parents)
#' @param n_female_founders Females in generation 0
#' @param n_gen Number of offspring generations to simulate
#' @param n_males_per_gen Target male offspring per generation
#' @param n_females_per_gen Target female offspring per generation
#' @param litter_size Offspring per dam per generation
#' @param seed Optional RNG seed for reproducibility
#'
#' @return A `data.frame` with columns `id`, `sire`, `dam`, `sex`, and `gen`.
#' Rows are ordered parents before offspring, which is required for additive
#' relationship matrix construction.
#' @export
#'
#' @examples
#' ped <- build_ped(n_gen = 1, seed = 1)
#' head(ped)
build_ped <- function(
    n_male_founders   = 5,
    n_female_founders = 5,
    n_gen             = 3,
    n_males_per_gen   = 10,
    n_females_per_gen = 10,
    litter_size       = 2,
    seed              = NULL
) {
  if (!is.null(seed)) set.seed(seed)

  ctr <- 0L
  nxt <- function(n = 1L) {
    ids <- sprintf("IND_%04d", ctr + seq_len(n))
    ctr <<- ctr + n
    ids
  }

  rows <- list()

  # --- Generation 0: founders (no parents) ---
  m0 <- nxt(n_male_founders)
  f0 <- nxt(n_female_founders)
  rows[[1]] <- data.frame(
    id   = c(m0, f0),
    sire = NA_character_,
    dam  = NA_character_,
    sex  = c(rep("M", n_male_founders), rep("F", n_female_founders)),
    gen  = 0L,
    stringsAsFactors = FALSE
  )
  sire_pool <- m0
  dam_pool  <- f0

  # --- Subsequent generations ---
  for (gen in seq_len(n_gen)) {
    n_tot   <- n_males_per_gen + n_females_per_gen
    n_lit   <- ceiling(n_tot / litter_size)
    m_dams  <- sample(dam_pool,  n_lit, replace = n_lit > length(dam_pool))
    m_sires <- sample(sire_pool, n_lit, replace = TRUE)

    gen_rows <- list()
    new_m    <- character(0)
    new_f    <- character(0)
    done     <- FALSE

    for (k in seq_len(n_lit)) {
      for (l in seq_len(litter_size)) {
        nm <- n_males_per_gen   - length(new_m)
        nf <- n_females_per_gen - length(new_f)
        if (nm <= 0 && nf <= 0) { done <- TRUE; break }

        sex <- if (nm <= 0) "F" else if (nf <= 0) "M" else sample(c("M","F"), 1)
        id  <- nxt(1L)

        gen_rows[[length(gen_rows) + 1]] <- list(
          id = id, sire = m_sires[k], dam = m_dams[k], sex = sex
        )
        if (sex == "M") new_m <- c(new_m, id) else new_f <- c(new_f, id)
      }
      if (done) break
    }

    rows[[length(rows) + 1]] <- data.frame(
      id   = vapply(gen_rows, `[[`, character(1), "id"),
      sire = vapply(gen_rows, `[[`, character(1), "sire"),
      dam  = vapply(gen_rows, `[[`, character(1), "dam"),
      sex  = vapply(gen_rows, `[[`, character(1), "sex"),
      gen  = as.integer(gen),
      stringsAsFactors = FALSE
    )
    sire_pool <- new_m
    dam_pool  <- new_f
  }

  ped_df <- do.call(rbind, rows)
  rownames(ped_df) <- NULL
  return(ped_df)
}


#' Compute the additive relationship matrix (A) from a pedigree
#'
#' A[i,j] = 2 * kinship(i,j).  Diagonal = 1 + F_i (inbreeding coefficient).
#' Founders are assumed non-inbred and unrelated to each other.
#'
#' Uses pedigreemm::getA() when the package is installed; otherwise falls back
#' to the Henderson (1976) tabular method (O(n^2), suitable for n < ~5000).
#'
#' @param ped_df data.frame with columns id, sire, dam (as returned by build_ped)
#' @return Named n-by-n matrix A
#' @export
#'
#' @examples
#' ped <- build_ped(n_gen = 1, seed = 1)
#' A <- build_a_matrix(ped)
#' dim(A)
build_a_matrix <- function(ped_df) {
  n   <- nrow(ped_df)
  ids <- as.character(ped_df$id)

  if (requireNamespace("pedigreemm", quietly = TRUE)) {
    ped <- pedigreemm::pedigree(
      sire  = as.character(ped_df$sire),
      dam   = as.character(ped_df$dam),
      label = ids
    )
    A <- as.matrix(pedigreemm::getA(ped))
    rownames(A) <- colnames(A) <- ids
    return(A)
  }

  message("pedigreemm not installed - using Henderson (1976) tabular method")
  message("  Install with: install.packages('pedigreemm')")

  # Tabular method:
  #   A[i,i] = 1 + 0.5 * A[sire, dam]
  #   A[i,j] = 0.5 * (A[sire_i, j] + A[dam_i, j])  for j < i
  lookup <- stats::setNames(seq_len(n), ids)
  si <- lookup[as.character(ped_df$sire)]  # NA for founders
  di <- lookup[as.character(ped_df$dam)]

  A <- diag(1.0, n)

  for (i in seq_len(n)) {
    s <- si[[i]]
    d <- di[[i]]
    if (!is.na(s) && !is.na(d)) A[i, i] <- 1.0 + A[s, d] * 0.5
    for (j in seq_len(i - 1L)) {
      aij <- 0.0
      if (!is.na(s)) aij <- aij + A[s, j] * 0.5
      if (!is.na(d)) aij <- aij + A[d, j] * 0.5
      A[i, j] <- aij
      A[j, i] <- aij
    }
  }

  rownames(A) <- colnames(A) <- ids
  return(A)
}
