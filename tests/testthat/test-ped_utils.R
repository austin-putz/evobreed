test_that("build_ped returns founders and offspring in parent-before-child order", {
  ped <- build_ped(
    n_male_founders = 2,
    n_female_founders = 3,
    n_gen = 2,
    n_males_per_gen = 2,
    n_females_per_gen = 2,
    litter_size = 2,
    seed = 1
  )

  expect_s3_class(ped, "data.frame")
  expect_named(ped, c("id", "sire", "dam", "sex", "gen"))
  expect_equal(nrow(ped), 2 + 3 + 2 * (2 + 2))
  expect_true(all(is.na(ped$sire[ped$gen == 0])))
  expect_true(all(is.na(ped$dam[ped$gen == 0])))
  expect_true(all(ped$sex %in% c("M", "F")))

  id_pos <- setNames(seq_len(nrow(ped)), ped$id)
  offspring <- ped$gen > 0
  expect_true(all(id_pos[ped$sire[offspring]] < which(offspring)))
  expect_true(all(id_pos[ped$dam[offspring]] < which(offspring)))
})

test_that("build_a_matrix returns a symmetric named relationship matrix", {
  ped <- build_ped(
    n_male_founders = 2,
    n_female_founders = 2,
    n_gen = 1,
    n_males_per_gen = 2,
    n_females_per_gen = 2,
    litter_size = 2,
    seed = 2
  )

  A <- build_a_matrix(ped)

  expect_type(A, "double")
  expect_equal(dim(A), c(nrow(ped), nrow(ped)))
  expect_equal(rownames(A), ped$id)
  expect_equal(colnames(A), ped$id)
  expect_equal(A, t(A), tolerance = 1e-12)
  expect_true(all(diag(A) >= 1))
  expect_equal(unname(diag(A[ped$gen == 0, ped$gen == 0, drop = FALSE])), rep(1, 4))
})
