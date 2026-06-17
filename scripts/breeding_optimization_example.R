# ============================================================================
# Evolutionary Algorithm for Optimal Contribution Selection in Breeding
# Example practice script with small population
# ============================================================================
#
# PROBLEM GEOMETRY:
#   This is a bi-objective optimization on the gain-diversity Pareto frontier.
#   - X-axis: mean kinship of next generation (MINIMIZE)  [c'Ac]
#   - Y-axis: genetic gain from selected group (MAXIMIZE) [sum(c_i * EBV_i)]
#
#   target_angle controls where on that frontier we aim:
#     0°  → pure genetic diversity: contributions spread equally across all animals
#    90°  → pure genetic gain: contributions concentrated on top 1-2 animals
#   45°  → equal weight to both objectives
#
#   The angle is converted to objective weights via sin/cos so the unit circle
#   keeps both weights on [0,1] summing to 1 at every angle.
#
# DECISION VARIABLES:
#   Each animal i gets a continuous contribution c_i >= 0.
#   Male contributions are normalized to sum to sex_ratio (0.5).
#   Female contributions are normalized to sum to 1 - sex_ratio (0.5).
#   This enforces 50/50 sex ratio exactly — no penalty term needed.
#
# OCS FORMULATION:
#   Genetic gain  = sum(c_i * EBV_i)     [weighted mean EBV]
#   Mean kinship  = c' A c               [standard OCS coancestry criterion]
#
# KEY BEHAVIOR BY ANGLE:
#   0°  → all animals get equal contribution (maximum diversity)
#  90°  → top male and top female get nearly all contribution (threshold selection)
# ============================================================================

library(GA)
library(Matrix)

utils_path <- file.path("R", "ped_utils.R")
if (!file.exists(utils_path)) {
  utils_path <- file.path("..", "R", "ped_utils.R")
}
source(utils_path)

set.seed(42)

# ============================================================================
# 1. BUILD PEDIGREE AND CANDIDATE POPULATION
# ============================================================================
# build_ped() produces a multi-generation pedigree used to compute realistic
# relationship coefficients via the A matrix.  The optimization runs on the
# last generation only (candidates for selection).

ped_full <- build_ped(
  n_male_founders   = 5,
  n_female_founders = 5,
  n_gen             = 3,
  n_males_per_gen   = 10,
  n_females_per_gen = 10,
  litter_size       = 2,
  seed              = 42
)

cat(sprintf("Full pedigree: %d animals across %d generations (0 = founders)\n\n",
  nrow(ped_full), max(ped_full$gen)))

# Candidates for selection = last generation
last_gen       <- max(ped_full$gen)
cand_df        <- ped_full[ped_full$gen == last_gen, ]
individual_ids <- cand_df$id
n_individuals  <- nrow(cand_df)

# EBVs for candidates (in practice these come from a BLUP evaluation)
population <- data.frame(
  id  = cand_df$id,
  sex = cand_df$sex,
  ebv = rnorm(n_individuals, mean = 100, sd = 15),
  row.names = cand_df$id
)

cat("=== CANDIDATE POPULATION (last generation) ===\n")
print(population)
cat("\nMales:", sum(population$sex == "M"), "| Females:", sum(population$sex == "F"), "\n\n")

# Sex indices (used throughout)
male_idx   <- which(population$sex == "M")
female_idx <- which(population$sex == "F")
n_males    <- length(male_idx)
n_females  <- length(female_idx)

# ============================================================================
# 2. BUILD A MATRIX FROM PEDIGREE
# ============================================================================

A_full   <- build_a_matrix(ped_full)
A_matrix <- A_full[individual_ids, individual_ids]

cat("=== RELATIONSHIP MATRIX — candidates (first 5x5) ===\n")
print(round(A_matrix[1:5, 1:5], 4))
cat("\n")

# ============================================================================
# 3. CONSTRAINT CONFIGURATION
# ============================================================================

constraints <- list(
  sex_ratio    = 0.5,   # 50% male, 50% female (enforced by parameterization)
  target_angle = 45     # degrees in [0, 90]
)

cat(sprintf("Target angle: %g° => gain weight = %.4f | diversity weight = %.4f\n",
  constraints$target_angle,
  sin(constraints$target_angle * pi / 180),
  cos(constraints$target_angle * pi / 180)))
cat("\n")

# ============================================================================
# 4. PARAMETERIZATION: RAW GA PARAMETERS → CONTRIBUTION VECTOR
# ============================================================================
# The GA optimizes n_males + n_females non-negative raw parameters.
# Male parameters are normalized to sum to sex_ratio; female to 1 - sex_ratio.
# This reparameterization makes the sex-ratio constraint automatically satisfied.

params_to_contributions <- function(params, male_idx, female_idx, sex_ratio) {
  n_m   <- length(male_idx)
  n_f   <- length(female_idx)
  n_tot <- n_m + n_f

  male_raw   <- params[seq_len(n_m)]
  female_raw <- params[seq_len(n_f) + n_m]

  c_vec <- numeric(n_tot)

  m_sum <- sum(male_raw)
  f_sum <- sum(female_raw)

  if (m_sum > 0) c_vec[male_idx]   <- male_raw   / m_sum * sex_ratio
  if (f_sum > 0) c_vec[female_idx] <- female_raw / f_sum * (1 - sex_ratio)

  return(c_vec)
}

# ============================================================================
# 5. EVALUATE SOLUTION
# ============================================================================
# Takes a contribution vector (not binary) and returns gain, kinship, and
# effective number of animals used (concentration metric).

evaluate_solution <- function(contributions, population, A_matrix, male_idx, female_idx) {
  c_vec <- contributions

  genetic_gain <- sum(c_vec * population$ebv)
  mean_kinship <- as.numeric(t(c_vec) %*% A_matrix %*% c_vec)

  # Effective number per sex: (sum c_i)^2 / sum(c_i^2)
  # = 1 when all weight on one animal, = n when weight is equal
  male_c   <- c_vec[male_idx]
  female_c <- c_vec[female_idx]
  eff_males   <- if (sum(male_c^2)   > 0) sum(male_c)^2   / sum(male_c^2)   else 0
  eff_females <- if (sum(female_c^2) > 0) sum(female_c)^2 / sum(female_c^2) else 0

  return(list(
    contributions = c_vec,
    genetic_gain  = genetic_gain,
    mean_kinship  = mean_kinship,
    eff_males     = eff_males,
    eff_females   = eff_females
  ))
}

# ============================================================================
# 6. NORMALIZATION REFERENCE POINTS
# ============================================================================
# Pre-compute objective extremes so both objectives normalize to [0, 1].

# Gain extremes: concentrate all weight on best (or worst) animal of each sex
best_m  <- male_idx[which.max(population$ebv[male_idx])]
best_f  <- female_idx[which.max(population$ebv[female_idx])]
worst_m <- male_idx[which.min(population$ebv[male_idx])]
worst_f <- female_idx[which.min(population$ebv[female_idx])]

max_gain <- population$ebv[best_m]  * 0.5 + population$ebv[best_f]  * 0.5
min_gain <- population$ebv[worst_m] * 0.5 + population$ebv[worst_f] * 0.5

# Kinship extremes: concentrated (best male + best female) vs uniform spread
c_conc    <- numeric(n_individuals)
c_conc[best_m] <- 0.5
c_conc[best_f] <- 0.5

c_uniform <- numeric(n_individuals)
c_uniform[male_idx]   <- (1 / n_males)   * 0.5
c_uniform[female_idx] <- (1 / n_females) * 0.5

max_kinship <- as.numeric(t(c_conc)    %*% A_matrix %*% c_conc)
min_kinship <- as.numeric(t(c_uniform) %*% A_matrix %*% c_uniform)

cat(sprintf("Gain reference:    min = %.3f | max = %.3f\n", min_gain, max_gain))
cat(sprintf("Kinship reference: min = %.6f (uniform) | max = %.6f (concentrated)\n",
  min_kinship, max_kinship))
cat("\n")

# ============================================================================
# 7. FITNESS FUNCTION
# ============================================================================
# Objectives are blended by the Pareto angle.
# Sex ratio is enforced by parameterization, not by penalty.

ga_fitness <- function(params, population, A_matrix, constraints,
                       male_idx, female_idx,
                       min_gain, max_gain, min_kinship, max_kinship) {

  c_vec <- params_to_contributions(params, male_idx, female_idx, constraints$sex_ratio)

  if (all(c_vec == 0)) return(-1e6)

  eval <- evaluate_solution(c_vec, population, A_matrix, male_idx, female_idx)

  angle_rad        <- constraints$target_angle * pi / 180
  weight_gain      <- sin(angle_rad)
  weight_diversity <- cos(angle_rad)

  gain_range    <- max_gain    - min_gain
  kinship_range <- max_kinship - min_kinship

  normalized_gain <- if (gain_range > 0)
    pmax(0, pmin(1, (eval$genetic_gain - min_gain) / gain_range)) else 0.5

  normalized_diversity <- if (kinship_range > 0)
    pmax(0, pmin(1, 1 - (eval$mean_kinship - min_kinship) / kinship_range)) else 0.5

  return(weight_gain * normalized_gain + weight_diversity * normalized_diversity)
}

# ============================================================================
# 8. RUN GENETIC ALGORITHM OPTIMIZATION
# ============================================================================

angle <- constraints$target_angle
cat("=== RUNNING GENETIC ALGORITHM OPTIMIZATION ===\n")
cat(sprintf("Target angle: %g° | Gain weight: %.4f | Diversity weight: %.4f\n",
  angle, sin(angle * pi / 180), cos(angle * pi / 180)))
cat(sprintf("Decision variables: %d males + %d females (real-valued contributions)\n\n",
  n_males, n_females))

ga_result <- ga(
  type    = "real-valued",
  fitness = function(x, ...) ga_fitness(
    x, population, A_matrix, constraints,
    male_idx, female_idx,
    min_gain, max_gain, min_kinship, max_kinship
  ),
  lower      = rep(0, n_males + n_females),
  upper      = rep(1, n_males + n_females),
  popSize    = 100,
  maxiter    = 300,
  run        = 50,
  pmutation  = 0.1,
  pcrossover = 0.8,
  seed       = 123,
  verbose    = FALSE
)

# ============================================================================
# 9. REPORT RESULTS
# ============================================================================

best_params <- ga_result@solution[1, ]
best_c      <- params_to_contributions(best_params, male_idx, female_idx, constraints$sex_ratio)
best_eval   <- evaluate_solution(best_c, population, A_matrix, male_idx, female_idx)

results_df              <- population
results_df$contribution <- round(best_c, 6)
results_df              <- results_df[order(-results_df$contribution), ]

cat("=== OPTIMIZATION RESULTS ===\n\n")
cat("CONTRIBUTION TABLE (sorted by contribution):\n")
print(results_df)

cat("\n--- SELECTION SUMMARY ---\n")
cat(sprintf("Effective # males:    %.2f / %d total males\n",   best_eval$eff_males,   n_males))
cat(sprintf("Effective # females:  %.2f / %d total females\n", best_eval$eff_females, n_females))
cat(sprintf("Male contrib total:   %.4f\n", sum(best_c[male_idx])))
cat(sprintf("Female contrib total: %.4f\n", sum(best_c[female_idx])))

cat("\n--- GENETIC GAIN ---\n")
cat(sprintf("Weighted mean EBV:  %.4f\n", best_eval$genetic_gain))
cat(sprintf("Max possible:       %.4f  (all weight on best M + best F)\n", max_gain))
cat(sprintf("Min possible:       %.4f  (all weight on worst M + worst F)\n", min_gain))

cat("\n--- KINSHIP ---\n")
cat(sprintf("Mean kinship (c'Ac): %.6f\n", best_eval$mean_kinship))
cat(sprintf("Min (uniform):       %.6f\n", min_kinship))
cat(sprintf("Max (concentrated):  %.6f\n", max_kinship))

cat("\n--- FITNESS ---\n")
cat(sprintf("Final fitness: %.6f\n", ga_result@fitnessValue))

# ============================================================================
# 10. SENSITIVITY ANALYSIS: Sweep target angle across the Pareto frontier
# ============================================================================
# Expected behavior:
#   0°  → eff_males ≈ n_males, eff_females ≈ n_females (equal contributions)
#  90°  → eff_males ≈ 1, eff_females ≈ 1 (concentrate on top animals)

cat("\n\n=== SENSITIVITY ANALYSIS: Pareto Frontier Sweep ===\n")
cat("Sweeping target angle from 0° (diversity) to 90° (gain)...\n\n")

sensitivity_results <- data.frame()

for (deg in c(0, 15, 30, 45, 60, 75, 90)) {
  constraints_test <- constraints
  constraints_test$target_angle <- deg

  ga_sens <- ga(
    type    = "real-valued",
    fitness = local({
      ct <- constraints_test
      function(x, ...) ga_fitness(
        x, population, A_matrix, ct,
        male_idx, female_idx,
        min_gain, max_gain, min_kinship, max_kinship
      )
    }),
    lower      = rep(0, n_males + n_females),
    upper      = rep(1, n_males + n_females),
    popSize    = 60,
    maxiter    = 200,
    run        = 30,
    seed       = 123,
    verbose    = FALSE
  )

  best_c_s <- params_to_contributions(
    ga_sens@solution[1, ], male_idx, female_idx, constraints_test$sex_ratio
  )
  eval_s <- evaluate_solution(best_c_s, population, A_matrix, male_idx, female_idx)

  sensitivity_results <- rbind(sensitivity_results, data.frame(
    angle_deg        = deg,
    weight_gain      = round(sin(deg * pi / 180), 3),
    weight_diversity = round(cos(deg * pi / 180), 3),
    eff_males        = round(eval_s$eff_males,     2),
    eff_females      = round(eval_s$eff_females,   2),
    genetic_gain     = round(eval_s$genetic_gain,  3),
    mean_kinship     = round(eval_s$mean_kinship,  6),
    fitness          = round(ga_sens@fitnessValue, 4)
  ))
}

print(sensitivity_results)

cat("\nInterpretation:\n")
cat("  0°  → eff_males ≈ 10, eff_females ≈ 10  (equal contributions, max diversity)\n")
cat("  90° → eff_males ≈ 1,  eff_females ≈ 1   (top animals only, max gain)\n")
cat("  Intermediate angles trace the efficient Pareto frontier.\n")
