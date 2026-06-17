# ============================================================================
# Optimal Contribution Selection with Practical Mating Constraints
# ============================================================================
#
# CONSTRAINTS:
#   Females: Exactly n_select_females = 5 selected; each mated exactly once.
#            Female selection is binary — either selected or not.
#            All selected females contribute equally: c_fi = 0.5 / 5 = 0.10
#
#   Males:   Integer matings per male in {0, 1, 2, 3}, max_male_matings = 3.
#            Total male matings = n_select_females = 5 (one per selected female).
#            Min 2 males used (3+2), max 5 males used (1+1+1+1+1).
#            Male contribution: c_mi = matings_i / 5 * 0.5 = matings_i / 10
#
# GA ENCODING (real-valued, decoded inside the fitness function):
#   First  n_males  params [0,1]: male mating weights
#     → proportional integer allocation (largest-remainder method, cap at 3)
#   Last   n_females params [0,1]: female selection scores
#     → top n_select_females females by score are selected
#
# OCS OBJECTIVES (blended by target_angle on the Pareto frontier):
#   Genetic gain: sum(c_i * EBV_i)    [maximize, weight = sin(angle)]
#   Mean kinship: c' A c               [minimize, weight = cos(angle)]
#
# EXPECTED BEHAVIOR BY ANGLE:
#   0°  → spread matings across 5 males (1+1+1+1+1), select diverse females
#  90°  → concentrate matings on top 2 males (3+2), select highest-EBV females
# ============================================================================

library(GA)
library(Matrix)
#library(pedigreemm)
library(dplyr)
library(ggplot2)

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
ped_full

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
  ebv = rnorm(n_individuals, mean = 100, sd = 10),
  row.names = cand_df$id
)

population |> arrange(sex, desc(ebv))

cat("=== CANDIDATE POPULATION (last generation) ===\n")
print(population)
cat("\nMales:", sum(population$sex == "M"), "| Females:", sum(population$sex == "F"), "\n\n")

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
  sex_ratio        = 0.5,  # 50% contribution from males, 50% from females
  target_angle     = 45,   # degrees [0, 90]: 0 = diversity only, 90 = gain only
  n_select_females = 5,    # exactly this many females selected (each mated once)
  max_male_matings = 3     # maximum number of matings per male
)

# Derived: total matings = n females selected (each female gets 1 mating)
n_total_matings <- constraints$n_select_females

cat(sprintf("n_select_females:   %d  (binary selection, 1 mating each)\n", constraints$n_select_females))
cat(sprintf("max_male_matings:   %d  per male\n",                           constraints$max_male_matings))
cat(sprintf("Total matings:      %d\n",                                      n_total_matings))
cat(sprintf("Male range:         %d to %d males used\n",
  ceiling(n_total_matings / constraints$max_male_matings),
  n_total_matings))
cat(sprintf("Target angle:       %g° | gain weight = %.4f | diversity weight = %.4f\n\n",
  constraints$target_angle,
  sin(constraints$target_angle * pi / 180),
  cos(constraints$target_angle * pi / 180)))

# ============================================================================
# 4. INTEGER MATING ALLOCATION (largest-remainder method)
# ============================================================================
# Converts continuous male weights → integer matings summing to n_total,
# each male capped at max_per. Sum is always exactly n_total.

allocate_matings <- function(weights, n_total, max_per) {
  if (sum(weights) == 0) weights <- rep(1, length(weights))

  raw   <- weights / sum(weights) * n_total
  alloc <- pmin(floor(raw), max_per)
  remain <- as.integer(n_total - sum(alloc))

  if (remain > 0) {
    # Distribute remainder by largest fractional part; skip males already at cap
    frac <- raw - floor(raw)
    frac[alloc >= max_per] <- -Inf
    for (i in order(frac, decreasing = TRUE)) {
      if (remain == 0L) break
      if (alloc[i] < max_per) {
        alloc[i] <- alloc[i] + 1L
        remain   <- remain - 1L
      }
    }
  }

  return(as.integer(alloc))
}

# ============================================================================
# 5. DECODE GA PARAMETERS → MATING PLAN
# ============================================================================
# GA has n_males + n_females real-valued parameters in [0, 1].
#   - First n_males parameters: male mating weights (decoded to integer matings)
#   - Last  n_females parameters: female selection scores (top n_select_females chosen)

params_to_mating_plan <- function(params, male_idx, female_idx, constraints) {
  n_m <- length(male_idx)
  n_f <- length(female_idx)

  male_params   <- params[seq_len(n_m)]
  female_params <- params[seq_len(n_f) + n_m]

  # Male integer mating allocation
  male_matings <- allocate_matings(
    weights = male_params,
    n_total = constraints$n_select_females,
    max_per = constraints$max_male_matings
  )

  # Female selection: top n_select_females by parameter score
  top_female_local <- order(female_params, decreasing = TRUE)[seq_len(constraints$n_select_females)]
  selected_females <- female_idx[top_female_local]  # global population row indices

  return(list(
    male_matings     = male_matings,    # integer vector, length = n_males
    selected_females = selected_females # global indices of selected females
  ))
}

# ============================================================================
# 6. MATING PLAN → CONTRIBUTION VECTOR
# ============================================================================
# Contributions sum to 1 (0.5 from males, 0.5 from females).

plan_to_contributions <- function(plan, n_individuals, male_idx, constraints) {
  c_vec <- numeric(n_individuals)

  # Male contribution proportional to matings
  c_vec[male_idx] <- plan$male_matings / constraints$n_select_females * constraints$sex_ratio

  # All selected females contribute equally
  c_vec[plan$selected_females] <- (1 / constraints$n_select_females) * (1 - constraints$sex_ratio)

  return(c_vec)
}

# ============================================================================
# 7. EVALUATE SOLUTION
# ============================================================================

evaluate_solution <- function(plan, population, A_matrix, male_idx, constraints) {
  c_vec <- plan_to_contributions(plan, nrow(population), male_idx, constraints)

  genetic_gain <- sum(c_vec * population$ebv)
  mean_kinship <- as.numeric(t(c_vec) %*% A_matrix %*% c_vec)

  n_males_active <- sum(plan$male_matings > 0)

  # Effective number of males: (sum c_m)^2 / sum(c_m^2) = n_males when equal
  male_c    <- c_vec[male_idx]
  eff_males <- if (sum(male_c^2) > 0) sum(male_c)^2 / sum(male_c^2) else 0

  return(list(
    contributions    = c_vec,
    male_matings     = plan$male_matings,
    selected_females = plan$selected_females,
    genetic_gain     = genetic_gain,
    mean_kinship     = mean_kinship,
    n_males_active   = n_males_active,
    eff_males        = eff_males
  ))
}

# ============================================================================
# 8. NORMALIZATION REFERENCE POINTS
# ============================================================================
# Compute gain and kinship for extreme feasible solutions.
# These are used to normalize both objectives to [0, 1] in the fitness function.

n_sf <- constraints$n_select_females

m_rank_desc <- order(population$ebv[male_idx],   decreasing = TRUE)
m_rank_asc  <- order(population$ebv[male_idx],   decreasing = FALSE)
f_rank_desc <- order(population$ebv[female_idx], decreasing = TRUE)
f_rank_asc  <- order(population$ebv[female_idx], decreasing = FALSE)

# Max gain: 3 matings to best male, 2 to 2nd best; top 5 females by EBV
max_matings                    <- rep(0L, n_males)
max_matings[m_rank_desc[1]]    <- 3L
max_matings[m_rank_desc[2]]    <- 2L
max_gain_plan <- list(
  male_matings     = max_matings,
  selected_females = female_idx[f_rank_desc[seq_len(n_sf)]]
)

# Min gain: 3 matings to worst male, 2 to 2nd worst; bottom 5 females by EBV
min_matings                   <- rep(0L, n_males)
min_matings[m_rank_asc[1]]    <- 3L
min_matings[m_rank_asc[2]]    <- 2L
min_gain_plan <- list(
  male_matings     = min_matings,
  selected_females = female_idx[f_rank_asc[seq_len(n_sf)]]
)

# Spread plan: 5 males × 1 mating, first 5 females — used as min kinship reference
spread_matings               <- rep(0L, n_males)
spread_matings[seq_len(n_sf)] <- 1L
spread_plan <- list(
  male_matings     = spread_matings,
  selected_females = female_idx[seq_len(n_sf)]
)

max_gain_eval <- evaluate_solution(max_gain_plan, population, A_matrix, male_idx, constraints)
min_gain_eval <- evaluate_solution(min_gain_plan, population, A_matrix, male_idx, constraints)
spread_eval   <- evaluate_solution(spread_plan,   population, A_matrix, male_idx, constraints)

ref_max_gain    <- max_gain_eval$genetic_gain
ref_min_gain    <- min_gain_eval$genetic_gain
ref_max_kinship <- max(max_gain_eval$mean_kinship, spread_eval$mean_kinship)
ref_min_kinship <- min(max_gain_eval$mean_kinship, spread_eval$mean_kinship)

cat(sprintf("Gain reference:    min = %.3f | max = %.3f\n", ref_min_gain, ref_max_gain))
cat(sprintf("Kinship reference: min = %.6f | max = %.6f\n", ref_min_kinship, ref_max_kinship))
cat("\n")

# ============================================================================
# 9. FITNESS FUNCTION
# ============================================================================

ga_fitness <- function(params, population, A_matrix, constraints,
                       male_idx, female_idx,
                       ref_min_gain, ref_max_gain,
                       ref_min_kinship, ref_max_kinship) {

  plan  <- params_to_mating_plan(params, male_idx, female_idx, constraints)
  c_vec <- plan_to_contributions(plan, nrow(population), male_idx, constraints)

  genetic_gain <- sum(c_vec * population$ebv)
  mean_kinship <- as.numeric(t(c_vec) %*% A_matrix %*% c_vec)

  angle_rad        <- constraints$target_angle * pi / 180
  weight_gain      <- sin(angle_rad)
  weight_diversity <- cos(angle_rad)

  gain_range    <- ref_max_gain    - ref_min_gain
  kinship_range <- ref_max_kinship - ref_min_kinship

  norm_gain <- if (gain_range > 0)
    pmax(0, pmin(1, (genetic_gain - ref_min_gain) / gain_range)) else 0.5

  norm_diversity <- if (kinship_range > 0)
    pmax(0, pmin(1, 1 - (mean_kinship - ref_min_kinship) / kinship_range)) else 0.5

  return(weight_gain * norm_gain + weight_diversity * norm_diversity)
}

# ============================================================================
# 10. RUN OPTIMIZATION
# ============================================================================

angle <- constraints$target_angle

cat("=== RUNNING GENETIC ALGORITHM OPTIMIZATION ===\n")
cat(sprintf("Target angle: %g° | Gain weight: %.4f | Diversity weight: %.4f\n",
  angle, sin(angle * pi / 180), cos(angle * pi / 180)))
cat(sprintf("GA parameters: %d male weights + %d female scores = %d total\n\n",
  n_males, n_females, n_males + n_females))

ga_result <- ga(
  type    = "real-valued",
  fitness = function(x, ...) ga_fitness(
    x, population, A_matrix, constraints,
    male_idx, female_idx,
    ref_min_gain, ref_max_gain,
    ref_min_kinship, ref_max_kinship
  ),
  lower      = rep(0, n_males + n_females),
  upper      = rep(1, n_males + n_females),
  popSize    = 200,
  maxiter    = 500,
  run        = 80,
  pmutation  = 0.15,
  pcrossover = 0.8,
  seed       = 123,
  verbose    = FALSE
)

# ============================================================================
# 11. REPORT RESULTS
# ============================================================================

best_params <- ga_result@solution[1, ]
best_plan   <- params_to_mating_plan(best_params, male_idx, female_idx, constraints)
best_eval   <- evaluate_solution(best_plan, population, A_matrix, male_idx, constraints)

# Male table: all males, showing mating count
male_df              <- population[male_idx, c("id", "sex", "ebv")]
male_df$matings      <- best_plan$male_matings
male_df$contribution <- round(best_eval$contributions[male_idx], 4)
male_df              <- male_df[order(-male_df$matings, -male_df$ebv), ]

# Female table: all females, marking selected ones
female_df              <- population[female_idx, c("id", "sex", "ebv")]
female_df$selected     <- ifelse(female_idx %in% best_plan$selected_females, "YES", "no")
female_df$contribution <- round(best_eval$contributions[female_idx], 4)
female_df              <- female_df[order(female_df$selected == "no", -female_df$ebv), ]

cat("=== OPTIMIZATION RESULTS ===\n\n")

cat("MALE MATING TABLE:\n")
print(male_df)

cat("\nFEMALE SELECTION TABLE:\n")
print(female_df)

active_matings <- sort(best_plan$male_matings[best_plan$male_matings > 0], decreasing = TRUE)

cat("\n--- MATING SUMMARY ---\n")
cat(sprintf("Males active:       %d  (mating split: %s)\n",
  best_eval$n_males_active,
  paste(active_matings, collapse = "+")))
cat(sprintf("Effective # males:  %.2f  (1 = fully concentrated, %d = fully spread)\n",
  best_eval$eff_males, n_total_matings))
cat(sprintf("Females selected:   %d / %d\n", length(best_plan$selected_females), n_females))
cat(sprintf("Total matings:      %d\n", sum(best_plan$male_matings)))

cat("\n--- GENETIC GAIN ---\n")
cat(sprintf("Weighted mean EBV:  %.4f\n", best_eval$genetic_gain))
cat(sprintf("Max feasible:       %.4f  (3+2 on top males, top 5 females)\n", ref_max_gain))
cat(sprintf("Min feasible:       %.4f  (3+2 on worst males, bottom 5 females)\n", ref_min_gain))

cat("\n--- KINSHIP ---\n")
cat(sprintf("Mean kinship (c'Ac): %.6f\n", best_eval$mean_kinship))
cat(sprintf("Min reference:       %.6f\n", ref_min_kinship))
cat(sprintf("Max reference:       %.6f\n", ref_max_kinship))

cat("\n--- FITNESS ---\n")
cat(sprintf("Final fitness: %.6f\n", ga_result@fitnessValue))

# ============================================================================
# 12. SENSITIVITY ANALYSIS: Sweep the Pareto frontier
# ============================================================================
# Expected pattern:
#   0°  → eff_males ≈ 5 (mating split 1+1+1+1+1), females chosen for diversity
#  90°  → eff_males ≈ 1.9 (mating split 3+2), females chosen for gain

cat("\n\n=== SENSITIVITY ANALYSIS: Pareto Frontier Sweep ===\n")
cat("Sweeping target angle from 0° (diversity) to 90° (gain)...\n\n")

sensitivity_results <- data.frame()

for (deg in c(0, 15, 30, 45, 60, 75, 90)) {
  constraints_test              <- constraints
  constraints_test$target_angle <- deg

  ga_sens <- ga(
    type    = "real-valued",
    fitness = local({
      ct <- constraints_test
      function(x, ...) ga_fitness(
        x, population, A_matrix, ct,
        male_idx, female_idx,
        ref_min_gain, ref_max_gain,
        ref_min_kinship, ref_max_kinship
      )
    }),
    lower      = rep(0, n_males + n_females),
    upper      = rep(1, n_males + n_females),
    popSize    = 150,
    maxiter    = 300,
    run        = 50,
    pmutation  = 0.15,
    pcrossover = 0.8,
    seed       = 123,
    verbose    = FALSE
  )

  plan_s <- params_to_mating_plan(ga_sens@solution[1, ], male_idx, female_idx, constraints_test)
  eval_s <- evaluate_solution(plan_s, population, A_matrix, male_idx, constraints_test)

  active_s    <- sort(plan_s$male_matings[plan_s$male_matings > 0], decreasing = TRUE)
  mating_str  <- paste(active_s, collapse = "+")

  sensitivity_results <- rbind(sensitivity_results, data.frame(
    angle_deg    = deg,
    weight_gain  = round(sin(deg * pi / 180), 3),
    weight_div   = round(cos(deg * pi / 180), 3),
    males_active = eval_s$n_males_active,
    eff_males    = round(eval_s$eff_males, 2),
    mating_split = mating_str,
    genetic_gain = round(eval_s$genetic_gain, 3),
    mean_kinship = round(eval_s$mean_kinship, 6),
    fitness      = round(ga_sens@fitnessValue, 4),
    stringsAsFactors = FALSE
  ))
}

print(sensitivity_results)

cat("\nInterpretation:\n")
cat("  0°  → matings spread: eff_males near 5, mating_split '1+1+1+1+1', lower kinship\n")
cat("  90° → matings concentrated: eff_males near 2, mating_split '3+2', highest gain\n")
cat("  eff_males = 1.92 for 3+2, 5.0 for 1+1+1+1+1 (theoretical bounds)\n")
