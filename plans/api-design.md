# evobreed API Design Notes

## Goal

Move the optimization logic currently embedded in `scripts/` into a small,
stable R package API. Users should be able to provide:

- a candidate table with animal IDs, sex, and EBV/index values
- either a pedigree or an already-built relationship matrix
- optimization options, such as target gain/diversity angle and mating constraints

The package should then return an object with contributions, selection/mating
details, objective values, and enough metadata to inspect or plot the solution.

This package does not need a large database-like editing API. The core user
workflow is closer to "prepare inputs -> solve OCS -> summarize/plot results".

## Suggested Public API

### Input preparation

```r
ocs_data <- prepare_ocs_data(
  candidates,
  id_col = "id",
  sex_col = "sex",
  value_col = "ebv",
  pedigree = NULL,
  relationship_matrix = NULL,
  candidate_ids = NULL
)
```

Purpose:

- standardize user-provided columns into `id`, `sex`, and `value`
- verify sex labels, duplicate IDs, missing EBVs, and matrix dimensions
- build and subset `A` when `pedigree` is supplied
- accept a precomputed relationship matrix when users already have `A` or `G`

Return:

An `evobreed_ocs_data` object, probably a simple list:

```r
list(
  candidates = data.frame(id, sex, value, ...),
  A = matrix,
  male_idx = integer(),
  female_idx = integer(),
  value_name = "ebv"
)
```

Notes:

- Keep column remapping here so solver functions can assume canonical columns.
- `relationship_matrix` should be allowed to be an additive relationship matrix,
  genomic relationship matrix, or any symmetric matrix on the same scale.
- `candidate_ids` is useful when the relationship matrix contains more animals
  than the candidate set.

### Continuous OCS solver

```r
solve_ocs(
  candidates,
  pedigree = NULL,
  relationship_matrix = NULL,
  target_angle = 45,
  sex_ratio = 0.5,
  value_col = "ebv",
  id_col = "id",
  sex_col = "sex",
  ga_control = ocs_ga_control(),
  seed = NULL
)
```

Purpose:

- solve the continuous contribution problem from
  `scripts/breeding_optimization_example.R`
- enforce male and female contribution totals by parameterization
- blend normalized gain and diversity using `sin(angle)` and `cos(angle)`

Return:

An `evobreed_ocs_result` object:

```r
list(
  mode = "continuous",
  call = match.call(),
  candidates = data.frame(...),
  contributions = named numeric,
  solution = data.frame(id, sex, value, contribution),
  objectives = list(
    genetic_gain = numeric,
    mean_kinship = numeric,
    fitness = numeric,
    weight_gain = numeric,
    weight_diversity = numeric
  ),
  diagnostics = list(
    effective_males = numeric,
    effective_females = numeric,
    male_contribution = numeric,
    female_contribution = numeric
  ),
  references = list(
    min_gain = numeric,
    max_gain = numeric,
    min_kinship = numeric,
    max_kinship = numeric
  ),
  ga = ga_result
)
```

Example:

```r
fit <- solve_ocs(
  candidates = population,
  pedigree = ped_full,
  target_angle = 45,
  value_col = "ebv",
  seed = 123
)

summary(fit)
contribution_table(fit)
```

### Constrained mating-plan solver

```r
solve_mating_plan(
  candidates,
  pedigree = NULL,
  relationship_matrix = NULL,
  target_angle = 45,
  sex_ratio = 0.5,
  n_select_females,
  max_male_matings,
  value_col = "ebv",
  id_col = "id",
  sex_col = "sex",
  ga_control = ocs_ga_control(pop_size = 200, maxiter = 500, run = 80),
  seed = NULL
)
```

Purpose:

- wrap the logic from `scripts/breeding_optimization_constrained.R`
- choose exactly `n_select_females`
- allocate integer male matings summing to `n_select_females`
- cap each male at `max_male_matings`
- return both contribution weights and practical mating counts

Return:

Another `evobreed_ocs_result`, with `mode = "mating_plan"` and extra fields:

```r
list(
  male_matings = named integer,
  selected_females = character(),
  male_table = data.frame(id, sex, value, matings, contribution),
  female_table = data.frame(id, sex, value, selected, contribution)
)
```

Example:

```r
plan <- solve_mating_plan(
  candidates = population,
  pedigree = ped_full,
  target_angle = 45,
  n_select_females = 5,
  max_male_matings = 3,
  seed = 123
)

mating_table(plan)
contribution_table(plan)
```

### Pareto/frontier sweep

```r
sweep_ocs(
  candidates,
  pedigree = NULL,
  relationship_matrix = NULL,
  angles = seq(0, 90, by = 15),
  mode = c("continuous", "mating_plan"),
  ...,
  seed = NULL
)
```

Purpose:

- turn the sensitivity loops in both scripts into a reusable function
- run `solve_ocs()` or `solve_mating_plan()` for each angle
- return compact metrics for plotting and the full fit objects for inspection

Return:

```r
list(
  results = data.frame(
    angle_deg,
    weight_gain,
    weight_diversity,
    genetic_gain,
    mean_kinship,
    fitness,
    effective_males,
    effective_females
  ),
  fits = list()
)
```

For mating plans, include `males_active` and `mating_split`.

### Control helper

```r
ocs_ga_control <- function(
  pop_size = 100,
  maxiter = 300,
  run = 50,
  pmutation = 0.1,
  pcrossover = 0.8,
  verbose = FALSE
) {
  list(
    popSize = pop_size,
    maxiter = maxiter,
    run = run,
    pmutation = pmutation,
    pcrossover = pcrossover,
    verbose = verbose
  )
}
```

Purpose:

- avoid exposing `GA::ga()` details in every solver signature
- keep defaults visible and documented
- allow advanced users to tune without making the high-level functions noisy

### Result helpers

```r
contribution_table <- function(x, sort = TRUE, digits = 6) {}
mating_table <- function(x, sort = TRUE) {}
objective_values <- function(x) {}
effective_number <- function(contributions) {}
```

S3 methods:

```r
print.evobreed_ocs_result <- function(x, ...) {}
summary.evobreed_ocs_result <- function(object, ...) {}
plot.evobreed_ocs_sweep <- function(x, ...) {}
```

Keep these helpers lightweight. The solver result should already contain the
data needed for users to make their own tables with base R, dplyr, or ggplot2.

## Internal Functions to Extract From Scripts

These should be internal helpers at first. Export later only if there is a clear
user need.

```r
validate_ocs_inputs()
subset_relationship_matrix()
angle_weights()
normalize_by_sex()
params_to_contributions()
evaluate_contributions()
continuous_reference_points()
continuous_fitness()
allocate_matings()
params_to_mating_plan()
plan_to_contributions()
evaluate_mating_plan()
mating_reference_points()
mating_plan_fitness()
run_ga()
```

Suggested file split:

- `R/input.R`: input validation and standardization
- `R/objectives.R`: gain, kinship, effective number, angle weights
- `R/continuous_ocs.R`: continuous solver and helpers
- `R/mating_plan.R`: constrained mating solver and helpers
- `R/sweep.R`: frontier sweeps
- `R/results.R`: S3 print/summary/table methods

## Wrapper Shape From Current Scripts

The scripts can become examples once the solvers exist.

Current script flow:

```r
ped_full <- build_ped(...)
candidates <- last_generation_with_ebv(ped_full)
A_full <- build_a_matrix(ped_full)
A_candidates <- A_full[candidates$id, candidates$id]
# many local helper functions
# GA call
# reporting tables
```

Package workflow:

```r
ped_full <- build_ped(...)
candidates <- subset(ped_full, gen == max(gen))
candidates$ebv <- rnorm(nrow(candidates), mean = 100, sd = 15)

fit <- solve_ocs(
  candidates = candidates,
  pedigree = ped_full,
  target_angle = 45,
  value_col = "ebv",
  seed = 123
)

summary(fit)
contribution_table(fit)
```

Constrained workflow:

```r
plan <- solve_mating_plan(
  candidates = candidates,
  pedigree = ped_full,
  target_angle = 45,
  n_select_females = 5,
  max_male_matings = 3,
  value_col = "ebv",
  seed = 123
)

summary(plan)
mating_table(plan)
```

Frontier workflow:

```r
frontier <- sweep_ocs(
  candidates = candidates,
  pedigree = ped_full,
  angles = c(0, 15, 30, 45, 60, 75, 90),
  mode = "continuous",
  value_col = "ebv",
  seed = 123
)

frontier$results
plot(frontier)
```

## API Design Principles

- Accept either `pedigree` or `relationship_matrix`, but require exactly one.
- Keep `build_ped()` and `build_a_matrix()` as utility functions, not required
  steps for every user.
- Treat EBV, index, or merit values generically as `value`; expose `value_col`
  for user data and display `value_name` in output.
- Return rich objects, not printed reports. Printing is a presentation layer.
- Keep defaults friendly for examples, but expose `ga_control` for real use.
- Do not export low-level GA encoders until users ask for them.
- Avoid penalty terms for constraints that can be enforced by parameterization.
- Keep all row order dependencies internal; users should interact by IDs.

## Validation Rules

Input validation should be strict because optimization failures can be opaque.

Candidate table:

- has unique, non-missing IDs
- has at least one male and one female
- has numeric, finite EBV/index values
- sex labels are either normalized to `M`/`F` or rejected with a clear message

Relationship matrix:

- square, numeric, finite
- row and column names present
- contains all candidate IDs
- symmetric within tolerance
- subset/reordered to match candidate order

Continuous OCS:

- `target_angle` is numeric in `[0, 90]`
- `sex_ratio` is numeric in `(0, 1)`

Mating plan:

- `n_select_females <= number of candidate females`
- `max_male_matings >= 1`
- male capacity is enough:
  `number_of_males * max_male_matings >= n_select_females`

## Tests to Add First

Core helpers:

- `angle_weights(0)`, `angle_weights(45)`, `angle_weights(90)`
- `params_to_contributions()` sums to 1 and enforces sex totals
- `effective_number()` returns 1 for concentrated and n for equal weights
- `allocate_matings()` always sums to `n_total` and respects `max_per`
- `plan_to_contributions()` sums to 1 and gives selected females equal weights

Solver smoke tests:

- `solve_ocs()` returns class `evobreed_ocs_result`
- continuous result contributions sum to 1
- male/female contribution totals match `sex_ratio`
- `solve_mating_plan()` returns exact number of selected females
- mating-plan male matings sum to `n_select_females`
- frontier sweep returns one row per angle

Use small populations and low GA iteration settings in tests to keep runtime
acceptable. Tests should verify invariants more than exact stochastic optima.

## Open Questions

1. Should the main value column be named `ebv` everywhere, or should the public
   API use a generic `value_col` because users may optimize on index values?

2. Should constrained female selection always mean "selected females contribute
   equally", or should the package eventually support unequal female
   contributions with min/max bounds?

3. Should `relationship_matrix` be explicitly called `A`, or should the package
   use a neutral name to support genomic relationship matrices?

4. Is the Pareto angle the preferred public control, or should the package also
   support direct objective weights like `gain_weight` and `diversity_weight`?

5. Should the primary OCS solver support a hard maximum mean-kinship constraint
   in addition to the current weighted objective?

6. What should be the default GA settings for realistic package use versus fast
   examples and tests?

7. Should `GA` move from `Suggests` to `Imports` once solver functions are
   exported, or should solvers check for `GA` at runtime and give an install
   message?

8. Should plotting helpers depend on `ggplot2`, or should the package return
   tidy data and leave plotting to users?

9. Should mating-plan output include explicit sire-dam pair assignments, or is
   male mating count plus selected females enough for the first version?

10. Should the package support single-sex or clonal contribution problems later,
    or should the first API assume two-sex livestock breeding throughout?

## Suggested Implementation Order

1. Add internal objective helpers and tests.
2. Add input preparation and relationship-matrix subsetting.
3. Implement `solve_ocs()` by extracting the continuous script logic.
4. Add result classes plus `print()`, `summary()`, and `contribution_table()`.
5. Implement `solve_mating_plan()` by extracting the constrained script logic.
6. Add `sweep_ocs()`.
7. Convert scripts into short examples that call package functions.
8. Update README quick start to use exported solvers instead of `source()`.

