/**
 * @file dllifting_c.h
 * @brief C ABI wrapper for DLLifting (usable from C and solver callbacks).
 */
#ifndef DLLIFTING_C_H
#define DLLIFTING_C_H

#ifdef __cplusplus
extern "C" {
#endif

/** Return codes for dllifting_lift_cover */
#define DLLIFTING_OK           0
#define DLLIFTING_ERR_ALLOC   -1
#define DLLIFTING_ERR_ARGS    -2

/**
 * Lift a cover inequality for a single knapsack row.
 *
 * Modifies @p coef in place and writes the lifted inequality rhs to @p rhs.
 *
 * @param n              Number of variables in the row
 * @param coef           Coefficients (input seed values, output lifted)
 * @param weight         Knapsack weights w_i
 * @param ub             Upper bounds u_i (use a large value for unbounded)
 * @param use_ub         1 if variable fixed at UB for lifting, else 0
 * @param cap            Knapsack capacity / demand b
 * @param is_subcap      1 if cap is residual subcapacity, 0 if global rhs
 * @param seed           Indices in the seed cover (length n_seed)
 * @param n_seed         Number of seed indices
 * @param lifting_order  Remaining variables in lifting order
 * @param n_order        Length of lifting_order
 * @param rhs            Output: lifted inequality rhs
 * @param is_leq         1 for sum w x <= b, 0 for sum w x >= b
 * @param threshold      < 100 prefers DL merges, > 100 prefers DP
 * @param x_frac         Optional fractional point (may be NULL)
 * @return DLLIFTING_OK on success
 */
int dllifting_lift_cover(
    int n,
    double* coef,
    const double* weight,
    const double* ub,
    const int* use_ub,
    double cap,
    int is_subcap,
    const int* seed,
    int n_seed,
    const int* lifting_order,
    int n_order,
    double* rhs,
    int is_leq,
    double threshold,
    const double* x_frac);

#ifdef __cplusplus
}
#endif

#endif /* DLLIFTING_C_H */
