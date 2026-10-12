#ifndef POSEIDON_SD_H
#define POSEIDON_SD_H

#include <stdbool.h>
#include "f251.h"

#define POSEIDON3_ROUNDS_SD 91
#define FULL_ROUND_SD 0xF
#define PARTIAL_ROUND_SD 0x1

/* Original C field, round constants, and MDS matrix. Initialize once and reuse.
 * Treat initialized fields as read-only. Field values use Montgomery form.
 * Configuration is public; this API supplies no randomness.
 */
typedef struct {
    felt_t multipliers[3];
    felt_t round_constants[POSEIDON3_ROUNDS_SD][3];
    felt_t mds_full[3][3];
    felt_t mds_partial[3][3];
    bool permute;
    uint8_t exponent;
    uint8_t partial_index;
} poseidon3_ctx_sd;

/* Uses Python's fifth power and first-element partial S-box, with C's field,
 * constants, 4/83/4 round schedule and base MDS. Multipliers are ordinary
 * integers; Python's example values are {{1,0,0,0},{1,0,0,0},{2,0,0,0}}.
 * permute=true rotates [s0,s1,s2] to [s2,s0,s1].
 * Returns -1 for any multiplier zero modulo p, leaving ctx unchanged; 0 on
 * success. All pointers must be non-null. The C original 91x3 constants are
 * scaled and permuted internally. x^5 is not invertible in this field;
 * this is an experimental construction, not the original Poseidon permutation.
 */
int poseidon3_init_sd(poseidon3_ctx_sd *ctx, const felt_t multipliers[3],
                      bool permute);
/* Experimental variant: exponent is 3 or 5; partial_index is 0..2 BEFORE
 * permutation. (3,2) permits comparison with the old C permutation; (5,0)
 * is the default. This does NOT reproduce Python's BLS12-381 hash.
 * Invalid parameters leave ctx unchanged and return -1.
 */
int poseidon3_init_custom_sd(poseidon3_ctx_sd *ctx, const felt_t multipliers[3],
                             bool permute, uint8_t exponent, uint8_t partial_index);

/* Ordinary state -> scaled/permuted Montgomery state. Apply exactly once. */
void prepare_state_3_sd(felt_t state[3], const poseidon3_ctx_sd *ctx);
/* Montgomery state -> ordinary unscaled/unpermuted state. */
void restore_state_3_sd(felt_t state[3], const poseidon3_ctx_sd *ctx);
/* Montgomery states. Return -1 for invalid mode/index, leaving state intact;
 * 0 on success. round_idx is a round number 0..90, not a constant offset.
 */
int mix_layer_3_sd(felt_t state[3], uint8_t round_mode,
                   const poseidon3_ctx_sd *ctx);
int round_3_sd(felt_t state[3], int round_idx, uint8_t round_mode,
               const poseidon3_ctx_sd *ctx);
/* Ordinary input/output. Includes scaling/permutation and 4/83/4 rounds.
 * Output remains scaled/permuted as in the uploaded Python.
 */
void permutation_3_sd(felt_t state[3], const poseidon3_ctx_sd *ctx);
/* Same, but undo scaling/permutation. With custom (3,2) S-box settings this
 * matches permutation_3, for every valid multiplier/permute setting.
 */
void permutation_3_restored_sd(felt_t state[3], const poseidon3_ctx_sd *ctx);
/* Python's input/output convention: initialize [3,a,b], return encoded state[1].
 * Still uses C parameters, not Python's BLS12-381 parameters.
 */
void poseidon_hash_3_sd(felt_t output, const felt_t preimage[2],
                        const poseidon3_ctx_sd *ctx);

#endif
