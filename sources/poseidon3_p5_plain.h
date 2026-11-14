#ifndef POSEIDON3_P5_PLAIN_H
#define POSEIDON3_P5_PLAIN_H

#include "f251.h"

#define FULL_ROUND_STANDARD 0xF
#define PARTIAL_ROUND_STANDARD 0x1
#define POSEIDON3_ROUNDS_STANDARD 91

/* Non-SD baseline matching the SD default parameters: C field, original
 * constants, base MDS, 4/83/4 rounds, x^5, first-element partial S-box.
 * No scaling, rotation, or calls to SD functions.
 * This differs from the older round_3_standard (cube, last-element S-box).
 * As in SD, x^5 is not invertible in this field; this is experimental.
 */

/* Initialize once before rounds; treat the context as read-only.
 * modulus is ordinary; constants and matrices are in Montgomery form.
 */
typedef struct {
    felt_t modulus;
    felt_t round_constants[POSEIDON3_ROUNDS_STANDARD][3];
    felt_t mds_full[3][3];
    felt_t mds_partial[3][3];
} poseidon3_ctx_standard;
void poseidon3_init_standard(poseidon3_ctx_standard *ctx);
/* Same general matrix loop as SD. Invalid mode returns -1 without mutation. */
int mix_layer_3_standard(felt_t state[3], uint8_t round_mode,
                         const poseidon3_ctx_standard *ctx);
/* Montgomery state, round index 0..90. Returns 0 on success or -1 for invalid
 * index/mode without modifying state. Adds all three constants in every round.
 */
int round_3_pow5_standard(felt_t state[3], int round_idx, uint8_t round_mode,
                          const poseidon3_ctx_standard *ctx);
/* Ordinary input/output; applies the complete 4/83/4 schedule. */
void permutation_3_pow5_standard(felt_t state[3], const poseidon3_ctx_standard *ctx);

#endif
