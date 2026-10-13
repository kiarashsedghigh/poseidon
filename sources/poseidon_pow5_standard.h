#ifndef POSEIDON_POW5_STANDARD_H
#define POSEIDON_POW5_STANDARD_H

#include "f251.h"

#define FULL_ROUND_STANDARD 0xF
#define PARTIAL_ROUND_STANDARD 0x1
#define POSEIDON3_ROUNDS_STANDARD 91

/* Non-SD baseline matching the SD default parameters: C field, original
 * constants, base MDS, 4/83/4 rounds, x^5, first-element partial S-box.
 * No scaling, rotation, context, or calls to SD functions.
 * This differs from the older round_3_standard (cube, last-element S-box).
 * As in SD, x^5 is not invertible in this field; this is experimental.
 */

/* C's base MDS on a Montgomery state, in place. */
void mix_layer_3_standard(felt_t state[3]);
/* Montgomery state, round index 0..90. Returns 0 on success or -1 for invalid
 * index/mode without modifying state. Adds all three constants in every round.
 */
int round_3_pow5_standard(felt_t state[3], int round_idx, uint8_t round_mode);
/* Ordinary input/output; applies the complete 4/83/4 schedule. */
void permutation_3_pow5_standard(felt_t state[3]);

#endif
