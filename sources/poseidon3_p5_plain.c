#include <string.h>
#include "poseidon3_p5_plain.h"
#include "poseidon_rc.h"

static void reduce_standard(felt_t z, const felt_t x)
{
    const uint64_t p3 = F251_MODULUS[3];
    uint64_t q = x[3] >> 59;
    felt_t sub = {q, 0, 0, q*p3}, tmp;
    uint64_t borrow = 0, carry = 0;
    for (int i = 0; i < 4; i++) {
        uint64_t a = x[i] - sub[i];
        uint64_t b = x[i] < sub[i];
        tmp[i] = a - borrow;
        borrow = b | (a < borrow);
    }
    uint64_t mask = (uint64_t)0 - borrow;
    for (int i = 0; i < 4; i++) {
        uint64_t a = tmp[i] + (F251_MODULUS[i] & mask);
        uint64_t c = a < tmp[i];
        z[i] = a + carry;
        carry = c | (z[i] < a);
    }
}

void poseidon3_init_standard(poseidon3_ctx_standard *ctx)
{
    const int matrix[3][3] = {{3, 1, 1}, {1, -1, 1}, {1, 1, -2}};
    memcpy(ctx->modulus, F251_MODULUS, sizeof(ctx->modulus));
    memcpy(ctx->round_constants, ROUND_CONST_MONT_P3, sizeof(ctx->round_constants));
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            int entry = matrix[i][j];
            felt_t base = {(uint64_t)(entry < 0 ? -entry : entry), 0, 0, 0};
            if (entry < 0) {
                f251_sub(base, ctx->modulus, base);
                reduce_standard(base, base);
            }
            f251_to_montgomery(base, base);
            reduce_standard(ctx->mds_full[i][j], base);
            memcpy(ctx->mds_partial[i][j], ctx->mds_full[i][j], sizeof(felt_t));
        }
    }
}

int mix_layer_3_standard(felt_t state[3], uint8_t round_mode, const poseidon3_ctx_standard *ctx)
{
    /* Previous optimized body, retained for reference:
    // M = [[3,1,1], [1,-1,1], [1,1,-2]].
    felt_t sum, twice, thrice;
    f251_sum_state_3(sum, state);
    f251_add(twice, state[0], state[0]);
    f251_add(state[0], sum, twice);
    // Use basic add/sub: the legacy assembly x-minus-c*y helper fails at y=0.
    f251_add(twice, state[1], state[1]);
    f251_sub(state[1], sum, twice);
    f251_add(twice, state[2], state[2]);
    f251_add(thrice, twice, state[2]);
    f251_sub(state[2], sum, thrice);
    */
    if (round_mode != FULL_ROUND_STANDARD && round_mode != PARTIAL_ROUND_STANDARD) return -1;
    felt_t input[3], product;
    const felt_t (*matrix)[3] = round_mode == FULL_ROUND_STANDARD
                                ? ctx->mds_full : ctx->mds_partial;
    /* Initialization guarantees matrix entries < p. Normalize each input
     * once, keeping both multiplication operands canonical to avoid the
     * field backend's near-2^256 overflow edge case. For x,y < p, REDC
     * returns < p + p^2/2^256 < 2p. Three products sum to < 6p < 2^256,
     * so only the completed row needs canonical reduction.
     */
    for (int j = 0; j < 3; j++) reduce_standard(input[j], state[j]);
    for (int i = 0; i < 3; i++) {
        f251_montgomery_mult(state[i], matrix[i][0], input[0]);
        for (int j = 1; j < 3; j++) {
            f251_montgomery_mult(product, matrix[i][j], input[j]);
            f251_add(state[i], state[i], product);
        }
        reduce_standard(state[i], state[i]);
    }
    return 0;
}

int round_3_pow5_standard(felt_t state[3], int round_idx, uint8_t round_mode,
                          const poseidon3_ctx_standard *ctx)
{
    if (round_idx < 0 || round_idx >= POSEIDON3_ROUNDS_STANDARD ||
        (round_mode != FULL_ROUND_STANDARD && round_mode != PARTIAL_ROUND_STANDARD))
        return -1;

    // AddRoundConstants: original constants on every element in every round.
    for (int i = 0; i < 3; i++) {
        f251_add(state[i], state[i], ctx->round_constants[round_idx][i]);
        reduce_standard(state[i], state[i]);
    }

    // SubWords: x^5 on all elements, or on the first element in partial rounds.
    if (round_mode == FULL_ROUND_STANDARD) {
        for (int i = 0; i < 3; i++)
            f251_montgomery_pow5(state[i], state[i]);
    } else {
        f251_montgomery_pow5(state[0], state[0]);
    }

    // MixLayer: the unmodified C base matrix, with no SD transformation.
    return mix_layer_3_standard(state, round_mode, ctx);
}

void permutation_3_pow5_standard(felt_t state[3], const poseidon3_ctx_standard *ctx)
{
    for (int i = 0; i < 3; i++) f251_to_montgomery(state[i], state[i]);
    for (int r = 0; r < POSEIDON3_ROUNDS_STANDARD; r++)
        round_3_pow5_standard(state, r,
            r < 4 || r >= 87 ? FULL_ROUND_STANDARD : PARTIAL_ROUND_STANDARD, ctx);
    for (int i = 0; i < 3; i++) f251_from_montgomery(state[i], state[i]);
}
