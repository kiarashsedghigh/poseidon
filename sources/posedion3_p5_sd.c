#include <string.h>
#include "posedion3_p5_sd.h"
#include "poseidon_rc.h"

const int poseidon3_permute_sd = SD_PERMUTE;
const int poseidon3_scale_sd = SD_SCALE;
#if SD_PERMUTE
static const int order_sd[3] = {2, 0, 1};
#define PARTIAL_INDEX_SD 1
#else
static const int order_sd[3] = {0, 1, 2};
#define PARTIAL_INDEX_SD 0
#endif

/* Reduce any 256-bit representative modulo the C field prime.
 * q=floor(x/2^251); subtract q*p, then conditionally add p after a borrow.
 * Since p=2^251+17*2^192+1, the provisional result is in (-p,p).
 */
static void reduce_sd(felt_t z, const felt_t x)
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

#if SD_SCALE
static void multiply_sd(felt_t z, const felt_t x, const felt_t y)
{
    felt_t a, b;
    reduce_sd(a, x);
    reduce_sd(b, y);
    f251_montgomery_mult(z, a, b);
    reduce_sd(z, z);
}

static int inverse_sd(felt_t z, const felt_t x)
{
    const felt_t exponent = {
        0xffffffffffffffffULL, 0xffffffffffffffffULL,
        0xffffffffffffffffULL, (F251_MODULUS_HIGH - 1ULL)
    };
    const felt_t one = {1, 0, 0, 0};
    felt_t result, product, a;
    reduce_sd(a, x);
    uint64_t nonzero = a[0] | a[1] | a[2] | a[3];
    f251_to_montgomery(result, one);
    for (int bit = 251; bit >= 0; bit--) {
        multiply_sd(result, result, result);
        multiply_sd(product, result, a);
        uint64_t mask = (uint64_t)0 - ((exponent[bit/64] >> (bit%64)) & 1);
        for (int i = 0; i < 4; i++)
            result[i] = (product[i] & mask) | (result[i] & ~mask);
    }
    memcpy(z, result, sizeof(result));
    return nonzero ? 0 : -1;
}

#endif

int poseidon3_init_sd(poseidon3_ctx_sd *ctx, const felt_t multipliers[3])
{
    poseidon3_ctx_sd tmp = {0};
#if SD_SCALE
    felt_t inverse[3], inverse_power[3];
#endif
    const int matrix[3][3] = {{3, 1, 1}, {1, -1, 1}, {1, 1, -2}};
    memcpy(tmp.modulus, F251_MODULUS, sizeof(tmp.modulus));
#if SD_SCALE
    for (int i = 0; i < 3; i++) {
        felt_t a;
        reduce_sd(a, multipliers[i]);
        f251_to_montgomery(tmp.multipliers[i], a);
        reduce_sd(tmp.multipliers[i], tmp.multipliers[i]);
        if (inverse_sd(inverse[i], tmp.multipliers[i])) return -1;
        f251_montgomery_pow5(inverse_power[i], inverse[i]);
    }
#else
    (void)multipliers;
    const felt_t one = {1, 0, 0, 0};
    for (int i = 0; i < 3; i++) {
        f251_to_montgomery(tmp.multipliers[i], one);
        reduce_sd(tmp.multipliers[i], tmp.multipliers[i]);
    }
#endif
    for (int r = 0; r < POSEIDON3_ROUNDS_SD; r++) {
        for (int i = 0; i < 3; i++) {
            int source = order_sd[i];
#if SD_SCALE
            multiply_sd(tmp.round_constants[r][i], ROUND_CONST_MONT_P3[r][source],
                         tmp.multipliers[source]);
#else
            memcpy(tmp.round_constants[r][i], ROUND_CONST_MONT_P3[r][source], sizeof(felt_t));
#endif
        }
    }
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            int row = order_sd[i];
            int col = order_sd[j];
            int entry = matrix[row][col];
            felt_t base = {(uint64_t)(entry < 0 ? -entry : entry), 0, 0, 0};
            if (entry < 0) {
                f251_sub(base, tmp.modulus, base);
                reduce_sd(base, base);
            }
            f251_to_montgomery(base, base);
#if SD_SCALE
            multiply_sd(base, base, tmp.multipliers[row]);
            // Full: D M D^-5. Partial: exponent 5 only on the first column.
            // PERMUTE conjugates both matrices by the [2,0,1] rotation.
            multiply_sd(tmp.mds_full[i][j], base, inverse_power[col]);
            multiply_sd(tmp.mds_partial[i][j], base,
                         col == 0 ? inverse_power[col] : inverse[col]);
#else
            reduce_sd(base, base);
            memcpy(tmp.mds_full[i][j], base, sizeof(felt_t));
            memcpy(tmp.mds_partial[i][j], base, sizeof(felt_t));
#endif
        }
    }
    memcpy(ctx, &tmp, sizeof(tmp));
    return 0;
}

int mix_layer_3_sd(felt_t state[3], uint8_t round_mode, const poseidon3_ctx_sd *ctx)
{
    if (round_mode != FULL_ROUND_SD && round_mode != PARTIAL_ROUND_SD) return -1;
#if SD_SCALE
    felt_t result[3] = {{0}}, product;
    const felt_t (*matrix)[3] = round_mode == FULL_ROUND_SD
                                ? ctx->mds_full : ctx->mds_partial;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            multiply_sd(product, matrix[i][j], state[j]);
            f251_add(result[i], result[i], product);
            reduce_sd(result[i], result[i]);
        }
    }
    memcpy(state, result, sizeof(result));
#else
    /* With unit scaling, both matrices are Pi*M*Pi^-1. Their rows are
     * sum+2*x0, sum-2*x1, sum-3*x2, at compile-time rotated positions.
     * Compute sum before overwriting any input. Use basic add/sub because
     * the legacy assembly x-minus-c*y helper has a zero-input edge case.
     */
    (void)ctx;
#if SD_PERMUTE
    enum { lane0 = 1, lane1 = 2, lane2 = 0 };
#else
    enum { lane0 = 0, lane1 = 1, lane2 = 2 };
#endif
    felt_t sum, twice, thrice;
    f251_sum_state_3(sum, state);
    f251_add(twice, state[lane0], state[lane0]);
    f251_add(state[lane0], sum, twice);
    f251_add(twice, state[lane1], state[lane1]);
    f251_sub(state[lane1], sum, twice);
    f251_add(twice, state[lane2], state[lane2]);
    f251_add(thrice, twice, state[lane2]);
    f251_sub(state[lane2], sum, thrice);
#endif
    return 0;
}

int round_3_sd(felt_t state[3], int round_idx, uint8_t round_mode,
               const poseidon3_ctx_sd *ctx)
{
    if (round_idx < 0 || round_idx >= POSEIDON3_ROUNDS_SD ||
        (round_mode != FULL_ROUND_SD && round_mode != PARTIAL_ROUND_SD)) return -1;

    // AddRoundConstants: scaled/permuted constants on every element.
    for (int i = 0; i < 3; i++) {
        f251_add(state[i], state[i], ctx->round_constants[round_idx][i]);
        reduce_sd(state[i], state[i]);
    }

    // SubWords: fifth power; the first element moves to index 1 if rotated.
    if (round_mode == FULL_ROUND_SD) {
        for (int i = 0; i < 3; i++) f251_montgomery_pow5(state[i], state[i]);
    } else {
        f251_montgomery_pow5(state[PARTIAL_INDEX_SD], state[PARTIAL_INDEX_SD]);
    }

    // MixLayer: full and partial rounds use different transformed matrices.
    return mix_layer_3_sd(state, round_mode, ctx);
}

void permutation_3_sd(felt_t state[3], const poseidon3_ctx_sd *ctx)
{
    // Convert, scale, and optionally rotate the input once before the rounds.
    felt_t scaled[3];
    for (int i = 0; i < 3; i++) {
        reduce_sd(scaled[i], state[i]);
        f251_to_montgomery(scaled[i], scaled[i]);
#if SD_SCALE
        multiply_sd(scaled[i], scaled[i], ctx->multipliers[i]);
#endif
    }
    for (int i = 0; i < 3; i++)
        memcpy(state[i], scaled[order_sd[i]], sizeof(felt_t));

    for (int r = 0; r < POSEIDON3_ROUNDS_SD; r++)
        round_3_sd(state, r, r < 4 || r >= 87 ? FULL_ROUND_SD : PARTIAL_ROUND_SD, ctx);

    // Return all three elements; scaling/rotation remain in the SD output.
    for (int i = 0; i < 3; i++) f251_from_montgomery(state[i], state[i]);
}
