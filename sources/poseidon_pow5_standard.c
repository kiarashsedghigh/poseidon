#include "poseidon_pow5_standard.h"
#include "poseidon_rc.h"

void mix_layer_3_standard(felt_t state[3])
{
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
}

int round_3_pow5_standard(felt_t state[3], int round_idx, uint8_t round_mode)
{
    if (round_idx < 0 || round_idx >= POSEIDON3_ROUNDS_STANDARD ||
        (round_mode != FULL_ROUND_STANDARD && round_mode != PARTIAL_ROUND_STANDARD))
        return -1;

    // AddRoundConstants: original constants on every element in every round.
    for (int i = 0; i < 3; i++)
        f251_add(state[i], state[i], ROUND_CONST_MONT_P3[round_idx][i]);

    // SubWords: x^5 on all elements, or on the first element in partial rounds.
    if (round_mode == FULL_ROUND_STANDARD) {
        for (int i = 0; i < 3; i++)
            f251_montgomery_pow5(state[i], state[i]);
    } else {
        f251_montgomery_pow5(state[0], state[0]);
    }

    // MixLayer: the unmodified C base matrix, with no SD transformation.
    mix_layer_3_standard(state);
    return 0;
}

void permutation_3_pow5_standard(felt_t state[3])
{
    for (int i = 0; i < 3; i++) f251_to_montgomery(state[i], state[i]);
    for (int r = 0; r < POSEIDON3_ROUNDS_STANDARD; r++)
        round_3_pow5_standard(state, r,
            r < 4 || r >= 87 ? FULL_ROUND_STANDARD : PARTIAL_ROUND_STANDARD);
    for (int i = 0; i < 3; i++) f251_from_montgomery(state[i], state[i]);
}
