import basic
import constants
from params import *
from basic import *
# poseidon hash function

def poseidon_hash(preimage):
    t = len(preimage) + 1
    if t not in basic.T_RANGE:
        print("error: the length of preimage is incorrect")
        exit()

    roundf = basic.ROUNDFULL
    roundp = basic.ROUNDPARTIAL[t]


    ############ CONSTANT
    round_constants = constants.generate_constants(t, roundf, roundp)

    for i in range(0, len(round_constants),3):
        round_constants[i] = basic.PrimeFieldOps.mul(round_constants[i], A0_CONST)
        round_constants[i + 1] = basic.PrimeFieldOps.mul(round_constants[i + 1], A1_CONST)
        round_constants[i + 2] = basic.PrimeFieldOps.mul(round_constants[i + 2], A2_CONST)

        a, b, c = round_constants[i], round_constants[i+1], round_constants[i+2]
        if PERMUTE:
            round_constants[i] = c
            round_constants[i+1] = a
            round_constants[i + 2] = b

    ############ STATE
    state = basic.PrimeFieldOps.init_state(preimage)

    for i in range(0, len(state),3):
        state[i] = basic.PrimeFieldOps.mul(state[i], A0_CONST)
        state[i + 1] = basic.PrimeFieldOps.mul(state[i + 1], A1_CONST)
        state[i + 2] = basic.PrimeFieldOps.mul(state[i + 2], A2_CONST)

        a, b, c = state[i], state[i+1], state[i+2]
        if PERMUTE:
            state[i] = c
            state[i+1] = a
            state[i + 2] = b

    ############ ROUNDS
    for i in range(int(roundf / 2)):
        state = basic.PrimeFieldOps.add_round_constants(
            state, round_constants[i * t : (i + 1) * t]
        )
        state = basic.PrimeFieldOps.s_boxes(state)
        state = basic.PrimeFieldOps.mds_mixing(state)


    for i in range(int(roundf / 2), int(roundf / 2 + roundp)):

        state = basic.PrimeFieldOps.add_round_constants(
            state, round_constants[i * t : (i + 1) * t]
        )
        # # state = [s2 s0 s1]
        if PERMUTE:
            state[1] = basic.PrimeFieldOps.s_box(state[1])
        else:
            state[0] = basic.PrimeFieldOps.s_box(state[0])

        state = basic.PrimeFieldOps.mds_mixing(state, "partial")


    for i in range(int(roundf / 2 + roundp), int(roundf + roundp)):
        state = basic.PrimeFieldOps.add_round_constants(
            state, round_constants[i * t : (i + 1) * t]
        )
        state = basic.PrimeFieldOps.s_boxes(state)
        state = basic.PrimeFieldOps.mds_mixing(state)

    print("=================")
    print(f"A0: {state[0]}")
    print(f"A1: {state[1]}")
    print(f"A2: {state[2]}")

    return state[1]


x = poseidon_hash([1,2])
