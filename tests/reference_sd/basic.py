import os

import basic
from params import *

# Poseidon Constants Parameters

P = 0x73EDA753299D7D483339D80809A1D80553BDA402FFFE5BFEFFFFFFFF00000001
M = 128
T_RANGE = {3, 5, 9, 12}
S_BOX_EXP = 5

ROUNDFULL = 8
ROUNDPARTIAL = {3: 55, 5: 56, 9: 57, 12: 57}


def pi_matrix(matrix):
    if PERMUTE:
        return [matrix[2], matrix[0], matrix[1]]
    return matrix
#
def matrix_piinv(matrix):
    if PERMUTE:
        new_matrix = [[0,0,0],[0,0,0],[0,0,0]]

        new_matrix[0][0] = matrix[0][2]
        new_matrix[1][0] = matrix[1][2]
        new_matrix[2][0] = matrix[2][2]

        new_matrix[0][1] = matrix[0][0]
        new_matrix[1][1] = matrix[1][0]
        new_matrix[2][1] = matrix[2][0]

        new_matrix[0][2] = matrix[0][1]
        new_matrix[1][2] = matrix[1][1]
        new_matrix[2][2] = matrix[2][1]

        return new_matrix
    return matrix





# prime field element arithmetic operation
class PrimeFieldOps:
    def add(op1, op2):
        return (op1 + op2) % P

    def mul(op1, op2):
        return (op1 * op2) % P

    def s_box(op1):
        tmp = 1
        for i in range(S_BOX_EXP):
            tmp = PrimeFieldOps.mul(tmp, op1)
        return tmp

    def s_boxes(op1_list):
        tmp = []
        for i in op1_list:
            tmp.append(PrimeFieldOps.s_box(i))
        return tmp

    def add_round_constants(state, constant):
        tmp = []
        for i in range(len(state)):
            tmp.append(PrimeFieldOps.add(state[i], constant[i]))
        return tmp

    def get_inverse(x):
        """get the modular inverse of x (mod p )"""
        """ the inverse of x equals to pow(x,p-2) mod p """
        p_bin = bin(P - 2)[2:]
        inverse = 1
        tmp = x % P
        for i in range(len(p_bin)):
            if p_bin[len(p_bin) - i - 1] == "1":
                inverse = (inverse * tmp) % P
            tmp = (tmp * tmp) % P

        return inverse


    def mds_mixing(state, mode="full"):
        t = len(state)
        new_state = []

        mds_matrix = []
        if mode == "full":
            mds_matrix = PrimeFieldOps.get_mds_matrix(t)
        else:
            mds_matrix = PrimeFieldOps.get_mds_matrix_partial(t)

        if PERMUTE:
            mds_matrix = pi_matrix(matrix_piinv(mds_matrix))

        for i in range(t):
            tmp = 0
            for j in range(t):
                tmp = PrimeFieldOps.add(
                    PrimeFieldOps.mul(state[j], mds_matrix[i][j]), tmp
                )
            new_state.append(tmp)

        return new_state

    def init_state(preimage):
        """get initial state from preimage"""
        domain_tag = [pow(2, len(preimage)) - 1]
        return domain_tag + preimage

    def get_mds_matrix(t):

        mds_x = range(0, t)
        mds_y = range(t, 2 * t)
        mds_matrix = []

        for x in mds_x:
            mds_vec = []
            for y in mds_y:
                mds_vec.append(PrimeFieldOps.get_inverse(x + y))

            mds_matrix.append(mds_vec)


        # Change the matrix here
        mds_matrix[0][0] = basic.PrimeFieldOps.mul(mds_matrix[0][0],
                                                   basic.PrimeFieldOps.mul(
                                                       1,
                                                       PrimeFieldOps.get_inverse(A0_CONST ** (S_BOX_EXP - 1))
                                                   )
                                                   )

        mds_matrix[0][1] = basic.PrimeFieldOps.mul(mds_matrix[0][1],
                                                   basic.PrimeFieldOps.mul(
                                                       A0_CONST,
                                                       PrimeFieldOps.get_inverse(A1_CONST ** (S_BOX_EXP))
                                                   )
                                                   )

        mds_matrix[0][2] = basic.PrimeFieldOps.mul(mds_matrix[0][2],
                                                   basic.PrimeFieldOps.mul(
                                                       A0_CONST,
                                                       PrimeFieldOps.get_inverse(A2_CONST ** (S_BOX_EXP))
                                                   )
                                                   )

        mds_matrix[1][0] = basic.PrimeFieldOps.mul(mds_matrix[1][0],
                                                   basic.PrimeFieldOps.mul(
                                                       A1_CONST,
                                                       PrimeFieldOps.get_inverse(A0_CONST ** (S_BOX_EXP))
                                                   )
                                                   )

        mds_matrix[1][1] = basic.PrimeFieldOps.mul(mds_matrix[1][1],
                                                   basic.PrimeFieldOps.mul(
                                                       1,
                                                       PrimeFieldOps.get_inverse(A1_CONST ** (S_BOX_EXP - 1))
                                                   )
                                                   )

        mds_matrix[1][2] = basic.PrimeFieldOps.mul(mds_matrix[1][2],
                                                   basic.PrimeFieldOps.mul(
                                                       A1_CONST,
                                                       PrimeFieldOps.get_inverse(A2_CONST ** (S_BOX_EXP))
                                                   )
                                                   )

        mds_matrix[2][0] = basic.PrimeFieldOps.mul(mds_matrix[2][0],
                                                   basic.PrimeFieldOps.mul(
                                                       A2_CONST,
                                                       PrimeFieldOps.get_inverse(A0_CONST ** (S_BOX_EXP))
                                                   )
                                                   )

        mds_matrix[2][1] = basic.PrimeFieldOps.mul(mds_matrix[2][1],
                                                   basic.PrimeFieldOps.mul(
                                                       A2_CONST,
                                                       PrimeFieldOps.get_inverse(A1_CONST ** (S_BOX_EXP))
                                                   )
                                                   )

        mds_matrix[2][2] = basic.PrimeFieldOps.mul(mds_matrix[2][2],
                                                   basic.PrimeFieldOps.mul(
                                                       1,
                                                       PrimeFieldOps.get_inverse(A2_CONST ** (S_BOX_EXP - 1))
                                                   )
                                                   )

        return mds_matrix

    def get_mds_matrix_partial(t):

        mds_x = range(0, t)
        mds_y = range(t, 2 * t)
        mds_matrix = []

        for x in mds_x:
            mds_vec = []
            for y in mds_y:
                mds_vec.append(PrimeFieldOps.get_inverse(x + y))

            mds_matrix.append(mds_vec)

        #                                            )

        # if PERMUTE:
        # Change the matrix here
        mds_matrix[0][0] = basic.PrimeFieldOps.mul(mds_matrix[0][0],
                                                   basic.PrimeFieldOps.mul(
                                                       1,
                                                       PrimeFieldOps.get_inverse(A0_CONST ** (S_BOX_EXP - 1))
                                                   )
                                                   )

        mds_matrix[0][1] = basic.PrimeFieldOps.mul(mds_matrix[0][1],
                                                   basic.PrimeFieldOps.mul(
                                                       A0_CONST,
                                                       PrimeFieldOps.get_inverse(A1_CONST)
                                                   )
                                                   )

        mds_matrix[0][2] = basic.PrimeFieldOps.mul(mds_matrix[0][2],
                                                   basic.PrimeFieldOps.mul(
                                                       A0_CONST,
                                                       PrimeFieldOps.get_inverse(A2_CONST)
                                                   )
                                                   )

        mds_matrix[1][0] = basic.PrimeFieldOps.mul(mds_matrix[1][0],
                                                   basic.PrimeFieldOps.mul(
                                                       A1_CONST,
                                                       PrimeFieldOps.get_inverse(A0_CONST ** S_BOX_EXP)
                                                   )
                                                   )


        mds_matrix[1][2] = basic.PrimeFieldOps.mul(mds_matrix[1][2],
                                                   basic.PrimeFieldOps.mul(
                                                       A1_CONST,
                                                       PrimeFieldOps.get_inverse(A2_CONST)
                                                   )
                                                   )

        mds_matrix[2][0] = basic.PrimeFieldOps.mul(mds_matrix[2][0],
                                                   basic.PrimeFieldOps.mul(
                                                       A2_CONST,
                                                       PrimeFieldOps.get_inverse(A0_CONST ** S_BOX_EXP)
                                                   )
                                                   )

        mds_matrix[2][1] = basic.PrimeFieldOps.mul(mds_matrix[2][1],
                                                   basic.PrimeFieldOps.mul(
                                                       A2_CONST,
                                                       PrimeFieldOps.get_inverse(A1_CONST)
                                                   )
                                                   )


            # Change the matrix here ###

            # M0 = mds_matrix[0][0]
            # M1 = mds_matrix[0][1]
            # M2 = mds_matrix[0][2]
            # M3 = mds_matrix[1][0]
            # M4 = mds_matrix[1][1]
            # M5 = mds_matrix[1][2]
            # M6 = mds_matrix[2][0]
            # M7 = mds_matrix[2][1]
            # M8 = mds_matrix[2][2]
            #
            # mds_matrix[0][0] = M8
            #
            # mds_matrix[0][1] = basic.PrimeFieldOps.mul(M6,
            #                                            basic.PrimeFieldOps.mul(
            #                                                A2_CONST,
            #                                                PrimeFieldOps.get_inverse(A0_CONST ** S_BOX_EXP)
            #                                            )
            #                                            )
            #
            # mds_matrix[0][2] = basic.PrimeFieldOps.mul(M7,
            #                                            basic.PrimeFieldOps.mul(
            #                                                A2_CONST,
            #                                                PrimeFieldOps.get_inverse(A1_CONST)
            #                                            )
            #                                            )
            #
            # mds_matrix[1][0] = basic.PrimeFieldOps.mul(M2,
            #                                            basic.PrimeFieldOps.mul(
            #                                                A0_CONST,
            #                                                PrimeFieldOps.get_inverse(A2_CONST)
            #                                            )
            #                                            )
            #
            # mds_matrix[1][1] = basic.PrimeFieldOps.mul(M0,
            #                                            basic.PrimeFieldOps.mul(
            #                                                1,
            #                                                PrimeFieldOps.get_inverse(A0_CONST ** (S_BOX_EXP - 1))
            #                                            )
            #                                            )
            #
            # mds_matrix[1][2] = basic.PrimeFieldOps.mul(M1,
            #                                            basic.PrimeFieldOps.mul(
            #                                                A0_CONST,
            #                                                PrimeFieldOps.get_inverse(A1_CONST)
            #                                            )
            #                                            )
            #
            #
            # mds_matrix[2][0] = basic.PrimeFieldOps.mul(M5,
            #                                            basic.PrimeFieldOps.mul(
            #                                                A1_CONST,
            #                                                PrimeFieldOps.get_inverse(A2_CONST)
            #                                            )
            #                                            )
            #
            # mds_matrix[2][1] = basic.PrimeFieldOps.mul(M3,
            #                                            basic.PrimeFieldOps.mul(
            #                                                A1_CONST,
            #                                                PrimeFieldOps.get_inverse(A0_CONST ** (S_BOX_EXP))
            #                                            )
            #                                            )
            # mds_matrix[2][2] = M4



        return mds_matrix





    def output_mds_matrixs():
        """generate mds_matrixs and write to files"""

        os.mkdir("mds_matrixs")
        os.chdir("mds_matrixs")
        for t in T_RANGE:
            mds_matrix = PrimeFieldOps.get_mds_matrix(t)
            fileobject = open("mds_matrix_{}.txt".format(t), "w")
            for mds_vec in mds_matrix:
                for element in mds_vec:
                    output_str = hex(element)
                    fileobject.write(output_str[2:] + "\n")

            fileobject.close()
