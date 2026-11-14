#ifndef F251_PARAMS_H
#define F251_PARAMS_H

/* Fixed StarkWare field: p = 2^251 + 17 * 2^192 + 1.
 * Shared by C and preprocessed assembly. Changing this definition alone
 * cannot change the field: reduction and Montgomery constants are specialized.
 */
#define F251_MODULUS_HIGH 0x0800000000000011
#define F251_MODULUS_INITIALIZER {1, 0, 0, F251_MODULUS_HIGH}

#endif
