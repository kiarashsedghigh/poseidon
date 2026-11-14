#ifndef POSEIDON_PROFILES_H
#define POSEIDON_PROFILES_H
#include "f251.h"
#include "poseidon3_sd_config.h"
#define POSEIDON_MAX_WIDTH 10
#define POSEIDON_MAX_ROUNDS 128

typedef struct { const char *name; int t, rf, rp; } poseidon_profile;
extern const poseidon_profile poseidon_profiles[7];
/* Contexts are initialized once and must remain read-only while used. */
typedef struct {
    int t, rf, rp;
    felt_t rc[POSEIDON_MAX_ROUNDS][POSEIDON_MAX_WIDTH];
    felt_t full[POSEIDON_MAX_WIDTH][POSEIDON_MAX_WIDTH];
    felt_t partial[POSEIDON_MAX_WIDTH][POSEIDON_MAX_WIDTH];
    felt_t multipliers[POSEIDON_MAX_WIDTH];
    int protected_version;
} poseidon_profile_ctx;
/* multipliers are ordinary values; protected_version selects the SD path. */
int poseidon_profile_init(poseidon_profile_ctx *, const poseidon_profile *, int protected_version, const felt_t *multipliers);
/* Execute all rounds on a prepared Montgomery state, without conversions. */
void poseidon_profile_rounds(const poseidon_profile_ctx *, felt_t *state);
/* Input/output are ordinary field values; protected output retains encoding. */
void poseidon_profile_permute(const poseidon_profile_ctx *, felt_t *state);
#endif
