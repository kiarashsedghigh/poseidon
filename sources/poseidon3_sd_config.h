#ifndef POSEIDON3_SD_CONFIG_H
#define POSEIDON3_SD_CONFIG_H

/* Compile every SD caller and implementation with the same options. */
#ifndef SD_PERMUTE
#define SD_PERMUTE 1
#endif
#ifndef SD_SCALE
#define SD_SCALE 1
#endif
#if (SD_PERMUTE != 0 && SD_PERMUTE != 1) || (SD_SCALE != 0 && SD_SCALE != 1)
#error "SD_PERMUTE and SD_SCALE must be 0 or 1"
#endif

#endif
