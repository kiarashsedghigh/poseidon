#define _POSIX_C_SOURCE 200809L
#include <errno.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "config.h"
#include "poseidon_profiles.h"

static volatile uint64_t sink;
static double elapsed(const struct timespec *a,const struct timespec *b) {
    return (double)(b->tv_sec-a->tv_sec)*1e9+(double)(b->tv_nsec-a->tv_nsec);
}
static int run(const poseidon_profile *p,uint64_t iterations,uint64_t warmup) {
    poseidon_profile_ctx plain,sd;
    felt_t multipliers[POSEIDON_MAX_WIDTH]={{0}}, input[POSEIDON_MAX_WIDTH]={{0}};
    felt_t prep[POSEIDON_MAX_WIDTH]={{0}}, state[POSEIDON_MAX_WIDTH];
    double total[2]={0,0};
    for(int i=0;i<p->t;i++){input[i][0]=(uint64_t)i+1; multipliers[i][0]=(uint64_t)(1+i%3);}
    if(poseidon_profile_init(&plain,p,0,NULL)||poseidon_profile_init(&sd,p,1,multipliers))return -1;
    for(int i=0;i<p->t;i++) f251_to_montgomery(prep[i],input[i]);
    const poseidon_profile_ctx *ctx[2]={&plain,&sd};
    for(int mode=0;mode<2;mode++) {
        for(uint64_t i=0;i<warmup;i++) {
            memcpy(state,prep,sizeof(state));if(mode)poseidon_profile_encode(&sd,state);poseidon_profile_rounds(ctx[mode],state);sink+=state[0][0];
        }
        for(uint64_t i=0;i<iterations;i++) {
            struct timespec a,b;
            memcpy(state,prep,sizeof(state));
            if(clock_gettime(CLOCK_MONOTONIC,&a))return -1;
            if(mode)poseidon_profile_encode(&sd,state);
            poseidon_profile_rounds(ctx[mode],state);
            if(clock_gettime(CLOCK_MONOTONIC,&b))return -1;
            total[mode]+=elapsed(&a,&b);sink+=state[0][0];
        }
    }
    printf("%-17s t=%2d RF=%d RP=%3d plain=%9.3f us protected=%9.3f us overhead=%+8.2f%%\n",
           p->name,p->t,p->rf,p->rp,total[0]/iterations/1e3,total[1]/iterations/1e3,
           total[0]>0?100.0*(total[1]/total[0]-1):0.0);
    return 0;
}
static int parse_count(const char *s,uint64_t *out,int allow_zero) {
    char *end;errno=0;if(!*s||*s=='-')return -1;
    unsigned long long n=strtoull(s,&end,10);
    if(errno||*end||(!allow_zero&&!n))return -1;
    *out=(uint64_t)n;return 0;
}
int main(int argc,char **argv) {
    const char *name="all";uint64_t iterations=ITER,warmup=WARMUP;
    if(argc>4 || (argc>2&&parse_count(argv[2],&iterations,0)) ||
       (argc>3&&parse_count(argv[3],&warmup,1))) {
        fprintf(stderr,"Usage: %s [all|PROFILE [ITERATIONS [WARMUP]]]\n",argv[0]);return 1;
    }
    if(argc>1)name=argv[1];
    int found=0;
    for(int i=0;i<7;i++)if(!strcmp(name,"all")||!strcmp(name,poseidon_profiles[i].name))found=1;
    if(!found){fprintf(stderr,"Unknown profile: %s\n",name);return 1;}
    printf("Plain: rounds; protected: input scaling/rotation + rounds; iterations=%" PRIu64 ", warmup=%" PRIu64 ", rotation=%d, scaling=%d\n",iterations,warmup,SD_PERMUTE,SD_SCALE);
    for(int i=0;i<7;i++)if(!strcmp(name,"all")||!strcmp(name,poseidon_profiles[i].name))
        if(run(&poseidon_profiles[i],iterations,warmup)){perror("benchmark initialization or clock_gettime");return 1;}
    printf("Checksum: %016" PRIx64 "\n",sink);return 0;
}
