#include <inttypes.h>
#include <stdio.h>
#include <string.h>
#include "poseidon_profiles.h"

int main(int argc,char **argv){
    if(argc!=2){fprintf(stderr,"Usage: %s PROFILE\nProfiles:\n",argv[0]);for(int i=0;i<7;i++)fprintf(stderr,"  %s\n",poseidon_profiles[i].name);return 1;}
    const poseidon_profile *p=NULL;
    for(int i=0;i<7;i++)if(!strcmp(argv[1],poseidon_profiles[i].name))p=&poseidon_profiles[i];
    if(!p){fprintf(stderr,"Unknown profile: %s\n",argv[1]);return 1;}
    poseidon_profile_ctx plain,protected_ctx;
    felt_t a[POSEIDON_MAX_WIDTH]={{0}}, baseline[POSEIDON_MAX_WIDTH]={{0}},sd[POSEIDON_MAX_WIDTH]={{0}};
    for(int i=0;i<p->t;i++){a[i][0]=(uint64_t)(1+i%3);baseline[i][0]=sd[i][0]=(uint64_t)(i+1);}
    if(poseidon_profile_init(&plain,p,0,NULL)||poseidon_profile_init(&protected_ctx,p,1,a)){fprintf(stderr,"Initialization failed\n");return 1;}
    poseidon_profile_permute(&plain,baseline);poseidon_profile_permute(&protected_ctx,sd);
    printf("%s: t=%d RF=%d RP=%d, x^5, SD_PERMUTE=%d SD_SCALE=%d\n",p->name,p->t,p->rf,p->rp,SD_PERMUTE,SD_SCALE);
    for(int i=0;i<p->t;i++){
        printf("[%d] plain=0x",i);for(int j=3;j>=0;j--)printf("%016" PRIx64,baseline[i][j]);
        printf(" protected=0x");for(int j=3;j>=0;j--)printf("%016" PRIx64,sd[i][j]);putchar('\n');
    }
    return 0;
}
