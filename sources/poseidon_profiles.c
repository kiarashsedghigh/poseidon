#include "poseidon_profiles.h"
#include <string.h>

const poseidon_profile poseidon_profiles[7] = {
    {"poseidon-128-t3",3,8,57}, {"poseidon-128-t5",5,8,60},
    {"poseidon-80-t3",3,8,33}, {"poseidon-80-t5",5,8,35},
    {"poseidon-256-t6",6,8,120}, {"poseidon-256-t10",10,8,120},
    {"poseidon-x",9,8,57}
};
/* Canonicalize arbitrary 256-bit field representatives. */
static void reduce(felt_t z, const felt_t x) {
    uint64_t q=x[3]>>59, sub[4]={q,0,0,q*F251_MODULUS[3]}, tmp[4];
    uint64_t borrow=0,carry=0;
    for(int i=0;i<4;i++){uint64_t a=x[i]-sub[i], b=x[i]<sub[i];tmp[i]=a-borrow;borrow=b|(a<borrow);}
    uint64_t mask=0-borrow;
    for(int i=0;i<4;i++){uint64_t a=tmp[i]+(F251_MODULUS[i]&mask),c=a<tmp[i];z[i]=a+carry;carry=c|(z[i]<a);}
}
static void mul(felt_t z,const felt_t x,const felt_t y){felt_t a,b;reduce(a,x);reduce(b,y);f251_montgomery_mult(z,a,b);reduce(z,z);}
static void add(felt_t z,const felt_t x,const felt_t y){f251_add(z,x,y);reduce(z,z);}
static void inv(felt_t out,const felt_t x){
    const felt_t exp={UINT64_MAX,UINT64_MAX,UINT64_MAX,F251_MODULUS_HIGH-1ULL};
    const felt_t one={1,0,0,0};felt_t v;f251_to_montgomery(v,one);
    for(int bit=251;bit>=0;bit--){mul(v,v,v);if((exp[bit/64]>>(bit%64))&1)mul(v,v,x);}
    memcpy(out,v,sizeof(v));
}
/* Grain LFSR, 252-bit prime field, x^5 S-box; rejection sampling. */
static unsigned step(unsigned char s[80]) {unsigned b=s[0]^s[13]^s[23]^s[38]^s[51]^s[62];memmove(s,s+1,79);s[79]=(unsigned char)b;return b;}
static void bits(unsigned char s[80],int *offset,unsigned v,int n){for(int k=n-1;k>=0;k--)s[(*offset)++]=(v>>k)&1;}
static int less_p(const felt_t x){for(int i=3;i>=0;i--)if(x[i]!=F251_MODULUS[i])return x[i]<F251_MODULUS[i];return 0;}
static void constants(poseidon_profile_ctx *c){
    unsigned char s[80];int o=0;bits(s,&o,1,2);bits(s,&o,1,4);bits(s,&o,252,12);bits(s,&o,c->t,12);bits(s,&o,c->rf,10);bits(s,&o,c->rp,10);while(o<80)s[o++]=1;
    for(int i=0;i<160;i++)step(s);
    for(int r=0;r<c->rf+c->rp;r++)for(int j=0;j<c->t;j++){
        felt_t x;
        do {memset(x,0,sizeof(x));int n=0;while(n<252){unsigned a=step(s),b=step(s);if(a){for(int k=3;k>0;k--)x[k]=(x[k]<<1)|(x[k-1]>>63);x[0]=(x[0]<<1)|b;n++;}}}while(!less_p(x));
        f251_to_montgomery(c->rc[r][j],x);reduce(c->rc[r][j],c->rc[r][j]);
    }
}
int poseidon_profile_init(poseidon_profile_ctx *c,const poseidon_profile *p,int protected_version,const felt_t *multipliers){
    if(!c||!p||p->t<1||p->t>10||p->rf!=8||p->rp<0||p->rf+p->rp>128|| (protected_version!=0&&protected_version!=1))return -1;
#if SD_SCALE
    if(protected_version&&!multipliers)return -1;
#endif
    poseidon_profile_ctx tmp;memset(&tmp,0,sizeof(tmp));tmp.t=p->t;tmp.rf=p->rf;tmp.rp=p->rp;tmp.protected_version=protected_version;
    felt_t inverse[10],inverse5[10];const felt_t one={1,0,0,0};
    for(int i=0;i<p->t;i++){
        felt_t value;memcpy(value,one,sizeof(value));
#if SD_SCALE
        if(protected_version)reduce(value,multipliers[i]);
#else
        (void)multipliers;
#endif
        if(!(value[0]|value[1]|value[2]|value[3]))return -1;
        f251_to_montgomery(tmp.multipliers[i],value);reduce(tmp.multipliers[i],tmp.multipliers[i]);
        inv(inverse[i],tmp.multipliers[i]);f251_montgomery_pow5(inverse5[i],inverse[i]);reduce(inverse5[i],inverse5[i]);
    }
    constants(&tmp);
    felt_t original[128][10];memcpy(original,tmp.rc,sizeof(original));
    for(int r=0;r<p->rf+p->rp;r++)for(int i=0;i<p->t;i++){
        int row=i;
#if SD_PERMUTE
        if(protected_version)row=(i+p->t-1)%p->t;
#endif
        mul(tmp.rc[r][i],original[r][row],tmp.multipliers[row]);
    }
    for(int i=0;i<p->t;i++)for(int j=0;j<p->t;j++){
        int row=i,col=j;
#if SD_PERMUTE
        if(protected_version){row=(i+p->t-1)%p->t;col=(j+p->t-1)%p->t;}
#endif
        felt_t x={(uint64_t)(p->t+row+col),0,0,0},base;
        f251_to_montgomery(x,x);inv(base,x);mul(base,base,tmp.multipliers[row]);
        mul(tmp.full[i][j],base,inverse5[col]);mul(tmp.partial[i][j],base,col==0?inverse5[col]:inverse[col]);
    }
    memcpy(c,&tmp,sizeof(tmp));return 0;
}
void poseidon_profile_permute(const poseidon_profile_ctx *c,felt_t *state){
    felt_t tmp[10],scaled[10];int t=c->t;
    for(int i=0;i<t;i++){reduce(tmp[i],state[i]);f251_to_montgomery(tmp[i],tmp[i]);mul(scaled[i],tmp[i],c->multipliers[i]);}
    for(int i=0;i<t;i++){int j=i;
#if SD_PERMUTE
        if(c->protected_version)j=(i+t-1)%t;
#endif
        memcpy(state[i],scaled[j],sizeof(felt_t));}
    for(int r=0;r<c->rf+c->rp;r++){
        int full=r<c->rf/2||r>=c->rf/2+c->rp;
        for(int i=0;i<t;i++)add(state[i],state[i],c->rc[r][i]);
        for(int i=0;i<t;i++){int active=0;
#if SD_PERMUTE
            if(c->protected_version)active=1%t;
#endif
            if(full||i==active)f251_montgomery_pow5(state[i],state[i]);}
        const felt_t (*m)[10]=full?c->full:c->partial;
        for(int i=0;i<t;i++){felt_t product;memset(tmp[i],0,sizeof(felt_t));for(int j=0;j<t;j++){mul(product,m[i][j],state[j]);add(tmp[i],tmp[i],product);}}
        memcpy(state,tmp,(size_t)t*sizeof(felt_t));
    }
    for(int i=0;i<t;i++)f251_from_montgomery(state[i],state[i]);
}
