import ctypes as c, sys, os
sys.path.insert(0,'tools')
from poseidon_bls_parameters import generate_constants
p=2**251+17*2**192+1
F=c.c_uint64*4
S=F*10
Matrix=S*10
RC=S*128
class Profile(c.Structure): _fields_=[('name',c.c_char_p),('t',c.c_int),('rf',c.c_int),('rp',c.c_int)]
class Context(c.Structure): _fields_=[('t',c.c_int),('rf',c.c_int),('rp',c.c_int),('rc',RC),('full',Matrix),('partial',Matrix),('multipliers',S),('protected_version',c.c_int)]
lib=c.CDLL(os.environ.get('POSEIDON_PROFILE_LIBRARY', '/tmp/poseidon_profiles.so')); lib.poseidon_profile_init.argtypes=[c.POINTER(Context),c.POINTER(Profile),c.c_int,c.POINTER(F)];lib.poseidon_profile_permute.argtypes=[c.POINTER(Context),c.POINTER(F)]
profiles=(Profile*7).in_dll(lib,'poseidon_profiles')
def integer(a):return sum(int(a[i])<<(64*i) for i in range(4))
def felt(x):return F(*[(x>>(64*i))&((1<<64)-1) for i in range(4)])
rotate=int(os.environ.get('TEST_PERMUTE','1'))
scale=int(os.environ.get('TEST_SCALE','1'))
for profile in profiles:
 t,rf,rp=profile.t,profile.rf,profile.rp
 const=generate_constants(t,rf,rp,p)
 x=list(range(1,t+1)); a=[1+i%3 for i in range(t)]
 for r in range(rf+rp):
  x=[(x[i]+const[r][i])%p for i in range(t)]
  for i in (range(t) if r<4 or r>=4+rp else (0,)):x[i]=pow(x[i],5,p)
  x=[sum(pow(t+i+j,-1,p)*x[j] for j in range(t))%p for i in range(t)]
 for protected in (0,1):
  ctx=Context(); mult=S(*[felt(v) for v in a]);assert lib.poseidon_profile_init(c.byref(ctx),c.byref(profile),protected,mult)==0
  state=S(*[felt(i+1) for i in range(t)]);lib.poseidon_profile_permute(c.byref(ctx),state)
  expected=([x[i]*a[i]%p for i in range(t)] if protected and scale else x[:])
  if protected and rotate: expected=expected[-1:]+expected[:-1]
  assert [integer(v) for v in state[:t]]==expected,(profile.name,protected)
 print(profile.name.decode(), 'passed')
