"""Small exact checks using a separate coefficient implementation."""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import permutations
from math import comb, factorial
from pathlib import Path
import json

def add(*polys):
    out = [F(0)] * max(map(len, polys))
    for p in polys:
        for i, c in enumerate(p):
            out[i] += c
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return tuple(out)

def scale(p, c):
    return tuple(c * a for a in p)

def mul(p, q):
    out = [F(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i+j] += a*b
    return tuple(out)

def value(p, x):
    return sum((c*x**i for i, c in enumerate(p)), F(0))

def choose(x, j):
    out = F(1)
    for i in range(j):
        out *= F(x-i, i+1)
    return out

def valid(s):
    return not s or (s[0] >= 2 and all(b-a >= 2 for a,b in zip(s,s[1:])))

@lru_cache(None)
def polynomial(s):
    if not s:
        return (F(1),)
    if not valid(s):
        return (F(0),)
    m, parent = s[-1], s[:-1]
    basis = (F(1),)
    for i in range(m-1):
        basis = scale(mul(basis, (-i, 1)), F(1,i+1))
    return add(scale(basis,value(polynomial(parent),m-1)),
               scale(polynomial(parent),-2), scale(polynomial(parent+(m-1,)),-1))

def transform(p, m):
    return tuple((-1)**i * comb(m,i) * value(p,m-i) for i in range(m+1))

def moment(d, y):
    return sum((c/(y+i) for i,c in enumerate(d)),F(0))

def convolution(d, t, r):
    return sum((c*comb(r,b)*(-t)**(r-b)*(1-t**(i+b+1))/F(i+b+1)
                for i,c in enumerate(d) for b in range(r+1)),F(0))

def sign(x):
    return (x>0)-(x<0)

def kernel(k, y, roots):
    t,u = (r*r for r in roots)
    out = F((-1)**k,2)*sum(((-1)**j*choose(y-1,j)*(t**j+u**j)
                            for j in range(k-1)), F(0))
    if k%2 == 0:
        out -= sum((r**(int(2*y)-2) for r in roots),F(0))
    return out

def all_sets(bound):
    for mask in range(1, 1<<(bound-1)):
        s=tuple(i+2 for i in range(bound-1) if mask&(1<<i))
        if valid(s):
            yield s

stats = Counter()
rows=[]
for s in all_sets(10):
    p=polynomial(s)
    m=s[-1]; parent=s[:-1]; a=parent[-1] if parent else 0; k=m-a
    f=polynomial(parent); d=transform(f,a); ds=transform(p,m)
    assert value(p,m)==0
    derivative=sum((i*c*m**(i-1) for i,c in enumerate(p) if i),F(0))
    assert derivative>0, ('derivative',s)
    assert value(ds,0)==value(ds,1)==0
    c=F(factorial(m),factorial(a)*factorial(k-1))
    # m+1 distinct values determine both degree-at-most-m polynomials exactly.
    for j in range(m+1):
        t=F(j,m)
        actual=value(ds,t)
        assert actual==value(ds,1-t), ('symmetry',s,t)
        predicted=c*(moment(d,F(k))-convolution(d,t,k-1)-convolution(d,1-t,k-1))
        assert actual==predicted, ('convolution identity',s,t)
        assert j in (0,m) or actual>0, ('D positivity',s,t)
        stats['transform_identity_points']+=1
    for y in sorted({F(1,1000),F(1,2),F(999,1000),F(1),F(1001,1000),
                     F(k)-F(1,1000),F(k),F(k)+F(1,1000),F(k)+F(1,2),F(3*k)}):
        x=a+y
        pref=F(1,factorial(a))
        for j in range(a+1): pref*=x-j
        assert value(f,x)==pref*moment(d,y), ('moment',s,x)
        integral=sum(((-1)**(k-2-j)*choose(y-1,j)*moment(d,F(j+1))
                      for j in range(k-1)),F(0))
        if k%2==0: integral-=2*moment(d,y)
        assert value(p,x)==pref*integral, ('main integral',s,x)
        assert sign(value(p,x))==sign(y-k), ('main sign',s,x)
        stats['rational_sign_moment_points']+=1
    stats['nonempty_sets']+=1
    stats['original_target_sets']+=len(s)>=2
    rows.append({'S':s,'monomial_ascending':[[v.numerator,v.denominator] for v in p],
                 'derivative_at_m':[derivative.numerator,derivative.denominator]})

for n in range(1,9):
    counts=Counter(tuple(i+1 for i in range(1,n-1) if p[i-1]<p[i]>p[i+1])
                   for p in permutations(range(n)))
    for s,count in counts.items():
        assert value(polynomial(s),n)*2**(n-len(s)-1)==count, ('permutation',n,s)
    stats['permutations']+=factorial(n)

for k in range(2,11):
    for p in range(1,4*k+1):
        y=F(p,2)
        for roots in ((F(3,5),F(4,5)),(F(5,13),F(12,13))):
            t=roots[0]**2
            v=kernel(k,y,roots)
            assert sign(v)==sign(y-k), ('kernel',k,y,t)
            rhs=choose(y-1,k-1)/F(2*k)*((y-k)*(t**k+(1-t)**k)-k*(t**(k-1)+(1-t)**(k-1)))
            assert kernel(k+2,y,roots)-v==rhs, ('two-step',k,y,t)
            stats['kernel_points']+=1
            stats['kernel_zero_points']+=y==k

out={'status':'PASS','scope':'exact finite falsification only; not a universal proof',
     'statistics':dict(stats),'failures':[], 'polynomials':rows}
Path('independent_results.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
print(json.dumps({k:v for k,v in out.items() if k!='polynomials'},indent=2))
