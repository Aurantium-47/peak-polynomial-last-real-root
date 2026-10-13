#!/usr/bin/env python3
"""Exact independent checks accompanying proof.tex/proof.pdf.
Python >=3.10 and SymPy are required. No floating point or numerical quadrature.
  python verify_peak.py --max 18 --kernel-max 24 --output data
Optional --input-zip audits all input polynomials, NOT all input root counts.
The finite tests are consistency checks, not premises of the universal proof.
"""
from __future__ import annotations
import argparse, gzip, hashlib, json, platform, sys, zipfile
from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache, reduce
from itertools import permutations
from math import comb, factorial, gcd
from pathlib import Path


def need(ok,msg):
    if not ok: raise AssertionError(msg)

def trim(c):
    c=list(c)
    while len(c)>1 and not c[-1]: c.pop()
    return tuple(c)

def add(c,d):
    out=[0]*max(len(c),len(d))
    for i,v in enumerate(c): out[i]+=v
    for i,v in enumerate(d): out[i]+=v
    return trim(out)

def mul(c,d):
    out=[0]*(len(c)+len(d)-1)
    for i,v in enumerate(c):
        for j,w in enumerate(d): out[i+j]+=v*w
    return trim(out)

def at(c,x):
    v=0
    for b in reversed(c): v=v*x+b
    return v

def reflect(c):
    return trim([(-1)**i*sum(c[j]*comb(j,i) for j in range(i,len(c))) for i in range(len(c))])

def valid(S):
    return not S or (S[0]>=2 and all(b-a>=2 for a,b in zip(S,S[1:])))

def sets(M,start=2):
    yield ()
    for f in range(start,M+1):
        for tail in sets(M,f+2): yield (f,)+tail

@lru_cache(None)
def newton(S):
    if not S: return (1,)
    if not valid(S): return (0,)
    m=S[-1]; f=newton(S[:-1]); g=newton(S[:-1]+(m-1,)); out=[0]*m
    for j,b in enumerate(f): out[j]-=2*b
    for j,b in enumerate(g): out[j]-=b
    out[m-1]+=sum(b*comb(m-1,j) for j,b in enumerate(f))
    return trim(out)

def ivalue(S,n):
    return sum(b*comb(n,j) for j,b in enumerate(newton(S)) if j<=n)

def rvalue(S,x):
    out=Q(0); choose=Q(1)
    for j,b in enumerate(newton(S)):
        if j: choose*=(x-j+1)/j
        out+=b*choose
    return out

@lru_cache(None)
def falling(d):
    return (1,) if not d else mul(falling(d-1),(1-d,1))

def primitive(c):
    c=trim(c); g=reduce(gcd,c,0)
    need(g!=0,'zero polynomial')
    if c[-1]<0: g=-g
    return tuple(v//g for v in c)

def ipoly(S):
    b=newton(S); d=len(b)-1; out=[0]*(d+1)
    for j,v in enumerate(b):
        for i,w in enumerate(falling(j)): out[i]+=v*(factorial(d)//factorial(j))*w
    return primitive(out)

def divide(c,r):
    need(len(c)>1 and at(c,r)==0,'inexact division')
    out=[0]*(len(c)-1); out[-1]=c[-1]
    for j in range(len(out)-2,-1,-1): out[j]=c[j+1]+r*out[j+1]
    return trim(out)

@lru_cache(None)
def transform(S):
    m=S[-1] if S else 0
    return trim([(-1)**h*comb(m,h)*ivalue(S,m-h) for h in range(m+1)])

def extend(c,a,m):
    c=tuple(Q(v) for v in c)
    for n in range(a+1,m+1):
        out=[Q(0)]*(len(c)+1)
        for h,v in enumerate(c):
            w=Q(n,h+1)*v; out[0]+=w; out[h+1]-=w
        c=trim(out)
    need(all(v.denominator==1 for v in c),'nonintegral extension')
    return tuple(int(v) for v in c)

def moment(c,z):
    need(z>0,'moment domain')
    return sum((Q(v)/(z+h) for h,v in enumerate(c)),Q(0))

def derivative(S):
    m=S[-1]
    return sum((b*comb(m,j)*sum((Q(1,m-i) for i in range(j)),Q(0))
                for j,b in enumerate(newton(S))),Q(0))

def frac(x):
    x=Q(x); return [x.numerator,x.denominator]

def records_file(path,rows):
    h=hashlib.sha256(); n=0
    with path.open('wb') as raw:
        with gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0) as stream:
            for row in rows:
                b=(json.dumps(row,sort_keys=True,separators=(',',':'))+'\n').encode()
                stream.write(b); h.update(b); n+=1
    return dict(file=path.name,rows=n,uncompressed_sha256=h.hexdigest())

def brute(output):
    rows=[]
    for n in range(1,9):
        counts=Counter()
        for p in permutations(range(1,n+1)):
            S=tuple(i+1 for i in range(1,n-1) if p[i-1]<p[i]>p[i+1]); counts[S]+=1
        for S,num in counts.items(): need(ivalue(S,n)*2**(n-len(S)-1)==num,f'brute {n,S}')
        rows.append(dict(n=n,permutations=factorial(n),counts=[dict(S=list(S),count=v) for S,v in sorted(counts.items())]))
    (output/'permutation_counts.json').write_text(json.dumps(rows,indent=2)+'\n')
    return dict(n_range=[1,8],total_permutations=sum(factorial(n) for n in range(1,9)),file='permutation_counts.json')

def kernel(k,y,l1,l2):
    need(l1*l1+l2*l2==1 and (2*y).denominator==1,'kernel test inputs')
    t=l1*l1; u=l2*l2; b=Q(1); out=Q(0)
    for j in range(k-1):
        if j: b*=-(y-j)/j
        out+=(-1)**k*b*(t**j+u**j)/2
    if k%2==0: out-=l1**(int(2*y)-2)+l2**(int(2*y)-2)
    return out

def kernel_checks(K,output):
    triples=[(3,4,5),(5,12,13),(7,24,25),(20,21,29)]; zeros=0
    def rows():
        nonlocal zeros
        for k in range(2,K+1):
            for p in range(1,4*k+1):
                y=Q(p,2)
                for u,v,w in triples:
                    l1,l2=Q(u,w),Q(v,w); t=l1*l1; h=kernel(k,y,l1,l2)
                    need((h>0)-(h<0)==(y>k)-(y<k),f'kernel sign {k,y,t}')
                    zeros+=y==k
                    b=Q(1)
                    for j in range(k-1): b*=(y-1-j)/(j+1)
                    rhs=b/(2*k)*((y-k)*(t**k+(1-t)**k)-k*(t**(k-1)+(1-t)**(k-1)))
                    need(kernel(k+2,y,l1,l2)-h==rhs,'kernel two-step identity')
                    yield dict(k=k,y=frac(y),t=frac(t),A=frac(h))
    data=records_file(output/'kernel_exact.jsonl.gz',rows())
    return dict(k_range=[2,K],y_rule='p/2 for 1<=p<=4k',triples=triples,zero_cases=zeros,
                two_step_identity_cases=data['rows'],**data)

def poly_checks(M,output):
    import sympy as sp
    x=sp.Symbol('x'); nonempty=0; applicable=0; values=0
    def rows():
        nonlocal nonempty,applicable,values
        for S in sets(M):
            if not S: continue
            m=S[-1]; T=S[:-1]; a=T[-1] if T else 0; k=m-a
            D=transform(S); F=transform(T); q=extend(F,a,m)
            predicted=add((ivalue(T,m),),tuple(-v for v in add(q,reflect(q))))
            need(D==predicted and D==reflect(D),f'transform identity {S}')
            need(at(D,0)==at(D,1)==0 and at(D,Q(1,2))>0,f'transform boundaries {S}')
            need(ivalue(S,m)==0,f'root at maximum {S}')
            for y in sorted({Q(1,2),Q(1),Q(3,2),Q(k)-Q(1,2),Q(k),Q(k)+Q(1,2),Q(2*k)+Q(1,3)}):
                xx=a+y; pref=Q(1,factorial(a))
                for j in range(a+1): pref*=xx-j
                need(rvalue(T,xx)==pref*moment(F,y),f'Mellin {S,y}')
                integ=Q(0); choose=Q(1)
                for j in range(k-1):
                    if j: choose*=(y-j)/j
                    integ+=(-1)**(k-2-j)*choose*moment(F,Q(j+1))
                if k%2==0: integ-=2*moment(F,y)
                val=rvalue(S,xx)
                need(val==pref*integ and (val>0)-(val<0)==(y>k)-(y<k),f'integral/sign {S,y}')
                values+=1
            der=derivative(S); need(der>0,f'derivative {S}')
            dp=sp.Poly.from_list(list(reversed(D)),x,domain=sp.ZZ)
            while dp.eval(0)==0: dp=dp.exquo(sp.Poly(x,x))
            while dp.eval(1)==0: dp=dp.exquo(sp.Poly(x-1,x))
            dc=int(dp.count_roots(0,1)) if dp.degree()>0 else 0
            need(dc==0,f'transform roots {S}')
            c=ipoly(S); roots=[]
            for r in range(-3,m+1):
                while len(c)>1 and at(c,r)==0: c=divide(c,r); roots.append(r)
            need(roots.count(m)==1 and not [r for r in roots if r>a and r!=m],f'integer roots {S}')
            c=primitive(c); pp=sp.Poly.from_list(list(reversed(c)),x,domain=sp.ZZ)
            tail=int(pp.count_roots(a,sp.oo)) if pp.degree()>0 else 0
            need(tail==0,f'nonintegral last roots {S}')
            nonempty+=1; applicable+=len(S)>=2
            if nonempty%1000==0: print('polynomials',nonempty,flush=True)
            yield dict(S=list(S),newton=list(newton(S)),D=list(D),derivative_at_maximum=frac(der),
                       integer_roots=roots,residual_ascending=list(c),residual_roots_above_a=tail,
                       transform_interior_roots=dc)
    data=records_file(output/'polynomials_exact.jsonl.gz',rows())
    return dict(maximum=M,nonempty=nonempty,conjecture_6_6_cases=applicable,
                singleton_cases=nonempty-applicable,exact_moment_and_sign_tests=values,
                transform_identity_cases=nonempty,exact_transform_root_counts=nonempty,
                exact_last_root_counts=nonempty,positive_derivative_checks=nonempty,
                sympy=sp.__version__,**data)

def input_audit(path):
    seen=set(); applicable=0
    with zipfile.ZipFile(path) as z:
        need(z.testzip() is None,'input CRC')
        manifest=json.loads(z.read('INPUT_MANIFEST.json'))
        for row in manifest['files']:
            b=z.read(row['path'])
            need(len(b)==row['bytes'] and hashlib.sha256(b).hexdigest()==row['sha256'],'input manifest')
        raw=gzip.decompress(z.read('exploration/results/peak_roots_max_24.jsonl.gz'))
        summary=json.loads(z.read('exploration/results/peak_summary_max_24.json'))
        need(hashlib.sha256(raw).hexdigest()==summary['uncompressed_jsonl_sha256'],'input data hash')
        for line in raw.splitlines():
            row=json.loads(line); S=tuple(row['S'])
            need(S and valid(S) and S[-1]<=24 and S not in seen,'input set')
            seen.add(S); c=tuple(row['residual_ascending_coefficients'])
            for r in row['integer_roots']: c=mul(c,(-r,1))
            need(primitive(c)==ipoly(S),f'input polynomial {S}')
            applicable+=len(S)>=2
        need(seen==set(sets(24))-{()},'input enumeration')
    return dict(input_zip_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),manifest_entries_checked=len(manifest['files']),
                algebraic_records_checked=len(seen),applicable_records=applicable,
                note='Input polynomial identities and enumeration checked; input root counts above the fresh bound NOT rerun.')

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--max',type=int,default=18); p.add_argument('--kernel-max',type=int,default=24)
    p.add_argument('--output',type=Path,default=Path('data')); p.add_argument('--input-zip',type=Path)
    args=p.parse_args(); need(args.max>=2 and args.kernel_max>=2,'invalid bounds'); args.output.mkdir(parents=True,exist_ok=True)
    result=dict(python=platform.python_version(),arithmetic='exact integer/rational, exact SymPy real-root counts')
    result['permutations']=brute(args.output)
    result['kernel']=kernel_checks(args.kernel_max,args.output)
    result['polynomials']=poly_checks(args.max,args.output)
    if args.input_zip: result['input_audit']=input_audit(args.input_zip)
    result['failures']=[]
    (args.output/'summary.json').write_text(json.dumps(result,indent=2)+'\n'); print(json.dumps(result,indent=2))

if __name__=='__main__': main()
