compa(x, y, N0) = {my(p1,q1,p2,q2,u,v,d,s1,r1,s2,r2,tp,tp22,tp2);
  p1 = numerator(x); 
  q1 = denominator(x); 
  p2 = numerator(y); 
  q2 = denominator(y);
  [u,v,d]=gcdext(p1,q1);
   s1 = u;r1 = -v; 
  [u,v,d]=gcdext(p2,q2);
   s2 = u;r2 = -v; 
   tp = 0;
 tp1 = gcd(q1*q2, N0); tp22 =q2*s1 - q1*s2;
if(tp22%tp1==0,tp=1);
 if (tp == 1, tp2 =tp22/ tp1;
[u,v,d]=gcdext((q1*q2) / tp1,(N0 / tp1));
u=-u*tp2;
v=v*tp2;
  s1 = s1 + q1*u;
 r1 = r1 + p1*u;
  return([tp, [p1, r1;q1, s1]*[p2, r2;q2, s2]^-1]));
 return([tp,0]);
}
compa0(x, N0) = {
    my(n = 0, op = 1, p, A, j, tmp);
    while(op > 0,
        n++;
        p = prime(n); 
        A = vector(p + 1); 
        op = 0;
              for(j = 1, p,
            tmp = compa(x, (x + j)/p, N0);  
            if(tmp[1] == 0,
                op++;
                break,
            A[j] = tmp[2]));   
            tmp = compa(x, x*p, N0);
        if(tmp[1] == 1,
            A[p + 1] = tmp[2],
            op++
        );
    );
    return(A); 
};

ph(tau,v)=
{my(s,q,q1,n);
s= 0.;
q=exp(2*I*Pi*tau); 
q1 = 1.; for (n = 1, 1500000, q1 *= q; s += v[n]/n*q1);
return(s);
}

cuspv(ell,N0)={my(B,N1,w0,n,v,N2,w2,w1,A,c,d,z0,z1,w,s);
e = ellinit(ell);
B=mfcusps(N0);N1=length(B);
w0=vector(N1);for(n=1,N1,w0[n]=compa0(B[n],N0));
v=ellan(e,1500000);
N2=length(w0);
w2=vector(N2);
for(n=1,N2,w=w0[n];N1=length(w);w1=vector(N1);A=vector(N1);\
for(m =1,N1,c=w[m][2,1];d=w[m][2,2];if(c!=0,z0=-d/c+1/abs(c)*I,z0=1+I);\
z1=(z0*w[m][1,1]+w[m][1,2])/(z0*w[m][2,1]+w[m][2,2]);A[m]=ph(z1,v)-ph(z0,v));\
s=0;for(m=1,N1,s=s+A[m]);w2[n]=[n,B[n],ellztopoint(e,s/(N1-v[N1-1])),s/(N1-v[N1-1])]);
return(w2);
}