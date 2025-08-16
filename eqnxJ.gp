mu(N)=
{my(A,N1,s,n);
 A=factor(N);N1=matsize(A)[1];
s=1;for(n=1,N1,s=s*(1+1/A[n,1]));
return(N*s);
}
fxJ(N,x,j) = 
{
  my(e,S,T,S1,T1,S2,M,M1,u,v,X,t1,m,n,r,v1);
  my(J,c,c1,N1,N2,tp,tp1,U,out,tp2,nw,i);
e = ellinit(concat(N,"a1"));S=mu(N);T=2*ellmoddegree(e);
S1=S+1;T1=T+1;S2=S1*T1;
M=max(S2+floor(S2/1),2*S+N*T+1);M1=floor(M/N);
[u,v]=elltaniyama(e,M+1);X=vector(S1);
X[1]=1;for(n=2,S1,X[n]=X[n-1]*u);
J=vector(T1);
v0=laurentseries(ellj,M1);v1=substpol(v0, x, x^N);
J[1]=1;for(n=2,T1,J[n]=J[n-1]*v1);
c=matrix(M,S2);
for(m=1,S1,for(n=1,T1,t1=Vec(X[m]*J[n]); \
N1=length(t1);tp=2*(S1-m)+N*(T1-n);N2=min(N1-tp,M-tp);\
tp1=T1*(m-1)+n;for(r=1,N2,c[r+tp,tp1] = t1[r])));
c[2*S+N*T+1,1]=1;
c1=matker(c,1);
U=matsize(c1)[2];
if(U!=0,out=vector(U);tp2=concat("P://fxJ",N);\
 for(n=1,U,out[n]=concat(concat(tp2,n),".c"));\
for(i=1,U,nw = fileopen(out[i], "w");s1=0;\
 for(m=1,S1,for(n=1,T1,s1=s1+c1[T1*(m-1)+n,i]*x^(m-1)*j^(n-1)));\
  filewrite(nw,s1); fileclose(nw)));
return(U);
}