mu(N)=
{my(A,N1,s,n);
 A=factor(N);N1=matsize(A)[1];
s=1;for(n=1,N1,s=s*(1+1/A[n,1]));
return(N*s);
}
fxj(N,x,j) = 
{
  my(S1,T1,S2,M,u,v,X,t1,m,n,r,v1);
  my(J,c,c1,N1,tp,tp1);
e = ellinit(concat(N,"a1"));
S1=mu(N)+1;T1=2*ellmoddegree(e)+1;S2=S1*T1;
M=S2+50;
[u,v]=elltaniyama(e,M+1);X=vector(S1);
X[1]=1;for(n=2,S1,X[n]=X[n-1]*u);
J=vector(T1);v1=laurentseries(ellj,M);
J[1]=1;for(n=2,T1,J[n]=J[n-1]*v1);
c=matrix(M+3,S2);
for(m=1,S1,for(n=1,T1,t1=Vec(X[m]*J[n]);\
N1=length(t1);tp=2*(S1-m)+(T1-n);\
tp1=T1*(m-1)+n;for(r=1,N1-tp,c[r+tp,tp1] = t1[r])));
c[2*(S1-1)+T1,1]=1;
rk=matrank(c);
c=c[1..M-10,];
c1=matker(c);
U=matsize(c1)[2];
if(U!=0, out=vector(U);tp2=concat("C://fxjj",N);\
 for(n=1,U,out[n]=concat(concat(tp2,n),".c"));\
for(i=1,U,nw = fileopen(out[i], "w");\
s1=0; for(m=1,S1,for(n=1,T1,s1=s1+c1[T1*(m-1)+n,i]*x^(m-1)*j^(n-1)));\
  filewrite(nw,s1); fileclose(nw)));
return(U);
}

