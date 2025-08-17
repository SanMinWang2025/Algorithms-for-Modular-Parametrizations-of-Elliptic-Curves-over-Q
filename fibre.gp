
jinv(j)=
{my(e,a,b,tp);
if(j==1728,return(I));
if(j==0,return(-1/2+I*sqrt(3)/2));
e=ellinit([1,0,0,-36*(j-1728)^-1,-(j-1728)^-1]);
a=e.omega[1];b=e.omega[2];tp=a/b;
if(imag(tp)<0,tp=b/a);
return(tp);
}

fibre0(x0,N)={my(j0,N2,j1,num,n,z1,op,m,
j20,N3,j2,w,N1,w1,z0,z2,k);
j0=polroots(fxj(x0,j));N2=length(j0);
j1=vector(N2);num=0;
for(n=1,N2,z1=j0[n];op=0;\
for(m=1,num,if(abs(j1[m]-z1)<10^-10,op++));\
if(op==0,num++;j1[num]=z1));
j1=j1[1..num];N2=num;
j20=polroots(fxJ(x0,j));N3=length(j20);
j2=vector(N3);num=0;
for(n=1,N3,z1=j20[n];op=0;\
for(m=1,num,if(abs(j2[m]-z1)<10^-10,op++));\
if(op==0,num++;j2[num]=z1));
j2=j2[1..num];N3=num;
w=mfcosets(N);N1=length(w);
w1=vector(2*N2);num=0;
for(n=1,N2,z0=jinv(j1[n]);\
for(m=1,N1,z1=(z0*w[m][1,1]+w[m][1,2])/(z0*w[m][2,1]+w[m][2,2]);\
z2=ellj(N*z1);\
for(k=1,N3,if(abs(j2[k]-z2)<10^-10,num++;w1[num]=[n,m,z1]))));
return(w1[1..num]);
}

ph(tau,v)=
{my(s,q,q1,n);
s= 0.;
q=exp(2*I*Pi*tau); 
q1 = 1.; for (n = 1, 1500000, q1 *= q; s += v[n]/n*q1);
return(s);
}

fibre(x0,y0,ell,N)={my(X,N1,w1,e,num,n,z0,z1,P);
X=fibre0(x0,N);
N1=length(X);w1=vector(N1);
e = ellinit(ell);
v=ellan(e,1500000);
num=0;
for(n=1,N1,z0=X[n][3];z1=ph(z0,v);P=ellztopoint(e,z1);
if(abs(P[1]-x0)<10^-10&&abs(P[2]-y0)<10^-10,
num++;w1[num]=z0));
return(w1[1..num]);
}
