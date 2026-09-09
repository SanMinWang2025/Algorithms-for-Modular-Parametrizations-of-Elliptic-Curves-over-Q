\\ certify_fxj389.gp
\\
\\ Purpose:
\\   Certify the identity
\\        F_389(x(q), j(q)) = 0
\\   by the finite q-expansion certificate used in the paper.
\\
\\ Input:
\\   The file P:/fxj389.gp must define a PARI/GP function
\\        fxj(x,j)
\\   for the candidate polynomial F_389(x,j).
\\
\\ Typical usage in GP:
\\   default(parisize,"7000G");       \\ adjust to the available machine
\\   read("P:/certify_fxj389.gp");
\\   R = certify_F389_fxj("P:/fxj389.gp");
\\
\\ If read("P:/fxj389.gp") gives "expression nested too deeply", rewrite
\\ fxj389.gp in a shallow form, for example as many lines of the form
\\        F += c*x^k*j^l;
\\ rather than as one enormous nested expression.

cert389_j_series(qprec)=
{
  my(E4,E6);
  E4 = Ser(vector(qprec+2, n, if(n==1, 1, 240*sigma(n-1,3))), q);
  E6 = Ser(vector(qprec+2, n, if(n==1, 1, -504*sigma(n-1,5))), q);
  return(1728*E4^3/(E4^3-E6^2));
};

cert389_x_series(qprec)=
{
  my(E,T,Xq);

  \\ 389a1 in the model y^2 + y = x^3 + x^2 - 2*x.
  E = ellinit([0,1,1,-2,0]);
  T = elltaniyama(E, qprec+2);
  Xq = subst(T[1], variable(T[1]), q);
  return(Xq);
};

cert389_report_series(Xq,Jq)=
{
  print("valuation(x(q)) = ", valuation(Xq,q));
  print("precision(x(q)) = ", serprec(Xq,q));
  print("valuation(j(q)) = ", valuation(Jq,q));
  print("precision(j(q)) = ", serprec(Jq,q));
};

certify_F389_fxj(polyfile="P:/fxj389.gp", qprec=63360)=
{
  my(N=389,d=40,mu=390,K=390,L=80,BF,minexp,Xq,Jq,R,v,sp);

  BF = 2*d*K + mu*L;
  minexp = -2*K - L;

  print("Reading candidate polynomial from: ", polyfile);
  read(polyfile);

  print("N = ", N);
  print("K = ", K, ", L = ", L, ", d = ", d, ", mu = ", mu);
  print("lowest possible exponent = ", minexp);
  print("finite certificate bound B_F = ", BF);
  print("requested input q-precision = ", qprec);

  Xq = cert389_x_series(qprec);
  Jq = cert389_j_series(qprec);
  cert389_report_series(Xq,Jq);

  if(valuation(Xq,q) != -2,
    error("Unexpected pole order for x(q).  Expected valuation -2.")
  );
  if(valuation(Jq,q) != -1,
    error("Unexpected pole order for j(q).  Expected valuation -1.")
  );

  print("Evaluating fxj(x(q),j(q)) ...");
  R = fxj(Xq,Jq);

  if(R == 0,
    print("Residual is exactly 0 as a PARI object.");
    print("F_389 certificate: PASS");
    return([1, BF, qprec, "exact zero"]);
  );

  v = valuation(R,q);
  sp = serprec(R,q);

  print("residual valuation = ", v);
  print("residual precision = ", sp);

  if(v > BF,
    print("F_389 certificate: PASS");
    print("All coefficients from q^", minexp, " through q^", BF, " vanish.");
    return([1, BF, qprec, v, sp]);
  );

  print("F_389 certificate: FAIL or insufficient precision.");
  print("The first visible residual term is at exponent ", v, ".");
  print("To certify the identity, the residual valuation must be > ", BF, ".");
  if(sp <= BF,
    print("The residual precision is not beyond B_F; increase qprec.")
  );
  return([0, BF, qprec, v, sp]);
};

