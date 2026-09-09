\\ ============================================================
\\ cusp_values_complete.gp
\\
\\ 按文中“计算 φ 在尖点处的精确值”的算法实现（PARI/GP）
\\
\\ 依赖：
\\   - 已有椭圆曲线 E = ellinit("...")
\\   - 已有 F_N(x,j), G_N(y,j) 的 GP 表达式，例如：
\\       F = fxj(Xsym, Jsym);
\\       G = fyj(Ysym, Jsym);
\\
\\ 主要入口：
\\   ZR = exact_cusp_values(N, E, F, G, Xsym, Ysym, Jsym);
\\   print_exact_cusp_values(ZR);
\\
\\ 若你的脚本中函数名就是 fxj/fyj，可用：
\\   ZR = run_fxj_fyj_cusp_values(N, E);
\\ ============================================================

default(realprecision, 80);

\\ -------------------- 基本数值工具 --------------------

close_complex(z1, z2, tol) = norml2(z1 - z2) < tol^2;

numeric_clean(z) =
{
  my(w = z, prec = getlocalprec());
  while(type(w) == "t_SER" || type(w) == "t_POL",
    w = subst(w, variable(w), 0);
  );
  if(type(w) == "t_INT" || type(w) == "t_FRAC" || type(w) == "t_REAL", return(precision(w + 0.*I, prec)));
  if(type(w) == "t_COMPLEX", return(precision(w, prec)));
  \\ 对少数代数/符号污染情形，强制进入复数数值域
  precision(real(w) + imag(w)*I, prec)
}

close_point(P, Q, tol) =
{
  if(#P != #Q, return(0));
  if(#P == 1, return(1));
  close_complex(P[1], Q[1], tol) && close_complex(P[2], Q[2], tol)
}

point_on_curve(E, x, y, tol = 1e-20) =
{
  my(lhs = y^2 + E.a1*x*y + E.a3*y,
     rhs = x^3 + E.a2*x^2 + E.a4*x + E.a6);
  norml2(lhs - rhs) < tol^2
}

reduce_pair(a, b) =
{
  if(b == 0, return([1, 0]));
  my(g = gcd(a, b));
  a /= g; b /= g;
  if(b < 0, a = -a; b = -b);
  [a, b]
}

frac_to_pair(q) =
{
  if(type(q) == "t_INT", return(reduce_pair(q, 1)));
  if(type(q) == "t_FRAC", return(reduce_pair(numerator(q), denominator(q))));
  error("frac_to_pair: unsupported cusp object")
}

pair_to_frac(v) =
{
  if(v[2] == 0, return("oo"));
  if(v[2] == 1, return(v[1]));
  v[1]/v[2]
}

mat_inverse_sl2(M) = [M[2,2], -M[1,2]; -M[2,1], M[1,1]];

mobius_pair(M, v) =
{
  my(a = v[1], b = v[2]);
  if(b == 0,
    if(M[2,1] == 0, return([M[1,1], 0]));
    return(reduce_pair(M[1,1], M[2,1]));
  );
  reduce_pair(M[1,1]*a + M[1,2]*b, M[2,1]*a + M[2,2]*b)
}

cusp_lift_matrix(a, c) =
{
  my(v = reduce_pair(a, c));
  a = v[1]; c = v[2];
  if(c == 0, return([1, 0; 0, 1]));
  my(uvg = gcdext(a, c));
  my(u = uvg[1], vv = uvg[2], g = uvg[3]);
  if(g != 1, error("cusp_lift_matrix: pair not primitive"));
  [a, -vv; c, u]   \\ determinant = a*u + c*vv = 1
}

solve_linear_congruence(a, b, m) =
{
  my(g = gcd(a, m));
  if(Mod(b, g) != 0, error("solve_linear_congruence: no solution"));
  my(a1 = a/g, b1 = b/g, m1 = m/g);
  if(m1 == 1, return(0));
  lift(Mod(b1, m1) / Mod(a1, m1))
}

\\ -------------------- 候选精确值集合 Z --------------------
\\ Z 的元素记为 [x_i, y_j, ix, iy]

vec_contains_point(V, P, tol = 1e-20) =
{
  for(i = 1, #V,
    if(close_point(V[i][1..2], P, tol), return(1));
  );
  0
}

unique_roots(R, tol = 1e-18) =
{
  my(L = List());
  for(i = 1, #R,
    my(hit = 0);
    for(j = 1, #L,
      if(close_complex(L[j], R[i], tol), hit = 1; break);
    );
    if(!hit, listput(~L, R[i]));
  );
  Vec(L)
}

leading_coeff_in_var(P, V) = polcoef(P, poldegree(P, V), V);

candidate_set_Z(E, F, G, Xvar, Yvar, Jvar, root_tol = 1e-18, curve_tol = 1e-16) =
{
  my(B = leading_coeff_in_var(F, Jvar));
  my(C = leading_coeff_in_var(G, Jvar));

  my(xroots = [], yroots = [], Z = []);

  if(type(B) == "t_POL", xroots = unique_roots(polroots(B), root_tol));
  if(type(C) == "t_POL", yroots = unique_roots(polroots(C), root_tol));

  for(ix = 1, #xroots,
    for(iy = 1, #yroots,
      if(point_on_curve(E, xroots[ix], yroots[iy], curve_tol),
        my(P = [xroots[ix], yroots[iy], ix, iy]);
        if(!vec_contains_point(Z, P, curve_tol),
          Z = concat(Z, [P]);
        );
      );
    );
  );

  [B, C, xroots, yroots, Z]
}

\\ -------------------- Eichler 积分与周期 --------------------

init_eichler_symbol(E) =
{
  my(mff = mffromell(E));
  my(mf = mff[1], f = mff[2]);
  [mf, f, mfsymbol(mf, f)]
}

gamma_from_symbol(fs, tau) = numeric_clean(2*Pi*I * mfsymboleval(fs, [oo, numeric_clean(tau)]));

period_of_matrix(fs, M, tau0 = 1.0*I) =
{
  gamma_from_symbol(fs, (M[1,1]*tau0 + M[1,2])/(M[2,1]*tau0 + M[2,2]))
  - gamma_from_symbol(fs, tau0)
}

nearest_period_data(E, z) =
{
  my(om = ellperiods(E), w1 = om[1], w2 = om[2]);
  my(det = real(w1)*imag(w2) - imag(w1)*real(w2));
  my(a = ( real(z)*imag(w2) - imag(z)*real(w2)) / det);
  my(b = (-real(z)*imag(w1) + imag(z)*real(w1)) / det);
  my(m = round(a), n = round(b));
  my(p = m*w1 + n*w2);
  [m, n, p, z - p]
}

\\ -------------------- 尖点等价与构造 M_j --------------------
\\ 使用 Cremona 2.2.3 的同余判别的可计算版：
\\ 只需枚举 y mod N, gcd(y,N)=1, 检查
\\   y*c' ≡ c (mod N),  且 gcd(c,N) | (y*a' - a)

cusp_equiv(N, v1, v2) =
{
  my(a = reduce_pair(v1[1], v1[2])[1], c = reduce_pair(v1[1], v1[2])[2]);
  my(ap = reduce_pair(v2[1], v2[2])[1], cp = reduce_pair(v2[1], v2[2])[2]);

  \\ 将 oo 类替换成 1/N 类，便于在 Γ0(N) 下比较
  if(c == 0, a = 1; c = N);
  if(cp == 0, ap = 1; cp = N);

  my(d = gcd(c, N));
  if(gcd(cp, N) != d, return(0));

  for(y = 0, N-1,
    if(gcd(y, N) == 1,
      if(Mod(y*cp - c, N) == 0,
        if(Mod(y*ap - a, d) == 0, return(1));
      );
    );
  );
  0
}

find_aux_prime_for_cusp(N, cusp, pbound = 1000) =
{
  my(v = reduce_pair(cusp[1], cusp[2]));
  my(s = v[1], r = v[2]);

  forprime(p = 2, pbound,
    if(gcd(p, N) == 1,
      my(ok = 1);
      if(!cusp_equiv(N, [s, r], [p*s, r]), ok = 0);
      if(ok,
        for(j = 0, p-1,
          if(!cusp_equiv(N, [s, r], [r*j + s, p*r]), ok = 0; break);
        );
      );
      if(ok, return(p));
    );
  );
  error("find_aux_prime_for_cusp: no prime found up to pbound")
}

construct_gamma0_matrix_between_cusps(N, src, dst) =
{
  my(vs = reduce_pair(src[1], src[2]));
  my(vt = reduce_pair(dst[1], dst[2]));
  my(u = vs[1], v = vs[2], s = vt[1], r = vt[2]);

  my(Au = cusp_lift_matrix(u, v));
  my(At = cusp_lift_matrix(s, r));

  my(du = Au[2,2], dt = At[2,2]);
  my(rhs = r*du - dt*v);
  my(coeff = r*v);
  my(n = solve_linear_congruence(coeff, rhs, N));
  my(Tn = [1, n; 0, 1]);
  my(M = At * Tn * mat_inverse_sl2(Au));

  if(Mod(M[2,1], N) != 0,
    error("construct_gamma0_matrix_between_cusps: lower-left entry not divisible by N")
  );
  M
}

cusp_sources_for_prime(cusp, p) =
{
  my(v = reduce_pair(cusp[1], cusp[2]));
  my(s = v[1], r = v[2]);
  my(V = []);
  for(j = 0, p-1,
    V = concat(V, [reduce_pair(r*j + s, p*r)]);
  );
  V = concat(V, [reduce_pair(p*s, r)]);
  V
}

\\ -------------------- 尖点的近似像与精确匹配 --------------------

cusp_gamma_value(N, E, fs, cusp, pbound = 1000, tau0 = 1.0*I) =
{
  my(v = reduce_pair(cusp[1], cusp[2]));
  my(s = v[1], r = v[2]);
  my(p = find_aux_prime_for_cusp(N, [s, r], pbound));
  my(srcs = cusp_sources_for_prime([s, r], p));
  my(sumP = 0.0 + 0.0*I);

  for(j = 1, #srcs,
    my(M = construct_gamma0_matrix_between_cusps(N, srcs[j], [s, r]));
    sumP += period_of_matrix(fs, M, tau0);
  );

  my(ap = ellak(E, p));
  [p, srcs, sumP/(p + 1 - ap)]
}

cusp_image_approx(N, E, fs, cusp, pbound = 1000, tau0 = 1.0*I, lattice_tol = 1e-30) =
{
  my(CG = cusp_gamma_value(N, E, fs, cusp, pbound, tau0));
  my(gam_raw = numeric_clean(CG[3]));

  \\ 先对周期格做归约
  my(PD = nearest_period_data(E, gam_raw));
  my(m = PD[1], n = PD[2], per = PD[3], rem = numeric_clean(PD[4]));

  \\ 如果已经极接近周期格点，则直接判为无穷远点 O
  if(norml2(rem) < lattice_tol^2,
    return([CG[1], CG[2], gam_raw, [0], [m, n, per, rem]]);
  );

  \\ 否则再送入 ellztopoint
  my(P = ellztopoint(E, rem));
  [CG[1], CG[2], gam_raw, P, [m, n, per, rem]]
}

match_point_to_Z(P, Z, tol = 1e-8) =
{
  if(#P == 1, return([0, "O"]));
  for(i = 1, #Z,
    if(close_point(P, Z[i][1..2], tol),
      return([i, Z[i]]);
    );
  );
  [0, "UNMATCHED"]
}

\\ -------------------- 主程序 --------------------
\\ 返回 [Zdata, results]
\\ 其中
\\   Zdata = [B(X), C(Y), xroots, yroots, Z]
\\   results 的每个元素为
\\   [cusp_pair, aux_prime_p, gamma_value, approx_point, exact_match]

exact_cusp_values(N, E, F, G, Xvar, Yvar, Jvar,pbound = 1000, tau0 = 1.0*I,root_tol = 1e-18, curve_tol = 1e-16, match_tol = 1e-8) =
{
  my(Zdata = candidate_set_Z(E, F, G, Xvar, Yvar, Jvar, root_tol, curve_tol));
  my(Z = Zdata[5]);
  my(es = init_eichler_symbol(E));
  my(fs = es[3]);
  my(cusps = mfcusps(N));
  my(R = []);

  print("#cusps = ", #cusps);
  print("#candidate finite values in Z = ", #Z);

  for(i = 1, #cusps,
    my(cpair = frac_to_pair(cusps[i]));
    my(A = cusp_image_approx(N, E, fs, cpair, pbound, tau0, 1e-30));
    my(p = A[1], gam = A[3], P = A[4], red = A[5]);
    my(ex = match_point_to_Z(P, Z, match_tol));
    R = concat(R, [[cpair, p, gam, P, ex]]);
  );
  [Zdata, R]
}

\\ -------------------- 便捷入口 --------------------

run_fxj_fyj_cusp_values(N, E,pbound = 1000, tau0 = 1.0*I,root_tol = 1e-18, curve_tol = 1e-16, match_tol = 1e-8) =
{
  my(Xsym = 'X, Ysym = 'Y, Jsym = 'J);
  my(F = fxj(Xsym, Jsym));
  my(G = fyj(Ysym, Jsym));
  exact_cusp_values(N, E, F, G, Xsym, Ysym, Jsym,
                    pbound, tau0, root_tol, curve_tol, match_tol)
}

\\ -------------------- 打印 --------------------

print_candidate_set(Zdata) =
{
  my(B = Zdata[1], C = Zdata[2], xr = Zdata[3], yr = Zdata[4], Z = Zdata[5]);
  print("B(X) = leading coeff of F(X,J) in J");
  print(B);
  print("C(Y) = leading coeff of G(Y,J) in J");
  print(C);
  print("#roots of B = ", #xr);
  print("#roots of C = ", #yr);
  print("#candidate points Z = ", #Z);
  for(i = 1, #Z,
    print("Z[", i, "] = (", Z[i][1], ", ", Z[i][2], ")  [ix=", Z[i][3], ", iy=", Z[i][4], "]");
  );
}

print_exact_cusp_values(ZR) =
{
  my(Zdata = ZR[1], R = ZR[2]);
  print_candidate_set(Zdata);
  print("========================================");
  print("Cusp values:");
  for(i = 1, #R,
    print("----------------------------------------");
    print("#", i);
    print("cusp rep       = [", R[i][1][1], "/", R[i][1][2], "]");
    print("aux prime p    = ", R[i][2]);
    print("gamma(cusp)    = ", R[i][3]);
    if(#R[i] >= 6,
      print("reduced gamma  = ", R[i][6][4]);
      print("nearest period = ", R[i][6][3]);
      );
    print("approx point   = ", R[i][4]);
    if(type(R[i][5][2]) == "t_STR",
      print("exact value     = ", R[i][5][2]);
    ,
      print("exact value idx = Z[", R[i][5][1], "]");
      print("exact value     = (", R[i][5][2][1], ", ", R[i][5][2][2], ")");
    );
  );
}
