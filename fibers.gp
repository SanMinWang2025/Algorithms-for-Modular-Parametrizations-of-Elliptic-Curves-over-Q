\\ ============================================================
\\ fiber_complete.gp
\\
\\ 计算 \varphi 的纤维（PARI/GP）
\\
\\ 主要入口：
\\   1. 通用入口（含尖点步骤）
\\      R = fiber_of_phi_with_cusps(N, E, P, F, H, G, Xsym, Ysym, jsym, Jsym);
\\      其中
\\         F = F_N(X,j),
\\         H = f_N(X,J),
\\         G = G_N(Y,j).
\\
\\   2. 若你的脚本名就是 fxj, fxJ, fyj，可直接用
\\      R = run_fxj_fxJ_fyj_fiber(N, E, P);
\\
\\   3. 若暂时不想做尖点补充，可用
\\      R = run_fxj_fxJ_fiber(N, E, P);
\\
\\ 约定：
\\   - 有限点 P 以 [alpha,beta] 传入；
\\   - 无穷远点 O 以 [0] 传入。
\\ ============================================================

default(realprecision, 80);

\\ -------------------- 基本工具 --------------------

numeric_clean(z) =
{
  my(w = z, prec = getlocalprec());
  while(type(w) == "t_SER" || type(w) == "t_POL",
    w = subst(w, variable(w), 0);
  );
  if(type(w) == "t_COMPLEX", return(precision(w, prec)));
  if(type(w) == "t_INT" || type(w) == "t_FRAC" || type(w) == "t_REAL",
    return(precision(w + 0.*I, prec))
  );
  precision(real(w) + imag(w)*I, prec)
}

safe_real_tol(t, dflt = 1e-15) =
{
  my(u = t);

  if(type(u) == "t_SER" || type(u) == "t_POL",
    u = subst(u, variable(u), 0);
  );

  if(type(u) == "t_INT" || type(u) == "t_FRAC" || type(u) == "t_REAL",
    return(abs(u + 0.0));
  );

  if(type(u) == "t_COMPLEX",
    return(abs(real(u) + 0.0));
  );

  dflt
}

cabs2(z) = norml2(numeric_clean(z));

close_complex(z1, z2, tol) =
{
  my(w1 = numeric_clean(z1), w2 = numeric_clean(z2));
  my(sc = max(max(1.0, sqrt(norml2(w1))), sqrt(norml2(w2))));
  norml2(w1 - w2) < (tol*sc)^2
}

point_close(P, Q, tol) =
{
  if(#P == 1 && #Q == 1, return(1));
  if(#P == 1 || #Q == 1, return(0));
  close_complex(P[1], Q[1], tol) && close_complex(P[2], Q[2], tol)
}

point_on_curve(E, x, y, tol = 1e-18) =
{
  my(lhs = y^2 + E.a1*x*y + E.a3*y,
     rhs = x^3 + E.a2*x^2 + E.a4*x + E.a6);
  close_complex(lhs, rhs, tol)
}

unique_roots(R, tol = 1e-18) =
{
  my(L = List());
  for(i = 1, #R,
    my(hit = 0);
    for(j = 1, #L,
      if(close_complex(L[j], R[i], tol), hit = 1; break);
    );
    if(!hit, listput(~L, numeric_clean(R[i])));
  );
  Vec(L)
}

root_in_set(z, V, tol = 1e-8) =
{
  for(i = 1, #V,
    if(close_complex(z, V[i], tol), return(1));
  );
  0
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

mobius(M, z) = (M[1,1]*z + M[1,2]) / (M[2,1]*z + M[2,2]);

mat_inverse_sl2(M) = [M[2,2], -M[1,2]; -M[2,1], M[1,1]];

matrix_in_gamma0(N, M) =
{
  if(type(M) != "t_MAT", return(0));
  if(matsize(M)[1] != 2 || matsize(M)[2] != 2, return(0));
  if(matdet(M) != 1 && matdet(M) != -1, return(0));
  if(Mod(M[2,1], N) != 0, return(0));
  1
}

stabilizer_matrices_from_j(j0, rooteps = 1e-18) =
{
  my(I2 = [1, 0; 0, 1]);
  my(S  = [0, -1; 1, 0]);
  my(ST = [0, -1; 1, 1]);

  if(close_complex(j0, 1728.0, rooteps), return([I2, S]));
  if(close_complex(j0, 0.0, rooteps), return([I2, ST, ST^2]));
  [I2]
}

same_x0_point_via_stabilizer(N, j0, M1, M2, rooteps = 1e-18) =
{
  my(St = stabilizer_matrices_from_j(j0, rooteps));
  my(M2inv = mat_inverse_sl2(M2));

  for(k = 1, #St,
    my(A = M1 * St[k] * M2inv);
    if(matrix_in_gamma0(N, A), return(1));
  );
  0
}

\\ -------------------- 右陪集代表 Γ0(N) \ SL2(Z) --------------------

vec_contains_exact(V, v) =
{
  for(i = 1, #V, if(V[i] == v, return(1)));
  0
}

canon_pair(c, d, N) =
{
  my(best = [N+1, N+1], got = 0, cu, du);
  for(u = 0, N-1,
    if(gcd(u, N) == 1,
      cu = lift(Mod(u*c, N));
      du = lift(Mod(u*d, N));
      if(!got || cu < best[1] || (cu == best[1] && du < best[2]),
        best = [cu, du];
        got = 1;
      );
    );
  );
  best
}

P1ZN(N) =
{
  my(V = []);
  for(c = 0, N-1,
    for(d = 0, N-1,
      if(gcd(gcd(c,d), N) == 1,
        my(rep = canon_pair(c, d, N));
        if(!vec_contains_exact(V, rep), V = concat(V, [rep]));
      );
    );
  );
  V
}

lift_primitive_pair(c0, d0, N) =
{
  my(t = 0);
  while(1,
    if(gcd(c0, d0 + t*N) == 1, return([c0, d0 + t*N]));
    if(gcd(c0 + t*N, d0) == 1, return([c0 + t*N, d0]));
    t++;
    if(t > 20*N, error("lift_primitive_pair failed"));
  );
}

pair_to_matrix(c, d, N) =
{
  my(v = lift_primitive_pair(c, d, N));
  my(c1 = v[1], d1 = v[2]);
  my(guv = gcdext(c1, d1));
  my(u = guv[1], vv = guv[2], g = guv[3]);
  if(g != 1, error("pair_to_matrix: pair not primitive"));
  [vv, -u; c1, d1]
}

gamma0_right_cosets(N) =
{
  my(P = P1ZN(N), reps = []);
  for(i = 1, #P,
    reps = concat(reps, [pair_to_matrix(P[i][1], P[i][2], N)]);
  );
  reps
}

\\ -------------------- 从 j 恢复 tau --------------------

curve_from_j_model(j0) =
{
  ellinit([1, 0, 0, -36/(j0 - 1728), -1/(j0 - 1728)])
}

tau_from_j(j0, rooteps = 1e-18) =
{
  my(z = numeric_clean(j0));
  if(close_complex(z, 1728.0, rooteps), return(1.0*I));
  if(close_complex(z, 0.0, rooteps), return((-1 + sqrt(-3.0))/2));
  my(Ej = curve_from_j_model(z));
  my(W = ellperiods(Ej));
  my(t1 = numeric_clean(W[1] / W[2]));
  if(imag(t1) > 0, return(t1));
  numeric_clean(W[2] / W[1])
}

\\ -------------------- Eichler 积分与 \varphi(tau) --------------------

init_eichler_symbol(E) =
{
  my(mff = mffromell(E));
  my(mf = mff[1], f = mff[2]);
  [mf, f, mfsymbol(mf, f)]
}

gamma_from_symbol(fs, tau) =
{
  numeric_clean(2*Pi*I * mfsymboleval(fs, [oo, numeric_clean(tau)]))
}

nearest_period_data(E, z) =
{
  my(om = ellperiods(E), w1 = om[1], w2 = om[2]);
  my(det = real(w1)*imag(w2) - imag(w1)*real(w2));
  my(a = ( real(z)*imag(w2) - imag(z)*real(w2)) / det);
  my(b = (-real(z)*imag(w1) + imag(z)*real(w1)) / det);
  my(m = round(a), n = round(b));
  my(p = m*w1 + n*w2);
  [m, n, p, numeric_clean(z - p)]
}

point_from_gamma(E, gam, lattice_tol = 1e-15) =
{
  my(pd = nearest_period_data(E, gam));
  my(rem = numeric_clean(pd[4]));
  my(lt = safe_real_tol(lattice_tol, 1e-15));

  if(norml2(rem) < lt^2,
    return([[0], pd]);
  );

  [ellztopoint(E, rem), pd]
}

point_from_tau(E, fs, tau, lattice_tol = 1e-15) =
{
  my(gam = gamma_from_symbol(fs, tau));
  my(lt = safe_real_tol(lattice_tol, 1e-15));
  my(R = point_from_gamma(E, gam, lt));
  [gam, R[1], R[2]]
}

\\ -------------------- 候选尖点值集合 Z --------------------

leading_coeff_in_var(P, V) = polcoef(P, poldegree(P, V), V);

candidate_set_Z(E, F, G, Xvar, Yvar, jvar, root_tol = 1e-18, curve_tol =1e-15) =
{
  my(B = leading_coeff_in_var(F, jvar));
  my(C = leading_coeff_in_var(G, jvar));
  my(xroots = [], yroots = [], Z = []);

  if(type(B) == "t_POL", xroots = unique_roots(polroots(B), root_tol));
  if(type(C) == "t_POL", yroots = unique_roots(polroots(C), root_tol));

  for(ix = 1, #xroots,
    for(iy = 1, #yroots,
      if(point_on_curve(E, xroots[ix], yroots[iy], curve_tol),
        my(hit = 0);
        for(k = 1, #Z,
          if(point_close([xroots[ix], yroots[iy]], Z[k][1..2], curve_tol), hit = 1; break);
        );
        if(!hit, Z = concat(Z, [[xroots[ix], yroots[iy], ix, iy]]));
      );
    );
  );
  [B, C, xroots, yroots, Z]
}

match_point_to_Z(P, Z, tol = 1e-8) =
{
  if(#P == 1, return([0, "O"]));
  for(i = 1, #Z,
    if(point_close(P, Z[i][1..2], tol), return([i, Z[i]]));
  );
  [0, "UNMATCHED"]
}

\\ -------------------- 算法 4：尖点值（修正版） --------------------

cusp_lift_matrix(a, c) =
{
  my(v = reduce_pair(a, c));
  a = v[1]; c = v[2];
  if(c == 0, return([1, 0; 0, 1]));
  my(uvg = gcdext(a, c));
  my(u = uvg[1], vv = uvg[2], g = uvg[3]);
  if(g != 1, error("cusp_lift_matrix: pair not primitive"));
  [a, -vv; c, u]
}

cusp_equiv(N, v1, v2) =
{
  my(a = reduce_pair(v1[1], v1[2])[1], c = reduce_pair(v1[1], v1[2])[2]);
  my(ap = reduce_pair(v2[1], v2[2])[1], cp = reduce_pair(v2[1], v2[2])[2]);

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

solve_linear_congruence(a, b, m) =
{
  my(g = gcd(a, m));
  if(Mod(b, g) != 0, error("solve_linear_congruence: no solution"));
  my(a1 = a/g, b1 = b/g, m1 = m/g);
  if(m1 == 1, return(0));
  lift(Mod(b1, m1) / Mod(a1, m1))
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

period_of_matrix(fs, M, tau0 = 1.0*I) =
{
  gamma_from_symbol(fs, (M[1,1]*tau0 + M[1,2])/(M[2,1]*tau0 + M[2,2]))
  - gamma_from_symbol(fs, tau0)
}

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
  [p, srcs, numeric_clean(sumP/(p + 1 - ap))]
}

cusp_image_approx(N, E, fs, cusp, pbound = 1000, tau0 = 1.0*I, lattice_tol = 1e-15) =
{
  my(CG = cusp_gamma_value(N, E, fs, cusp, pbound, tau0));
  my(gam_raw = numeric_clean(CG[3]));
  my(PR = point_from_gamma(E, gam_raw, lattice_tol));
  [CG[1], CG[2], gam_raw, PR[1], PR[2]]
}

exact_cusp_values(N, E, F, G, Xvar, Yvar, jvar,pbound = 1000, tau0 = 1.0*I, root_tol = 1e-18, curve_tol = 1e-15,match_tol = 1e-8, lattice_tol = 1e-15) =
{
  my(Zdata = candidate_set_Z(E, F, G, Xvar, Yvar, jvar, root_tol,curve_tol));
  my(Z = Zdata[5]);
  my(es = init_eichler_symbol(E));
  my(fs = es[3]);
  my(cusps = mfcusps(N));
  my(R = []);

  for(i = 1, #cusps,
    my(cpair = frac_to_pair(cusps[i]));
    my(A = cusp_image_approx(N, E, fs, cpair, pbound, tau0, lattice_tol));
    my(ex = match_point_to_Z(A[4], Z, match_tol));
    R = concat(R, [[cpair, A[1], A[3], A[4], ex, A[5]]]);
  );
  [Zdata, R]
}

\\ -------------------- P = O 时的非尖点纤维（即有限极点） --------------------

fiber_noncuspidal_O(N, E, F, Xvar, jvar,root_tol = 1e-18, lattice_tol = 1e-15, tau_tol = 1e-12) =
{
  my(K = poldegree(F, Xvar));
  my(AK = polcoef(F, K, Xvar));
  my(jroots = [], reps, es, fs, Inv = []);
  my(pole_tol = max(safe_real_tol(lattice_tol, 1e-15), 1e-12));

  \\ A_K(J) 必须是真正的关于 j 的多项式
  if(type(AK) != "t_POL",
    if(AK != 0,
      print("A_K(j) is a nonzero constant, so there are no noncuspidal poles.");
      return([[], [], []]);
    ,
      error("fiber_noncuspidal_O: leading coefficient A_K is zero")
    )
  );

  jroots = unique_roots(polroots(AK), root_tol);

  \\ 这里直接沿用 poles.gp 的做法
  reps = mfcosets(N);

  es = init_eichler_symbol(E);
  fs = es[3];

  print("#roots of A_K(j) = ", #jroots);
  print("#right cosets = ", #reps);

  for(i = 1, #jroots,
    my(taui = tau_from_j(jroots[i], root_tol));

    for(n = 1, #reps,
      my(tau = numeric_clean(mobius(reps[n], taui)));
      my(gam = numeric_clean(gamma_from_symbol(fs, tau)));
      my(pd  = nearest_period_data(E, gam));
      my(rem = numeric_clean(pd[4]));

      \\ 对 O 的判定：看 gamma(tau) 是否足够接近周期格
      if(norml2(rem) < pole_tol^2,
        if(!record_contains_fiber_point(Inv, tau, jroots[i], reps[n], N, tau_tol, root_tol),
          Inv = concat(Inv, [[tau, [0], gam, pd[3], rem, jroots[i], "POLE", reps[n]]]);
        );
      );
    );
  );

  \\ 与普通 fiber_noncuspidal 保持同样的返回格式
  \\ 这里第二项 Jroots 对 O 情形没有意义，置为空向量
  [jroots, [], Inv]
}

\\ -------------------- 有限点纤维搜索（CM 椭圆点按双陪集去重） --------------------

record_contains_tau(V, tau, tol = 1e-12) =
{
  for(i = 1, #V,
    if(type(V[i][1]) == "t_STR", next());
    if(close_complex(V[i][1], tau, tol), return(1));
  );
  0
}

record_contains_fiber_point(V, tau, jroot, M, N, tau_tol = 1e-12, rooteps = 1e-18) =
{
  for(i = 1, #V,
    if(type(V[i][1]) == "t_STR", next());

    \\ 数值上已经是同一个 tau
    if(close_complex(V[i][1], tau, tau_tol), return(1));

    \\ 对 CM 椭圆点，按 Γ0(N)\SL2(Z)/Stab(tau0) 去重
    if(close_complex(V[i][6], jroot, rooteps),
      if(same_x0_point_via_stabilizer(N, jroot, M, V[i][8], rooteps), return(1));
    );
  );
  0
}

fiber_noncuspidal(N, E, P, F, H, Xvar, jvar, Jvar,  root_tol = 1e-18, J_tol = 1e-8, point_tol = 1e-8,  lattice_tol = 1e-15, tau_tol = 1e-12) =
{
  my(alpha, beta);

  \\ =========================
  \\ 关键修补：
  \\ 若 P = [0]，不要再把 alpha 伪装成 0 去解 F(alpha,j), H(alpha,J)
  \\ 而是直接搜索所有有限极点
  \\ =========================
  if(#P == 1,
    return(fiber_noncuspidal_O(N, E, F, Xvar, jvar,
             root_tol, lattice_tol, tau_tol));
  );

  alpha = P[1];
  beta  = P[2];

  my(Fa = subst(F, Xvar, alpha));
  my(Ha = subst(H, Xvar, alpha));
  my(jroots = unique_roots(polroots(Fa), root_tol));
  my(Jroots = unique_roots(polroots(Ha), root_tol));
  my(reps = gamma0_right_cosets(N));
  my(es = init_eichler_symbol(E));
  my(fs = es[3]);
  my(Inv = []);

  print("#roots of F(alpha,j) = ", #jroots);
  print("#roots of f(alpha,J) = ", #Jroots);
  print("#right cosets = ", #reps);

  for(i = 1, #jroots,
    my(taui = tau_from_j(jroots[i], root_tol));
    for(n = 1, #reps,
      my(tau = numeric_clean(mobius(reps[n], taui)));
      my(Jval = numeric_clean(ellj(N*tau)));
      if(!root_in_set(Jval, Jroots, J_tol), next());

      my(PT = point_from_tau(E, fs, tau, lattice_tol));
      my(gam = PT[1], Q = PT[2], pd = PT[3]);

      if(point_close(Q, P, point_tol),
        if(!record_contains_fiber_point(Inv, tau, jroots[i], reps[n], N, tau_tol, root_tol),
          Inv = concat(Inv, [[tau, Q, gam, pd[3], pd[4], jroots[i], Jval, reps[n]]]);
        );
      );
    );
  );

  [jroots, Jroots, Inv]
}

fiber_cusps_from_exact_values(P, ZR, point_tol = 1e-8) =
{
  my(R = ZR[2], C = []);
  for(i = 1, #R,
    my(cusp = R[i][1], approxP = R[i][4], ex = R[i][5]);
    if(type(ex[2]) == "t_STR",
      if(ex[2] == "O" && #P == 1,
        C = concat(C, [["cusp", cusp, approxP, ex]]);
      ,
        if(point_close(approxP, P, point_tol),
          C = concat(C, [["cusp", cusp, approxP, ex]]);
        );
      );
    ,
      if(point_close(ex[2][1..2], P, point_tol) || point_close(approxP, P, point_tol),
        C = concat(C, [["cusp", cusp, ex[2][1..2], ex]]);
      );
    );
  );
  C
}

fiber_of_phi_with_cusps(N, E, P, F, H, G, Xvar, Yvar, jvar, Jvar,root_tol = 1e-18, J_tol = 1e-8,point_tol = 1e-8, lattice_tol = 1e-15,tau_tol = 1e-12, pbound = 1000, tau0 = 1.0*I,cusp_match_tol = 1e-8) =
{
  my(A = fiber_noncuspidal(N, E, P, F, H, Xvar, jvar, Jvar,root_tol,J_tol, point_tol, lattice_tol, tau_tol));
  my(ZR = exact_cusp_values(N, E, F, G, Xvar, Yvar, jvar,pbound, tau0, root_tol, point_tol, cusp_match_tol, lattice_tol));
  my(C = fiber_cusps_from_exact_values(P, ZR, point_tol));
  [A[1], A[2], A[3], C, ZR]
}

\\ -------------------- 便捷入口 --------------------

run_fxj_fxJ_fiber(N, E, P,root_tol = 1e-18, J_tol = 1e-8,point_tol =1e-8, lattice_tol = 1e-15,tau_tol = 1e-12) =
{
  my(Xsym = 'X, jsym = 'j, Jsym = 'J);
  my(F = fxj(Xsym, jsym));
  my(H = fxJ(Xsym, Jsym));
  fiber_noncuspidal(N, E, P, F, H, Xsym, jsym, Jsym,root_tol, J_tol, point_tol, lattice_tol, tau_tol)
}

run_fxj_fxJ_fyj_fiber(N, E, P,root_tol = 1e-18, J_tol = 1e-8,point_tol = 1e-8, lattice_tol = 1e-15,tau_tol = 1e-12, pbound = 1000,tau0 = 1.0*I, cusp_match_tol = 1e-8) =
{
  my(Xsym = 'X, Ysym = 'Y, jsym = 'j, Jsym = 'J);
  my(F = fxj(Xsym, jsym));
  my(H = fxJ(Xsym, Jsym));
  my(G = fyj(Ysym, jsym));
  fiber_of_phi_with_cusps(N, E, P, F, H, G, Xsym, Ysym, jsym, Jsym,root_tol, J_tol, point_tol, lattice_tol,tau_tol, pbound, tau0, cusp_match_tol)
}

\\ -------------------- 打印 --------------------

print_fiber_noncuspidal(R) =
{
  my(jroots = R[1], Jroots = R[2], Inv = R[3]);
  print("========================================");
  print("Noncuspidal fiber points:");
  print("#candidate j-roots = ", #jroots);
  print("#candidate J-roots = ", #Jroots);
  print("#points found = ", #Inv);
  for(i = 1, #Inv,
    print("----------------------------------------");
    print("#", i);
    print("tau            = ", Inv[i][1]);
    print("phi(tau)        = ", Inv[i][2]);
    print("gamma(tau)      = ", Inv[i][3]);
    print("nearest period  = ", Inv[i][4]);
    print("remainder       = ", Inv[i][5]);
    print("root j          = ", Inv[i][6]);
    print("J = j(N tau)    = ", Inv[i][7]);
    print("matrix M        = ", Inv[i][8]);
  );
}

print_fiber_with_cusps(R) =
{
  my(jroots = R[1], Jroots = R[2], Inv = R[3], C = R[4]);
  print_fiber_noncuspidal([jroots, Jroots, Inv]);
  print("========================================");
  print("Cuspidal fiber points:");
  print("#cusps found = ", #C);
  for(i = 1, #C,
    print("----------------------------------------");
    print("#", i);
    print("cusp rep       = [", C[i][2][1], "/", C[i][2][2], "]");
    print("value          = ", C[i][3]);
    print("match          = ", C[i][4]);
  );
}
