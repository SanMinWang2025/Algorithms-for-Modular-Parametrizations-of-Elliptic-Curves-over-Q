\\ ============================================================
\\ branch_points.gp
\\
\\ 按 Arxiv1.tex 第 6 节“寻找分歧点的算法”实现（PARI/GP）
\\
\\ 依赖：请先 read("fibers.gp");
\\ 说明：fibers.gp 已经包含
\\   - 算法 2：fiber_noncuspidal / fiber_of_phi_with_cusps
\\   - 算法 4：exact_cusp_values / fiber_cusps_from_exact_values
\\
\\ 本文件补上：
\\   (1) 由四个模多项式构造候选多项式 U(x), V(y)；
\\   (2) 枚举 E 上候选点并调用算法 2 验证；
\\   (3) 尖点部分：先求尖点在椭圆曲线上的像 P，
\\       再计算纤维 phi^{-1}(P)；若纤维大小小于模参数化次数 d，
\\       则 P 为分歧点，否则不是。
\\
\\ 主要入口：
\\   BR = branch_points_complete(N, E, F, H, G, K, Xsym, Ysym, jsym, Jsym);
\\
\\ 其中
\\   F = F_N(X,j),
\\   H = f_N(X,J),
\\   G = G_N(Y,j),
\\   K = g_N(Y,J).
\\
\\ 若你的脚本中函数名就是 fxj, fxJ, fyj, fyJ，可直接用：
\\   BR = run_fxj_fxJ_fyj_fyJ_branch(N, E);
\\
\\ 返回对象 BR 的结构：
\\   [UxData, VyData, CandData, ZR, FiniteData, CuspData, BranchValues]
\\ 其中
\\   UxData = [R1, R2, U]
\\   VyData = [R3, R4, V]
\\   CandData = [xroots, yroots, candidate_points]
\\   ZR     = exact_cusp_values(...) 的原始输出
\\   FiniteData 的元素：
\\      [P, fiber_size, moddeg, is_branch, FiberRecord]
\\   CuspData 的元素：
\\      [cusp_pair, image_point, fiber_size, is_branch, note, FiberRecord]
\\   BranchValues = 所有分歧值（去重后）的点向量。
\\ ============================================================

default(realprecision, 80);

\\ -------------------- 多项式工具 --------------------

bp_poly_primpart(P, V) =
{
  if(type(P) != "t_POL", return(P));
  my(c = content(P));
  if(c == 0, return(P));
  P/c
}

bp_poly_monic(P, V) =
{
  if(type(P) != "t_POL", return(P));
  my(Q = bp_poly_primpart(P, V));
  if(poldegree(Q, V) <= 0, return(Q));
  my(lc = polcoef(Q, poldegree(Q, V), V));
  if(lc == 0, return(Q));
  Q/lc
}

bp_poly_squarefree_part(P, V) =
{
  if(type(P) != "t_POL", return(P));
  my(Q = bp_poly_primpart(P, V));
  if(poldegree(Q, V) <= 0, return(Q));
  my(D = deriv(Q, V));
  if(D == 0, return(bp_poly_monic(Q, V)));
  bp_poly_monic(Q / gcd(Q, D), V)
}

\\ 计算公共平方自由因子，但全程不做因式分解。
\\ 理论上：common sqfree factor = sqfree(gcd(P,Q)).
\\ 这通常比先分别做 sqfree(P), sqfree(Q) 再 gcd 更省。
bp_common_squarefree_factor(P, Q, V) =
{
  if(type(P) != "t_POL" || type(Q) != "t_POL", return(1));
  my(G = gcd(P, Q));
  if(type(G) != "t_POL", return(G));
  bp_poly_squarefree_part(G, V)
}

bp_candidate_x_data(F, H, Xvar, jvar, Jvar) =
{
  my(R1 = polresultant(F, deriv(F, jvar), jvar));
  my(R2 = polresultant(H, deriv(H, Jvar), Jvar));
  my(U  = bp_common_squarefree_factor(R1, R2, Xvar));
  [R1, R2, U]
}

bp_candidate_y_data(G, K, Yvar, jvar, Jvar) =
{
  my(R3 = polresultant(G, deriv(G, jvar), jvar));
  my(R4 = polresultant(K, deriv(K, Jvar), Jvar));
  my(V  = bp_common_squarefree_factor(R3, R4, Yvar));
  [R3, R4, V]
}

\\ -------------------- 候选点 --------------------

bp_add_unique_point(S, P, tol = 1e-10) =
{
  for(i = 1, #S,
    if(point_close(S[i], P, tol), return(S));
  );
  concat(S, [P])
}

bp_candidate_points_on_curve(E, U, V, Xvar, Yvar, root_tol = 1e-18, curve_tol = 1e-12) =
{
  my(xroots = [], yroots = [], Cand = []);

  if(type(U) == "t_POL" && poldegree(U, Xvar) >= 1,
    xroots = unique_roots(polroots(U), root_tol);
  );
  if(type(V) == "t_POL" && poldegree(V, Yvar) >= 1,
    yroots = unique_roots(polroots(V), root_tol);
  );

  for(i = 1, #xroots,
    for(j = 1, #yroots,
      if(point_on_curve(E, xroots[i], yroots[j], curve_tol),
        Cand = bp_add_unique_point(Cand, [numeric_clean(xroots[i]), numeric_clean(yroots[j])], curve_tol);
      );
    );
  );
  [xroots, yroots, Cand]
}

\\ -------------------- 纤维大小与分歧判定 --------------------
\\ 注意：
\\ 当 P = [0] 时，fiber_noncuspidal(...) 已在 fibers.gp 中自动切换为
\\ “非尖点极点搜索”模式，不再使用 F(alpha,j), H(alpha,J) 代入法。

bp_fiber_size(FR) = #FR[3] + #FR[4];

bp_is_branch_from_fiber(FR, d) = bp_fiber_size(FR) < d;

bp_build_fiber_record_from_ZR(N, E, P, F, H, G, Xvar, Yvar, jvar, Jvar, ZR,  root_tol = 1e-18, J_tol = 1e-8, point_tol = 1e-8, lattice_tol = 1e-15, tau_tol = 1e-12) =
{
  my(A = fiber_noncuspidal(N, E, P, F, H, Xvar, jvar, Jvar,                          root_tol, J_tol, point_tol, lattice_tol, tau_tol));
  my(C = fiber_cusps_from_exact_values(P, ZR, point_tol));
  [A[1], A[2], A[3], C, ZR]
}

\\ -------------------- 非尖点分歧点坐标乘积 --------------------

bp_finite_branch_points(BR) =
{
  my(FD = BR[5], Out = []);
  for(i = 1, #FD,
    if(FD[i][4], Out = concat(Out, [FD[i][1]]));
  );
  Out
}

bp_coordinate_product_poly(Pts, idx, V) =
{
  my(S = 1);
  for(i = 1, #Pts,
    if(type(Pts[i]) == "t_VEC" && #Pts[i] >= idx,
      S *= (V - Pts[i][idx]);
    );
  );
  S
}

bp_poly_bestappr(P, V) =
{
  if(type(P) != "t_POL", return(bestappr(P)));
  my(d = poldegree(P, V));
  if(d < 0, return(bestappr(P)));
  my(S = 0);
  for(n = 0, d,
    S += bestappr(polcoef(P, n, V)) * V^n;
  );
  S
}

\\ -------------------- 尖点：由纤维大小判定是否分歧 --------------------

bp_cusp_exact_point(rec) =
{
  my(ex = rec[5]);
  if(type(ex[2]) == "t_STR",
    if(ex[2] == "O", return([0]));
    return(rec[4]);
  );
  ex[2][1..2]
}

bp_cusp_exact_status(rec) =
{
  my(ex = rec[5]);
  if(type(ex[2]) == "t_STR", return(ex[2]));
  "MATCHED_Z"
}

bp_find_cached_point(Cache, P, tol = 1e-8) =
{
  for(i = 1, #Cache,
    if(point_close(Cache[i][1], P, tol), return(i));
  );
  0
}

bp_cusp_branch_note(status, cnt, d) =
{
  Str("FIBER(status=", status, ", size=", cnt, "/", d, ")")
}

bp_cusp_branch_data_via_fibers(N, E, ZR, F, H, G, Xvar, Yvar, jvar, Jvar, d,  root_tol = 1e-18, J_tol = 1e-8, point_tol = 1e-8,  lattice_tol = 1e-15, tau_tol = 1e-12) =
{
  my(R = ZR[2]);
  my(Cache = []);
  my(Out = []);

  for(i = 1, #R,
    my(cusp = R[i][1]);
    my(P0 = bp_cusp_exact_point(R[i]));
    my(status = bp_cusp_exact_status(R[i]));
    my(idx = bp_find_cached_point(Cache, P0, point_tol));
    my(cnt, isb, note, FR);

    if(idx == 0,
      FR = bp_build_fiber_record_from_ZR(N, E, P0, F, H, G,             Xvar, Yvar, jvar, Jvar, ZR,             root_tol, J_tol, point_tol, lattice_tol, tau_tol);
      cnt = bp_fiber_size(FR);
      isb = bp_is_branch_from_fiber(FR, d);
      note = bp_cusp_branch_note(status, cnt, d);
      Cache = concat(Cache, [[P0, cnt, isb, note, FR]]);
      idx = #Cache;
    );

    Out = concat(Out, [[cusp, Cache[idx][1], Cache[idx][2], Cache[idx][3], Cache[idx][4], Cache[idx][5]]]);
  );
  Out
}

\\ -------------------- 主程序 --------------------

branch_points_complete(N, E, F, H, G, K, Xvar, Yvar, jvar, Jvar,  root_tol = 1e-18, J_tol = 1e-8, point_tol = 1e-8,  curve_tol = 1e-12, lattice_tol = 1e-15, tau_tol = 1e-12,  pbound = 1000, tau0 = 1.0*I, cusp_match_tol = 1e-8,  assume_twist_minimal = -1) =
{
  my(UxData = bp_candidate_x_data(F, H, Xvar, jvar, Jvar));
  my(VyData = bp_candidate_y_data(G, K, Yvar, jvar, Jvar));
  my(U = UxData[3], V = VyData[3]);
  my(CandData = bp_candidate_points_on_curve(E, U, V, Xvar, Yvar, root_tol, curve_tol));
  my(Candidates = CandData[3]);
  my(ZR = exact_cusp_values(N, E, F, G, Xvar, Yvar, jvar,                            pbound, tau0, root_tol, curve_tol,                            cusp_match_tol, lattice_tol));

  my(d = ellmoddegree(E));
  my(FiniteData = []);

  print("#candidate x-roots = ", #CandData[1]);
  print("#candidate y-roots = ", #CandData[2]);
  print("#candidate points on E = ", #Candidates);
  print("modular degree = ", d);

  for(i = 1, #Candidates,
    my(P = Candidates[i]);
    print("[finite candidate ", i, "/", #Candidates, "] P = ", P);
    my(FR = bp_build_fiber_record_from_ZR(N, E, P, F, H, G,                Xvar, Yvar, jvar, Jvar, ZR,                root_tol, J_tol, point_tol, lattice_tol, tau_tol));
    my(cnt = bp_fiber_size(FR));
    my(isb = bp_is_branch_from_fiber(FR, d));
    FiniteData = concat(FiniteData, [[P, cnt, d, isb, FR]]);
  );

  my(CuspData = bp_cusp_branch_data_via_fibers(N, E, ZR, F, H, G,                Xvar, Yvar, jvar, Jvar, d,                root_tol, J_tol, point_tol, lattice_tol, tau_tol));

  my(BranchValues = []);
  for(i = 1, #FiniteData,
    if(FiniteData[i][4], BranchValues = bp_add_unique_point(BranchValues, FiniteData[i][1], point_tol));
  );
  for(i = 1, #CuspData,
    if(CuspData[i][4], BranchValues = bp_add_unique_point(BranchValues, CuspData[i][2], point_tol));
  );

  [UxData, VyData, CandData, ZR, FiniteData, CuspData, BranchValues]
}

\\ -------------------- 便捷入口 --------------------

run_fxj_fxJ_fyj_fyJ_branch(N, E,  root_tol = 1e-18, J_tol = 1e-8, point_tol = 1e-8,  curve_tol = 1e-12, lattice_tol = 1e-15, tau_tol = 1e-12,  pbound = 1000, tau0 = 1.0*I, cusp_match_tol = 1e-8,  assume_twist_minimal = -1) =
{
  my(Xsym = 'X, Ysym = 'Y, jsym = 'j, Jsym = 'J);
  my(F = fxj(Xsym, jsym));
  my(H = fxJ(Xsym, Jsym));
  my(G = fyj(Ysym, jsym));
  my(K = fyJ(Ysym, Jsym));
  branch_points_complete(N, E, F, H, G, K, Xsym, Ysym, jsym, Jsym,    root_tol, J_tol, point_tol, curve_tol, lattice_tol, tau_tol,    pbound, tau0, cusp_match_tol, assume_twist_minimal)
}

\\ -------------------- 打印辅助 --------------------

print_branch_polynomials(BR) =
{
  print("================ U(x) side ================");
  \\print("R1(x) = Res_j(F, dF/dj)");
  \\print(BR[1][1]);
  \\print("R2(x) = Res_J(H, dH/dJ)");
  \\print(BR[1][2]);
  print("U(x)  = common squarefree factor");
  print(BR[1][3]);

  print("================ V(y) side ================");
  \\print("R3(y) = Res_j(G, dG/dj)");
  \\print(BR[2][1]);
  \\print("R4(y) = Res_J(K, dK/dJ)");
  \\print(BR[2][2]);
  print("V(y)  = common squarefree factor");
  print(BR[2][3]);
}

print_branch_candidates(BR) =
{
  my(C = BR[3][3]);
  print("#candidate points on E = ", #C);
  for(i = 1, #C,
    print("Cand[", i, "] = ", C[i]);
  );
}

print_finite_branch_data(BR) =
{
  my(FD = BR[5]);
  print("============== finite candidate check ==============");
  for(i = 1, #FD,
    print("----------------------------------------");
    print("#", i, "  P = ", FD[i][1]);
    print("fiber size = ", FD[i][2], " / modular degree = ", FD[i][3]);
    print("is branch? = ", FD[i][4]);
  );
}

print_cusp_branch_data(BR) =
{
  my(CD = BR[6]);
  print("============== cusp check ==============");
  for(i = 1, #CD,
    print("----------------------------------------");
    print("#", i, "  cusp = [", CD[i][1][1], "/", CD[i][1][2], "]");
    print("image      = ", CD[i][2]);
    print("fiber size = ", CD[i][3]);
    print("is branch? = ", CD[i][4]);
    print("note       = ", CD[i][5]);
  );
}

print_branch_coordinate_products(BR) =
{
  my(Pts = bp_finite_branch_points(BR));
  my(x = 'x, y = 'y);
  my(sx = bp_coordinate_product_poly(Pts, 1, x));
  my(sy = bp_coordinate_product_poly(Pts, 2, y));

  print("============== finite branch coordinate products ==============");
  print("#finite non-cuspidal branch points = ", #Pts);
  print("s_x(x) = prod_i (x - x_i)");
  print(sx);
  print("bestappr(s_x) = ");
  print(bp_poly_bestappr(sx, x));
  print("s_y(y) = prod_i (y - y_i)");
  print(sy);
  print("bestappr(s_y) = ");
  print(bp_poly_bestappr(sy, y));
}

print_branch_summary(BR) =
{
  my(BV = BR[7]);
  print("============== branch values ==============");
  print("#branch values = ", #BV);
  for(i = 1, #BV,
    print("B[", i, "] = ", BV[i]);
  );
}

print_branch_report(BR) =
{
  print_branch_polynomials(BR);
  print_branch_candidates(BR);
  print_finite_branch_data(BR);
  print_cusp_branch_data(BR);
  print_branch_coordinate_products(BR);
  print_branch_summary(BR);
}
