clear
clc

%% =====================================================================
%  Linear-strength vortex panel method, closed-form ("XFOIL-style")
%  influence coefficients.
%  SVD/Tikhonov-regularized solve
%  RK4 streamline 

%% =====================================================================

%% -------------------- airfoil geometry input ---------------------
use_generated_airfoil = true;
naca_code = '0006';     % try '2412', '6412', '0012', ... when use_generated_airfoil = true

N_points = 80;           % 60+ recommended
beta = linspace(0,pi,N_points);
x_share = 0.5*(1 - cos(beta))';
x_share(1) = 0; x_share(end) = 1;

if use_generated_airfoil
    [m_c, p_c, t_c] = naca4_params(naca_code);
    [xu_true, yu_true, xl_true, yl_true] = naca4_offset(x_share, m_c, p_c, t_c);
    xu_true(1) = 0; xl_true(1) = 0; yu_true(1) = 0; yl_true(1) = 0;   % pin LE node
    % dense point cloud + arc-length splines - for the camber-line
    % VISUALIZATION search only - not used for panel nodes.
    xg = 0.5*(1-cos(linspace(0,pi,400)))';
    [xu, yu, xl, yl] = naca4_offset(xg, m_c, p_c, t_c);
    xu(1)=0; xl(1)=0; yu(1)=0; yl(1)=0;
    upper_init = [xu, yu]; lower_init = [xl, yl];
else
    coord_init = readmatrix('NACA_6412.txt');
    [~, LE_idx] = min(coord_init(:,1));
    [~, TE_idx] = max(coord_init(:,1));
    upper_init = coord_init([TE_idx:LE_idx],:);   % TE -> LE
    lower_init = coord_init([LE_idx:end],:);      % LE -> TE
    xu = upper_init(:,1); yu = upper_init(:,2);
    xl = lower_init(:,1); yl = lower_init(:,2);
end

%% -- arc-length parametrization of the surfaces (LE -> TE), 
%     used for the camber-line search (Fix 1) and, for FILE-based input only
%     for placing panel nodes too 
if use_generated_airfoil
    up = [xu, yu];             % already LE -> TE from naca4_offset
    lo = [xl, yl];
else
    up = flipud(upper_init);   % reorder to LE -> TE
    lo = lower_init;           % already LE -> TE
end
s_up = [0; cumsum(hypot(diff(up(:,1)), diff(up(:,2))))]; s_up = s_up/s_up(end);
s_lo = [0; cumsum(hypot(diff(lo(:,1)), diff(lo(:,2))))]; s_lo = s_lo/s_lo(end);
pp_ux = pchip(s_up, up(:,1));  pp_uy = pchip(s_up, up(:,2));
pp_lx = pchip(s_lo, lo(:,1));  pp_ly = pchip(s_lo, lo(:,2));

%% -- true mean camber line (kept for visualization/analysis only) --
xc = x_share;
yc = 0.5*(interp1(up(:,1), up(:,2), xc, 'linear', 'extrap') + ...
          interp1(lo(:,1), lo(:,2), xc, 'linear', 'extrap'));
n_iterate = 6;
for outer_it = 1:n_iterate
    dyc = gradient(yc, xc);
    theta = atan(dyc);
    xu_cam = zeros(N_points,1); yu_cam = zeros(N_points,1);
    xl_cam = zeros(N_points,1); yl_cam = zeros(N_points,1);
    for i = 1:N_points
        d = [cos(theta(i)), sin(theta(i))];   % local tangential direction
        C = [xc(i), yc(i)];
        % Well-conditioned scalar root: find arc-length t such that the
        % offset (surface_point(t) - C) is purely NORMAL to the camber
        % tangent d. No ray shooting through an x-indexed spline.
        g_u = @(t) (ppval(pp_ux,t)-C(1))*d(1) + (ppval(pp_uy,t)-C(2))*d(2);
        g_l = @(t) (ppval(pp_lx,t)-C(1))*d(1) + (ppval(pp_ly,t)-C(2))*d(2);
        t0_u = nearest_s_guess(up(:,1), s_up, xc(i));
        t0_l = nearest_s_guess(lo(:,1), s_lo, xc(i));
        try
            tu = fzero(g_u, t0_u);
        catch
            tu = t0_u;
        end
        try
            tl = fzero(g_l, t0_l);
        catch
            tl = t0_l;
        end
        tu = min(max(tu,0),1); tl = min(max(tl,0),1);
        xu_cam(i) = ppval(pp_ux, tu); yu_cam(i) = ppval(pp_uy, tu);
        xl_cam(i) = ppval(pp_lx, tl); yl_cam(i) = ppval(pp_ly, tl);
    end
    xc = 0.5*(xu_cam + xl_cam);
    yc = 0.5*(yu_cam + yl_cam);
end
xc(end) = 1; yc(end) = 0;

figure()
plot(up(:,1), up(:,2), 'b.', lo(:,1), lo(:,2), 'r.');
hold on
plot(xc, yc, 'ko', 'LineWidth', 1.4)
axis equal
grid on
legend('Orig upper','Orig lower','Camber (arc-length normal-intersection)')
title(sprintf('True camber line, NACA %s (analysis only)', naca_code))
hold off

%% -- panel-method surface nodes --
if ~use_generated_airfoil
    xu_true = zeros(N_points,1); yu_true = zeros(N_points,1);
    xl_true = zeros(N_points,1); yl_true = zeros(N_points,1);
    for i = 1:N_points
        t0_u = nearest_s_guess(up(:,1), s_up, x_share(i));
        t0_l = nearest_s_guess(lo(:,1), s_lo, x_share(i));
        gu = @(t) ppval(pp_ux,t) - x_share(i);
        gl = @(t) ppval(pp_lx,t) - x_share(i);
        try, tu = fzero(gu, t0_u); catch, tu = t0_u; end
        try, tl = fzero(gl, t0_l); catch, tl = t0_l; end
        tu = min(max(tu,0),1); tl = min(max(tl,0),1);
        xu_true(i) = ppval(pp_ux, tu); yu_true(i) = ppval(pp_uy, tu);
        xl_true(i) = ppval(pp_lx, tl); yl_true(i) = ppval(pp_ly, tl);
    end
    xu_true(1) = 0; xl_true(1) = 0; yu_true(1) = 0; yl_true(1) = 0;
end


utrue = flipud([xu_true,yu_true]);
ltrue = [xl_true,yl_true];

%% Linear-strength vortex panel method (closed-form coefficients)
nodes = [utrue; ltrue(2:end,:)];
Nn = size(nodes,1);
Np = Nn - 1;

cl_pts    = zeros(Np,2);
panel_len = zeros(Np,1);
panel_nrm = zeros(Np,2);

for i = 1:Np
    dx = nodes(i+1,1) - nodes(i,1);
    dz = nodes(i+1,2) - nodes(i,2);
    L  = sqrt(dx^2 + dz^2);
    panel_len(i) = L;
    panel_nrm(i,:) = [dz, -dx]/L;
    cl_pts(i,:) = 0.5*(nodes(i,:) + nodes(i+1,:));
end

fprintf('Min/median/max panel length: %.6g / %.6g / %.6g (chord units)\n', ...
        min(panel_len), median(panel_len), max(panel_len));

U_inf = 10;
alpha_deg = 5;
alpha = alpha_deg*pi/180;   % radians, used below for the field/streamline plots

%% build geometry-only matrix + SVD ONCE, then solve for this alpha --
[Usvd, sv, Vsvd, lambda] = build_panel_matrix(nodes, panel_len, panel_nrm, Nn, Np, cl_pts);
[CL, GAMMA] = solve_CL(alpha_deg, U_inf, panel_len, panel_nrm, Nn, Np, ...
                        Usvd, sv, Vsvd, lambda, N_points);

%% velocity field for the quiver overview plot
Nf = 40;
xg1 = -0.6; xg2 = 1.6;
zg1 = -1.1; zg2 = 1.1;
xp = linspace(xg1,xg2,Nf)';
zp = linspace(zg1,zg2,Nf)';
[xp,zp] = meshgrid(xp,zp);
svp = size(xp,1);
vp = [xp(:), zp(:)];
vi = zeros(svp^2,2);

for i = 1:svp^2
    [uu, ww] = panel_induced_velocity(vp(i,1), vp(i,2), nodes, panel_len, GAMMA);
    vi(i,:) = [uu, ww];
end
vi(:,1) = vi(:,1) + U_inf*cos(alpha);
vi(:,2) = vi(:,2) + U_inf*sin(alpha);

inside_pts = inpolygon(vp(:,1), vp(:,2), nodes(:,1), nodes(:,2));

% distance-based stand-off buffer + speed cap, purely for a clean-looking
% quiver plot (velocity legitimately grows large arbitrarily close to a
% vortex sheet -- that's expected potential-flow behavior, not noise,
% but it's not useful/plottable at quiver scale)

buffer = 0.01;
cap    = 5*U_inf;
dist_to_body = zeros(svp^2,1);
for i = 1:svp^2
    dist_to_body(i) = min_dist_to_panels(vp(i,:), nodes);
end
near_body = dist_to_body < buffer;
vi(inside_pts | near_body, :) = NaN;
Vmag = hypot(vi(:,1), vi(:,2));
valid = ~isnan(Vmag) & Vmag > 0;
scale = ones(svp^2,1);
scale(valid) = min(1, cap./Vmag(valid));
vi = vi .* scale;

figure()
quiver(vp(:,1),vp(:,2),vi(:,1),vi(:,2)); hold on
plot(nodes(:,1),nodes(:,2),'-m','LineWidth',1.5)
axis equal
xlabel('x/c'); ylabel('z/c');
title(sprintf('Velocity field, NACA %s, \\alpha=%.0f deg, N=%d panels', naca_code, alpha_deg, Np));
hold off

%% -- streamlines via exact-field RK4 tracer (no background grid) --
n_seed = 25;
sx = xg1*ones(n_seed,1);
sy = linspace(-0.5, 0.5, n_seed)';

ds        = 0.008;
max_steps = 4000;
core_r    = max(3*min(panel_len), 1e-4);

figure()
hold on
fill(nodes(:,1), nodes(:,2), [0.8 0.8 0.8], 'EdgeColor','k','LineWidth',1.2);
for s_i = 1:n_seed
    pts = trace_streamline([sx(s_i), sy(s_i)], nodes, panel_len, GAMMA, ...
                            U_inf, alpha, ds, max_steps, ...
                            [xg1 xg2], [zg1 zg2], core_r);
    plot(pts(:,1), pts(:,2), 'b-', 'LineWidth', 0.8);
end
axis equal
xlim([xg1 xg2]); ylim([zg1 zg2]);
xlabel('x/c'); ylabel('z/c');
title(sprintf('Streamlines (exact-field RK4 tracer), NACA %s, \\alpha=%.0f deg, N=%d panels', naca_code, alpha_deg, Np));
hold off

%% ===================== CL-vs-alpha sweep =====================

plot_alpha_deg = [-10,-5,0,5,10];
clplot_n = numel(plot_alpha_deg);
CL_sweep = zeros(clplot_n,1);
for i = 1:clplot_n
    CL_sweep(i) = solve_CL(plot_alpha_deg(i), U_inf, panel_len, panel_nrm, Nn, Np, ...
                            Usvd, sv, Vsvd, lambda, N_points);
end

dat = readmatrix('NACA_0006_dat.txt');

NACA_CL = dat(:,2);
NACA_AoA = dat(:,1);

figure()
plot(plot_alpha_deg, CL_sweep, '-o', 'LineWidth', 1.5, ...
     'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k')
hold on
plot(NACA_AoA,NACA_CL,'-b')
grid on
xlabel('\alpha (deg)'); ylabel('C_L');
title(sprintf('C_L vs \\alpha, NACA %s, N=%d panels', naca_code, Np));
hold off

%% ===================== local functions ===============================

function [m,p,t] = naca4_params(code)
% Put NACA 4-digit code string, e.g. '2412' -> m=0.02, p=0.4, t=0.12
m = str2double(code(1))/100;
p = str2double(code(2))/10;
t = str2double(code(3:4))/100;
end

function [xu, yu, xl, yl] = naca4_offset(x, m, p, t)
% Closed-form NACA-4 thickness + camber offset surfaces, evaluated
% DIRECTLY at the given camber-line x-values -- no spline, no x-indexed
% lookup, so it is immune to the LE non-monotonicity quirk
x = x(:);
yt = (t/0.2)*(0.2969*sqrt(x) - 0.1260*x - 0.3516*x.^2 + 0.2843*x.^3 - 0.1015*x.^4);
yc = zeros(size(x)); dyc = zeros(size(x));
if p > 0
    m1 = x < p;
    yc(m1)  = (m/p^2)*(2*p*x(m1) - x(m1).^2);
    dyc(m1) = (2*m/p^2)*(p - x(m1));
end
if p < 1
    m2 = x >= p;
    yc(m2)  = (m/(1-p)^2)*((1-2*p) + 2*p*x(m2) - x(m2).^2);
    dyc(m2) = (2*m/(1-p)^2)*(p - x(m2));
end
theta = atan(dyc);
xu = x - yt.*sin(theta);  yu = yc + yt.*cos(theta);
xl = x + yt.*sin(theta);  yl = yc - yt.*cos(theta);
end

function t0 = nearest_s_guess(x_orig, s_orig, x_target)
% Monotonic-SAFE initial guess for the arc-length parameter s
% corresponding to a target camber-line x. Works even if x_orig loops
% backward near the LE (unlike interp1(x_orig, s_orig, x_target), which
% requires x_orig to be monotonic and silently misbehaves otherwise).
[~, idx] = min(abs(x_orig - x_target));
t0 = s_orig(idx);
end

function [xi, eta, ct, st] = panel_local_coords(nodeA, nodeB, field_pt)
% Rotate a field point into the panel-local (xi,eta) frame, where the
% panel runs from (0,0) to (Lp,0) along xi. Shared by build_panel_matrix
% (collocation points) and panel_induced_velocity (arbitrary field
% points), so the rotation math can't drift out of sync between them.
x1 = nodeA(1); z1 = nodeA(2);
x2 = nodeB(1); z2 = nodeB(2);
th = atan2(z2-z1, x2-x1);
ct = cos(th); st = sin(th);
dxp = field_pt(1)-x1; dzp = field_pt(2)-z1;
xi  =  dxp*ct + dzp*st;
eta = -dxp*st + dzp*ct;
end

function [I1,I2,J1,J2] = panel_ij(xi, eta, Lp)
% Exact closed-form linear-strength vortex panel integrals in the
% panel-local (xi,eta) frame. The atan(.) form below (NOT atan2) was
% re-verified against independent numerical quadrature to ~1e-14 --
% do not "fix" this to atan2, that was tried and is WRONG by pi
% whenever eta<0 (field point on the far side of the sheet).
r1sq = xi^2 + eta^2;
r2sq = (xi-Lp)^2 + eta^2;
a1 = atan(xi/eta);
a2 = atan((xi-Lp)/eta);
I1 = ((Lp-xi)*(a1-a2) + (eta/2)*log(r1sq/r2sq)) / Lp;
I2 = (xi*(a1-a2)      - (eta/2)*log(r1sq/r2sq)) / Lp;
J1 = -((Lp-xi)/(2*Lp))*log(r1sq/r2sq) - 1 + (eta/Lp)*(a1-a2);
J2 = -(xi/(2*Lp))*log(r1sq/r2sq)      + 1 - (eta/Lp)*(a1-a2);
end

function [Usvd, sv, Vsvd, lambda] = build_panel_matrix(nodes, panel_len, panel_nrm, Nn, Np, cl_pts)
% geometry-ONLY influence matrix + its SVD. Does not depend on
% alpha in any way, so it only needs to be built/factorized ONCE, no
% matter how many angles of attack are later solved via solve_CL.
A = zeros(Nn);

for i = 1:Np
    for j = 1:Np
        if i == j
            % self-induced normal velocity, closed-form limit (exact):
            A(i,j)   = A(i,j)   - 1/(2*pi);
            A(i,j+1) = A(i,j+1) + 1/(2*pi);
        else
            Lp = panel_len(j);
            [xi, eta, ct, st] = panel_local_coords(nodes(j,:), nodes(j+1,:), cl_pts(i,:));
            [I1,I2,J1,J2] = panel_ij(xi, eta, Lp);

            u1 = (I1/(2*pi))*ct - (J1/(2*pi))*st;   w1 = (I1/(2*pi))*st + (J1/(2*pi))*ct;
            u2 = (I2/(2*pi))*ct - (J2/(2*pi))*st;   w2 = (I2/(2*pi))*st + (J2/(2*pi))*ct;

            A(i,j)   = A(i,j)   + u1*panel_nrm(i,1) + w1*panel_nrm(i,2);
            A(i,j+1) = A(i,j+1) + u2*panel_nrm(i,1) + w2*panel_nrm(i,2);
        end
    end
end

A(Nn,1) = 1;
A(Nn,Nn) = 1;

[Usvd, Ssvd, Vsvd] = svd(A);
sv = diag(Ssvd);
lambda = 1e-4 * sv(1);
end

function [CL, GAMMA] = solve_CL(alpha_deg, U_inf, panel_len, panel_nrm, Nn, Np, ...
                                  Usvd, sv, Vsvd, lambda, N_points)

alpha = alpha_deg*pi/180;

RHS = zeros(Nn,1);
for i = 1:Np
    RHS(i) = -(U_inf*cos(alpha)*panel_nrm(i,1) + U_inf*sin(alpha)*panel_nrm(i,2));
end
RHS(Nn) = 0;

sv_reg = sv ./ (sv.^2 + lambda^2);
GAMMA = Vsvd * (sv_reg .* (Usvd' * RHS));

fprintf('cond(A) = %.3e   max|GAMMA| = %.3f\n', sv(1)/sv(end), max(abs(GAMMA)));

Gamma_total = 0;
for k = 1:Np
    Gamma_total = Gamma_total + 0.5*(GAMMA(k)+GAMMA(k+1))*panel_len(k);
end
CL = 2*Gamma_total/U_inf;
fprintf('N_points=%d  alpha=%.1f deg  CL=%.4f\n', N_points, alpha_deg, CL);
end

function [u, w] = panel_induced_velocity(xq, zq, nodes, panel_len, GAMMA)
% Exact closed-form velocity induced by ALL linear-strength vortex
% panels at a single field point (xq, zq). Freestream NOT included.
Np_ = numel(panel_len);
u = 0; w = 0;
for j = 1:Np_
    Lp = panel_len(j);
    [xi, eta, ct, st] = panel_local_coords(nodes(j,:), nodes(j+1,:), [xq, zq]);
    if abs(eta) < 1e-9 && xi > -1e-6 && xi < Lp+1e-6
        continue  % essentially on the panel itself
    end
    [I1,I2,J1,J2] = panel_ij(xi, eta, Lp);
    ul = (GAMMA(j)*I1 + GAMMA(j+1)*I2)/(2*pi);
    wl = (GAMMA(j)*J1 + GAMMA(j+1)*J2)/(2*pi);
    u = u + ul*ct - wl*st;
    w = w + ul*st + wl*ct;
end
end

function d = min_dist_to_panels(P, nodes)
% Minimum point-to-segment distance from P to any panel (vectorized).
A_ = nodes(1:end-1,:); B_ = nodes(2:end,:);
AB = B_ - A_;
denom = sum(AB.^2,2);
t = sum((P - A_).*AB,2)./denom;
t = min(max(t,0),1);
proj = A_ + t.*AB;
d = min(hypot(P(1)-proj(:,1), P(2)-proj(:,2)));
end

function v = unit_velocity_field(P, nodes, panel_len, GAMMA, U_inf, alpha)
[uu, ww] = panel_induced_velocity(P(1), P(2), nodes, panel_len, GAMMA);
uu = uu + U_inf*cos(alpha);
ww = ww + U_inf*sin(alpha);
mag = hypot(uu, ww);
if mag < 1e-8
    v = [0, 0];
else
    v = [uu, ww]/mag;
end
end

function pts = trace_streamline(P0, nodes, panel_len, GAMMA, U_inf, alpha, ...
                                 ds, max_steps, xlim_, zlim_, core_r)
% Marches a streamline by fixed ARCLENGTH steps using RK4 on the exact
% closed-form panel velocity, normalized to unit speed, so near-singular
% speed magnitudes close to the body never destabilize the integration
% the trace just terminates once it enters the core_r stand-off.
pts = zeros(max_steps,2);
pts(1,:) = P0;
P = P0;
for k = 2:max_steps
    if min_dist_to_panels(P, nodes) < core_r
        pts = pts(1:k-1,:); return;
    end
    if P(1)<xlim_(1) || P(1)>xlim_(2) || P(2)<zlim_(1) || P(2)>zlim_(2)
        pts = pts(1:k-1,:); return;
    end
    k1 = unit_velocity_field(P,            nodes, panel_len, GAMMA, U_inf, alpha);
    k2 = unit_velocity_field(P+0.5*ds*k1,  nodes, panel_len, GAMMA, U_inf, alpha);
    k3 = unit_velocity_field(P+0.5*ds*k2,  nodes, panel_len, GAMMA, U_inf, alpha);
    k4 = unit_velocity_field(P+ds*k3,      nodes, panel_len, GAMMA, U_inf, alpha);
    P = P + ds/6*(k1 + 2*k2 + 2*k3 + k4);
    pts(k,:) = P;
end
end

