clear
clc
%% =====================================================================
%  Vortex Lattice Method (specifically Weissinger)
%  Cosine spaced camber line
%  Geometry visual
%  AIC build
%  CL and CM calculation

%% =====================================================================
%% geometry local
coord_init = readmatrix('NACA_2412.txt');
% find LE TE
[~, LE_idx] = min(coord_init(:,1));
[~, TE_idx] = max(coord_init(:,1));
% split U L
upper_init = coord_init([TE_idx:LE_idx],:);
lower_init = coord_init([LE_idx:end],:);
% start cos space
N_points = 20;
beta = linspace(0,pi,N_points);
x_share = 0.5*(1 - cos(beta))';

% find real cos space (normal from camber line)
xu = upper_init(:,1);
xl = lower_init(:,1);
yu = upper_init(:,2);
yl = lower_init(:,2);
% smooth U and L
pp_u = pchip(xu, yu);
pp_l = pchip(xl, yl);
% initial guess camber
xc = x_share;
xc(1) = 0;
xc(end) = 1;
yc = 0.5*(ppval(pp_u, xc) + ppval(pp_l, xc));
n_iterate = 6;

for outer_it = 1:n_iterate      
    % find grad to find norm later
    dyc = gradient(yc, xc);
    theta = atan(dyc);
    % init real U and L
    xu_true = zeros(N_points,1); 
    yu_true = zeros(N_points,1);
    xl_true = zeros(N_points,1); 
    yl_true = zeros(N_points,1);
    for i = 1:N_points
        st = sin(theta(i));
        ct = cos(theta(i));
        t_guess = 0.5*(ppval(pp_u, xc(i)) - ppval(pp_l, xc(i)));
        % para f for U and L
        f_upper = @(s) ppval(pp_u, xc(i) - s*st) - (yc(i) + s*ct);
        f_lower = @(s) ppval(pp_l, xc(i) - s*st) - (yc(i) + s*ct);
        try
            su = fzero(f_upper, t_guess);
        catch
            su = t_guess; % fallback near LE/TE if fzero struggles
        end
        try
            sl = fzero(f_lower, -t_guess);
        catch
            sl = -t_guess;
        end

        xu_true(i) = xc(i) - su*st;  yu_true(i) = yc(i) + su*ct;
        xl_true(i) = xc(i) - sl*st;  yl_true(i) = yc(i) + sl*ct;
    end

    % camber midpoint of these normal intercepts
    xc = 0.5*(xu_true + xl_true);
    yc = 0.5*(yu_true + yl_true);
end
xc(end) = 1;
xl_true(end) = 1;
yc(end) = 0;
yl_true(end) = 0;
xu_true(end) = 1;
yu_true(end) = 0;

utrue = flipud([xu_true,yu_true]);
ltrue = [xl_true,yl_true];
coord_true = [utrue([1:end-1],:);ltrue([1:end],:)]; % remove double count of LE TE

figure()
plot(xu, yu, 'b.', xl, yl, 'r.'); 
hold on
plot(xc, yc, 'ko', 'LineWidth', 1.4)
plot(xu_true, yu_true, 'go', xl_true, yl_true, 'mo')
axis equal; 
grid on
legend('Orig upper','Orig lower','Camber (normal)','Surface pts (upper)','Surface pts (lower)')
hold off

%% 3d plot
Nhspan = 15;
hspan = 5;
Rchord = 1;
AR = 1;
Tchord = Rchord/AR;
Rtwist = deg2rad(10);
Ttwist = deg2rad(10); 
sweep = deg2rad(0);
dihed = deg2rad(0);
Vmod = 10;
AoA_deg = 10;
AoA = deg2rad(AoA_deg);

Vinf = Vmod*[cos(AoA); -sin(AoA); 0];

rho = 1.225;

Rot_centre = 0.25;

cf = @(z) Rchord + (Tchord - Rchord)*z/hspan;
tf = @(z) Rtwist + (Ttwist - Rtwist)*z/hspan;

camber = @(x) spline(xc,yc,x);

pp_camber  = pchip(xc, yc);
pp_camber_d = fnder(pp_camber, 1);
camber_slope_34 = ppval(pp_camber_d, 0.75);  % slope at the 3/4-chord control point
camber_angle    = atan(camber_slope_34);      % same for every station (2D section is just
                                               % scaled/rotated, not re-shaped, span-wise)

Whspan = hspan/(Nhspan - 1); % width of sections
zspan = 0;
zero_col = zeros(height(coord_true),1);
coord_tran = coord_true;
coord_tran(:,1) = coord_tran(:,1) - Rot_centre; % centre of rotation
coord_rot = [coord_tran,zero_col];
Nhalf = height(utrue); % split upper lower for graphic
Vp = [0.25 - Rot_centre,camber(0.25),0];
Cp = [0.75 - Rot_centre,camber(0.75),0];
Vpi = zeros(Nhspan,3);
chords = zeros(Nhspan,1);
Cpi = zeros(Nhspan - 1,3);
Mpi = zeros(Nhspan - 1,3);
Nor = zeros(Nhspan - 1,3);

dS = zeros(Nhspan-1,1);

figure()
for i = 1:Nhspan
    chord = cf(zspan);
    chords(i) = chord;
    ch_next = cf(zspan+ Whspan);
    twist = tf(zspan);

    rotation = [cos(twist),-sin(twist);sin(twist),cos(twist)];
    coordsnew = coord_tran*rotation*chord;
    % vortex bounds
    Vpi(i,[1,2]) = Vp(1:2)*rotation*chord;
    Vpi(i,1) = Vpi(i,1) + zspan*(tan(sweep));
    Vpi(i,3) = zspan;
    Vpi(i,2) = Vpi(i,2) + zspan*tan(dihed);

    if i < Nhspan
        Cpz = zspan + 0.5*Whspan;
        tw = tf(Cpz);
        ch = cf(Cpz);
        
        rot = [cos(tw),-sin(tw);sin(tw),cos(tw)];
        % control points
        Cpi(i,[1,2]) = Cp(1:2)*rot*ch;
        Cpi(i,1) = Cpi(i,1) + Cpz*tan(sweep);
        Cpi(i,3) = Cpz;
        Cpi(i,2) = Cpi(i,2) + Cpz*tan(dihed);

        Mpi(i,[1,2]) = Vp(1:2)*rot*ch;
        Mpi(i,1) = Mpi(i,1) + Cpz*tan(sweep);
        Mpi(i,3) = Cpz;
        Mpi(i,2) = Mpi(i,2) + Cpz*tan(dihed);

        tw_eff = tw + camber_angle;

        Nor(i,:) = [sin(tw_eff), cos(tw_eff)*cos(dihed), -cos(tw_eff)*sin(dihed)];

        dS(i) = 0.5*(chord + ch_next)*Whspan;
    else
    end
    
    coord_rot(:,1:2) = coordsnew;
    coord_rot(:,3) = zero_col;

    coord_rot(:,1) = coord_rot(:,1) + zspan*(tan(sweep));
    coord_rot(:,2) = coord_rot(:,2) + zspan*tan(dihed);

    if zspan == 0
    else
        reflect = [1,0,0;0,1,0;0,0,-1];
        coord_flip = coord_rot*reflect;
        plot3(coord_flip(1:Nhalf,3),coord_flip(1:Nhalf,1),coord_flip(1:Nhalf,2),'.b');
        plot3(coord_flip(Nhalf + 1:end,3),coord_flip(Nhalf + 1:end,1),coord_flip(Nhalf + 1:end,2),'.r');
    end

    plot3(coord_rot(1:Nhalf,3),coord_rot(1:Nhalf,1),coord_rot(1:Nhalf,2),'.b'); hold on
    plot3(coord_rot(Nhalf + 1:end,3),coord_rot(Nhalf + 1:end,1),coord_rot(Nhalf + 1:end,2),'.r');

    zero_col = zero_col + Whspan;
    zspan = zspan + Whspan;
end
plot3(Cpi(:,3),Cpi(:,1),Cpi(:,2),'xr'); hold on
plot3(Vpi(:,3),Vpi(:,1),Vpi(:,2),'-b')
ylabel('y')
xlabel('x')
zlabel('z')
grid on; axis equal
hold off

%% construct AIC
ex = [1, 0, 0]; 

AIC = zeros(Nhspan - 1);
RHS = zeros(Nhspan - 1,1);

dy = diff(Vpi(:,3));

horseshoe = @(Cpc, Left, Right) local_horseshoe_velocity(Cpc, Left, Right, ex);

for i = 1:Nhspan - 1
    Cpc = Cpi(i,:);
    NorC = Nor(i,:);

    for j = 1:Nhspan - 1
        Left  = Vpi(j,:);
        Right = Vpi(j+1,:);

        Vij_real = horseshoe(Cpc, Left, Right);

        LeftImage  = [Right(1), Right(2), -Right(3)];
        RightImage = [Left(1),  Left(2),  -Left(3)];
        Vij_image = horseshoe(Cpc, LeftImage, RightImage);

        Vij = Vij_real + Vij_image;

        AIC(i,j) = dot(NorC, Vij);
    end
    RHS(i) = -dot(Vinf, NorC);
end

GAMMA = AIC \ RHS;

dL = rho*Vmod.*GAMMA.*dy;

% total wing (both semi-spans, by symmetry) lift and reference area
L_total = 2*sum(dL);
S_total = 2*sum(dS);
CL_approx = L_total/(0.5*rho*Vmod^2*S_total);

%% pitching moment (Cm)
ref_frac = 0.25;                        % 0 = LE, 0.25 = quarter chord, etc.
x_ref    = (ref_frac - Rot_centre) * Rchord;

q     = 0.5*rho*Vmod^2;
C_ref = mean(chords);                   

arm_x = Mpi(:,1) - x_ref;               % Nx1, chordwise moment arm per panel

M_total = -2*sum(arm_x .* dL);

Cm = M_total/(q*S_total*C_ref);

fprintf('Total lift  = %.4f N\n', L_total);
fprintf('Approx CL   = %.4f\n', CL_approx);
fprintf('Total Reference Area = %.4f\n', S_total);
fprintf('Moment ref point (x, frac of root chord) = %.4f (%.4f)\n', x_ref, ref_frac);
fprintf('Approx Cm   = %.4f\n', Cm);

figure()
plot(Cpi(:,3), GAMMA, '-r')
hold on
plot(Cpi(:,3), dL, '-b')
xlabel('Spanwise position z (m)')
ylabel('\Gamma (m^2/s), dL (N/m)')
title('Spanwise circulation distribution (modelled semi-span)')
grid on

function V = local_horseshoe_velocity(Cp, Left, Right, ex)
    r1 = Cp - Left;
    r2 = Cp - Right;
    r0 = r1 - r2;

    cr1r2 = cross(r1, r2);
    bound = cr1r2 / dot(cr1r2, cr1r2) * dot(r0, r1/norm(r1) - r2/norm(r2));

    leg1 = cross(r1, ex) / (norm(r1) * (norm(r1) - dot(r1, ex)));
    leg2 = cross(r2, ex) / (norm(r2) * (norm(r2) - dot(r2, ex)));

    V = 1/(4*pi) * (bound + leg1 - leg2);
end
