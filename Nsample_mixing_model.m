function [Nm, Ck, Nm_hist] = Nsample_mixing_model(P0, L, lambda, rho, Zm, Eseg, dtseg, E_init)
% Nsample_mixing_model
% Analytic mixed-layer (thickness Zm) surface concentration, Eq. 8 style.
% K - number of exponential production terms; S - number of time segments
%
% UNITS:
%   P0     [1xK]  atoms g^-1 yr^-1   (surface amplitudes)
%   L      [1xK]  g cm^-2            (attenuation lengths, mass-depth)
%   lambda [1x1]  yr^-1              (decay)
%   rho    [1x1]  g cm^-3            (density)
%   Zm     [1x1]  cm                 (mixed-layer thickness; set 0 for no mixing)
%   Eseg   [1xS]  cm yr^-1           (erosion per segment)
%   dtseg  [1xS]  yr                 (duration per segment)
%   E_init [1x1]  cm yr^-1           (prehistory erosion; not used in Zm==0 fast path)
%
% Richard Ott, 2025

K = numel(P0);
S = numel(Eseg);

% ---- Fast path: Zm == 0 -> EXACTLY match master branch (advected cutoff) ----
if Zm == 0
    Nm = Nsurface_eq8_sumexp(P0, L, lambda, rho, Eseg, dtseg);
    Ck = []; 
    Nm_hist = Nm; % for interface consistency
    return
end

% ---- Mixed-layer path (Zm > 0) ----

% depth factor at z = Zm (mass-depth form)
exp_fac = exp(-(rho*Zm) ./ L);

% initialize substrate coefficients Ck at the start (steady state at E_init)
beta0 = lambda + (rho*E_init) ./ L;
Ck    = P0 ./ beta0;  % N(z,0) = sum_k Ck * exp(-rho z / L_k)

% initial Nm: steady mixed-layer for the first segment
E1   = Eseg(1);
r1   = lambda + E1/Zm;
Pbar = meanP_box(P0, L, rho, Zm);
beta1 = lambda + (rho*E1) ./ L;
Nb_ss = sum( exp_fac .* (P0 ./ beta1) );
Nm    = (Pbar + (E1/Zm)*Nb_ss) / r1;

Nm_hist = zeros(S,1);

% segment loop
for s = 1:S
    E  = Eseg(s);
    dt = dtseg(s);

    % per-term betas (yr^-1) this segment
    beta = lambda + (rho*E) ./ L;

    % box-mean production (atoms g^-1 yr^-1)
    Pbar = meanP_box(P0, L, rho, Zm);

    % N_b(t) at z = Zm = B0 + sum_k Bk exp(-beta_k t)
    B0 = sum( exp_fac .* (P0 ./ beta) );
    Bk =        exp_fac .* (Ck - (P0 ./ beta)); % [1xK]

    % mixed-layer ODE: dNm/dt = Pbar - lambda*Nm - (E/Zm)*(Nm - Nb(t))
    r  = lambda + E/Zm;
    S0 = Pbar + (E/Zm)*B0;
    Ak = (E/Zm) .* Bk;

    % exact update over dt
    Nm = Nm * exp(-r*dt) + (S0/r) * (1 - exp(-r*dt)) ...
       + sum( Ak .* phi_term(r, beta, dt) );

    Nm_hist(s) = Nm;

    % substrate coefficient update (Eq. 8 recursion)
    Ck = Ck .* exp(-beta*dt) + (P0 ./ beta) .* (1 - exp(-beta*dt));
end
end

% ------------ helpers ---------------

function Pbar = meanP_box(P0, L, rho, Zm)
% Mean production over 0..Zm for P(z) = sum_k P0_k * exp(-(rho*z)/L_k), z in cm
Pbar = sum( P0 .* (L./(rho*Zm)) .* (1 - exp(-(rho*Zm)./L)) );
end

function phi = phi_term(r, beta, dt)
% phi_k = (1 - exp(-(r - beta_k)*dt)) / (r - beta_k), safe near r≈beta
x = (r - beta) * dt;
phi = zeros(size(beta));
tol = 1e-12;
mask = abs(x) < tol;
phi(~mask) = (1 - exp(-x(~mask))) ./ (r - beta(~mask));
phi(mask)  = dt;
end

function Ns = Nsurface_eq8_sumexp(P0, L, lambda, rho, Eseg, dtseg)
% Eq. 8 surface update (no mixing) WITH advected-depth cutoff like master.
% UNITS: same as main function (Eseg in cm/yr; L in g/cm^2; rho in g/cm^3)

S = numel(Eseg);

% cumulative exhumation above present surface after each segment
Zcum = zeros(1,S);
if S > 1
    Zcum(1:S-1) = fliplr( cumsum( fliplr(Eseg(2:end).*dtseg(2:end)) ) );
end

% segment-by-segment sum (matches master logic)
Csum = zeros(size(P0)); % per-term accumulator at the surface
for s = 1:S
    b   = lambda + (rho*Eseg(s))./L;      % [1xK]
    fac = exp( -(rho*Zcum(s)) ./ L );     % advected-depth cutoff at segment s
    add = (P0./b) .* fac .* (1 - exp(-b*dtseg(s)));
    Csum = Csum .* exp(-lambda*dtseg(s)) + add;
end
Ns = sum(Csum);
end
