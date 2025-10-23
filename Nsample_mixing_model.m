
function [Nm, Ck, Nm_hist] = Nsample_mixing_model(P0, L, lambda, rho, Zm, Eseg, dtseg, E_init)
% Nmixedbox_eq8
% Analytic mixed-layer (thickness Zm) surface concentration.
% K - production pathways, S -  number of time steps
%
% Inputs
%   P0      [1xK]  surface production (atoms g^-1 yr^-1) for each exponential term
%   L       [1xK]  attenuation lengths (g cm^-2) 
%   lambda  [1x1]  decay constant (yr^-1)
%   rho     [1x1]  density (g cm^-3)
%   Zm      [1x1]  mixed layer thickness (cm) 
%   Eseg    [1xS]  erosion rate per segment (cm/yr)
%   dtseg   [1xS]  duration per segment (yr)
%   E_init  [1x1]  erosion rate before first segment (cm/yr) (for steady-state prehistory)
%
% Outputs
%   Nm      [1x1]  mixed-layer concentration at the end of the last segment
%   Ck      [1xK]  substrate coefficients at segment end
%   Nm_hist [Sx1]  Nm after each segment
%
% Notes
% - Depth is handled in mass depth via exp(-rho*z/L). z in meters.
%
% Richard Ott, 2025

K = numel(P0);      % number of erosion segments
S = numel(Eseg);    % number of time steps

% Handle Zm == 0 as the classic Eq. 8 surface (fast path)
if Zm == 0
    Nm = Nsurface_eq8_sumexp(P0, L, lambda, rho, Eseg, dtseg, E_init);
    Ck = []; Nm_hist = Nm; % for interface consistency
    return
end

% --- precompute depth factors for z = Zm
exp_fac = exp(-(rho*Zm) ./ L);      

% --- initialize substrate coefficients Ck at the start (steady state at E_init)
beta0 = lambda + (rho*E_init) ./ L; 
Ck = P0 ./ beta0;                   % N(z,0) = sum_k Ck * e^{-rho z/L_k}

% --- initial Nm
% mixed-layer steady state for the *first* segment
E1 = Eseg(1);
r1 = lambda + E1/Zm;
Pbar1 = meanP_box(P0, L, rho, Zm);
% N_b steady for first segment at z=Zm using same E1:
beta1 = lambda + (rho*E1) ./ L;
Nb_ss = sum( exp_fac .* (P0 ./ beta1) );
Nm = (Pbar1 + (E1/Zm)*Nb_ss) / r1;


Nm_hist = zeros(S,1);

% --- segment loop
for s = 1:S
    E  = Eseg(s);
    dt = dtseg(s);

    % per-term betas in this segment
    beta = lambda + (rho*E) ./ L; 

    % box-mean production over 0..Zm
    Pbar = meanP_box(P0, L, rho, Zm);

    % N_b(t) at z = Zm is sum of constant + exponentials:
    % B0 = sum exp_fac .* (P0./beta)
    % Bk = exp_fac .* (Ck - P0./beta)   with time factor e^{-beta t}
    B0 = sum( exp_fac .* (P0 ./ beta) );
    Bk =        exp_fac .* (Ck - (P0 ./ beta)); % [1xK]

    % ODE: dNm/dt = Pbar - lambda*Nm - (E/Zm)*(Nm - Nb(t))
    r  = lambda + E/Zm;
    S0 = Pbar + (E/Zm)*B0;
    Ak = (E/Zm) .* Bk; % per-term amplitudes for exponentials

    % exact closed-form update over dt:
    Nm = Nm * exp(-r*dt) + (S0/r) * (1 - exp(-r*dt)) ...
       + sum( Ak .* phi_term(r, beta, dt) );

    Nm_hist(s) = Nm;

    % update substrate coefficients Ck for the next segment (Eq. 8 recursion):
    % Ck_next = Ck*e^{-beta*dt} + (P0./beta) * (1 - e^{-beta*dt})
    Ck = Ck .* exp(-beta*dt) + (P0 ./ beta) .* (1 - exp(-beta*dt));
end
end

% ------------ helpers ---------------

function Pbar = meanP_box(P0, L, rho, Zm)
% Mean production over 0..Zm for sum of exponentials in mass depth
% P(z) = sum_k P0_k * exp(-(rho*z)/L_k)
Pbar = sum( P0 .* (L./(rho*Zm)) .* (1 - exp(-(rho*Zm)./L)) );
end

function phi = phi_term(r, beta, dt)
% phi_k = (1 - exp(-(r - beta_k)*dt)) / (r - beta_k), with safe r≈beta
x = (r - beta) * dt; % vectorized
phi = zeros(size(beta));
% safe evaluation
tol = 1e-10;
mask = abs(x) < tol;
phi(~mask) = (1 - exp(-x(~mask))) ./ (r - beta(~mask));
phi( mask) = dt; % limit as (r-beta)->0
end

function Ns = Nsurface_eq8_sumexp(P0, L, lambda, rho, Eseg, dtseg, E_init)
% Classic Eq. 8 surface update for a sum of exponentials (no mixing), for comparison.
K = numel(P0); S = numel(Eseg);
beta = @(E) lambda + (rho*E)./L;
Ck = P0 ./ beta(E_init); % steady-state prehistory
for s = 1:S
    E  = Eseg(s); dt = dtseg(s);
    b  = beta(E);
    % surface: exp(-rho*0/L)=1
    Ns_term = Ck .* exp(-lambda*dt) + (P0./b) .* (1 - exp(-b*dt));
    Ck = Ns_term;  % coefficients at surface after this segment
end
Ns = sum(Ck);
end
