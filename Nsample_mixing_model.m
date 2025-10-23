function [Nm, Ck, Nm_hist] = Nsample_mixing_model(P0, L, lambda, rho, Zm, Eseg, dtseg)
% Nsample_mixing_model
% Analytic mixed-layer (thickness Zm) surface concentration, Eq. 8 style.
% UNITS:
%   P0: atoms g^-1 yr^-1  |  L: g cm^-2  |  lambda: yr^-1  |  rho: g cm^-3
%   Zm: cm                |  Eseg: cm yr^-1                |  dtseg: yr

S = numel(Eseg);

% ---------- Zm == 0: EXACTLY match master ----------
if Zm == 0
    Nm = Nsurface_eq8_sumexp(P0, L, lambda, rho, Eseg, dtseg);
    Ck = []; 
    Nm_hist = Nm;
    return
end

% ---------- Zm > 0: mixed-layer ----------
exp_fac = exp(-(rho*Zm) ./ L);  % depth factor at z = Zm (mass depth)

Nm_hist = zeros(S,1);

% If the first segment is prehistory (Inf), initialize at steady state for E1 and skip s=1
start_idx = 1;
if isinf(dtseg(1))
    E1    = Eseg(1);
    beta1 = lambda + (rho*E1) ./ L;                 % yr^-1
    Ck    = P0 ./ beta1;                            % substrate coefficients at steady state
    % mixed-layer steady state for segment 1
    Pbar1 = meanP_box(P0, L, rho, Zm);
    Nb_ss = sum( exp_fac .* (P0 ./ beta1) );
    r1    = lambda + E1/Zm;
    Nm    = (Pbar1 + (E1/Zm)*Nb_ss) / r1;
    start_idx = 2;                                  % proceed with next segments
else
    % If no Inf prehistory is provided, start from zeros
    Ck = zeros(size(P0));
    Nm = 0;
end

for s = start_idx:S
    E  = Eseg(s);
    dt = dtseg(s);

    beta = lambda + (rho*E) ./ L;                   % [1xK] yr^-1
    Pbar = meanP_box(P0, L, rho, Zm);               % atoms g^-1 yr^-1

    % N_b(t) = B0 + sum_k Bk e^{-beta_k t}, at z = Zm
    B0 = sum( exp_fac .* (P0 ./ beta) );
    Bk =        exp_fac .* (Ck - (P0 ./ beta));

    r  = lambda + E/Zm;
    S0 = Pbar + (E/Zm)*B0;
    Ak = (E/Zm) .* Bk;

    Nm = Nm * exp(-r*dt) + (S0/r) * (1 - exp(-r*dt)) ...
       + sum( Ak .* phi_term(r, beta, dt) );

    Nm_hist(s) = Nm;

    % substrate coefficients update (Eq. 8 recursion)
    Ck = Ck .* exp(-beta*dt) + (P0 ./ beta) .* (1 - exp(-beta*dt));
end
end

% --------- helpers ---------
function Pbar = meanP_box(P0, L, rho, Zm)
Pbar = sum( P0 .* (L./(rho*Zm)) .* (1 - exp(-(rho*Zm)./L)) );
end

function phi = phi_term(r, beta, dt)
x = (r - beta) * dt;
phi = zeros(size(beta));
tol = 1e-12;
mask = abs(x) < tol;
phi(~mask) = (1 - exp(-x(~mask))) ./ (r - beta(~mask));
phi(mask)  = dt;
end

function Ns = Nsurface_eq8_sumexp(P0, L, lambda, rho, Eseg, dtseg)
% Eq. 8 surface update (no mixing) with advected-depth cutoff and prehistory = dtseg(1)=Inf.

S = numel(Eseg);

% cumulative exhumation BELOW present surface after each segment s:
% Zcum(s) = sum_{q=s+1..S} E(q)*dt(q)
Zcum = zeros(1,S);
if S > 1
    Zcum(1:S-1) = fliplr( cumsum( fliplr(Eseg(2:end).*dtseg(2:end)) ) );
end

Csum = zeros(size(P0)); % per-term accumulator (atoms g^-1)

for s = 1:S
    b   = lambda + (rho*Eseg(s))./L;            % [1xK] yr^-1
    fac = exp( -(rho*Zcum(s)) ./ L );           % advected-depth cutoff
    add = (P0./b) .* fac .* (1 - exp(-b*dtseg(s)));
    Csum = Csum .* exp(-lambda*dtseg(s)) + add; % inter-segment decay
end
Ns = sum(Csum);
end
