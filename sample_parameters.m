function [sp] = sample_parameters(lat,lon,altitude,consts)
% get top spallation production from Cronus v3 for 10Be, 14C, 26Al
% multi sample input possible (provide variables as vectors)
% Richard Ott, 2024

addpath('C:\Users\rott\OneDrive - UvA\Richard\PhD_ETH\matlab\CRONUS cosmo calculation\cronus 3.0\online-calculators-v3')
addpath('C:\Users\rott\OneDrive - UvA\Richard\PhD_ETH\matlab\CRONUS cosmo calculation\cronus 3.0\online-calculators-v3\data')

sp.pressure = ERA40atm(lat, lon, altitude);  % atmospheric pressure

if isrow(sp.pressure)   % make sure its a column vector
    sp.pressure = sp.pressure';
end

P10_SLHL = consts.refP_St(4);                % Stone scaling ref spallation production
P14_SLHL = consts.refP_St(5);
P26_SLHL = consts.refP_St(7);

sf_St = (stone2000(lat,sp.pressure,1))';        % scaling factor

sp.P10spal = sf_St .* P10_SLHL;              % site specific spallation production
sp.P14spal = sf_St .* P14_SLHL;
sp.P26spal = sf_St .* P26_SLHL;

%% define soil mixing parameters (maybe should be moved out of this function)
sp.mix = false; % set soil mixing per default as false. it gets upated later in main script
sp.mixing = 0;
sp.dzmix = 1;  % put a 1 here that no error are being produced, but the zero for sp.mixing will result in no mixing

%% muon production and attenuation approximated by fitting Model 1A

mindepth = 0; % g/cm2
maxdepth = 7800; % g/cm2

% Muons 10Be (values taken from Hippe et al. 2021, NatComm code) ----------
f_star=0.00191; %Model 1A, alpha=1;
Natoms = 2.006e22; %Oxygen atoms pr gram Quartz
sigma0 = 0.280e-30; % model 1A, alpha=1;

for i = 1:length(sp.pressure)
    p_muons=p_rate_calc2(f_star,Natoms,sigma0,sp.pressure(i),mindepth,maxdepth);

    sp.L_fm_10(i) = p_muons.L(1); % L_eff for fast muons
    sp.L_nm_10(i) = p_muons.L(2); % L_eff for muon capture
    sp.P10_fm(i) = p_muons.P(1);   % Production rate at surface
    sp.P10_nm(i) = p_muons.P(2);   % Surface production rate
end


% Muons 14C ---------------------------------------------------------------
f_star=0.137; % Model 1A, alpha=1;
Natoms = 2.006e22; % Oxygen atoms pr gram Quartz
sigma0 = 2.37e-30; % model 1A, alpha=1;

for i = 1:length(sp.pressure)
    p_muons=p_rate_calc2(f_star,Natoms,sigma0,sp.pressure(i),mindepth,maxdepth);

    sp.L_fm_14(i) = p_muons.L(1); % L_eff for fast muons
    sp.L_nm_14(i) = p_muons.L(2); % L_eff for muon capture
    sp.P14_fm(i) = p_muons.P(1);       % surface production fast muons
    sp.P14_nm(i) = p_muons.P(2);       % surface production muon capture
end

%Muons 26Al
f_star=0.0133; %Model 1A, alpha=1;
Natoms = 1.003e22; %Si atoms pr gram Quartz
sigma0 = 3.89e-30; % model 1A, alpha=1;

for i = 1:length(sp.pressure)
    p_muons=p_rate_calc2(f_star,Natoms,sigma0,sp.pressure(i),mindepth,maxdepth);
    
    sp.L_fm_26(i) = p_muons.L(1); 
    sp.L_nm_26(i) = p_muons.L(2); 
    sp.P26_fm(i) = p_muons.P(1); 
    sp.P26_nm(i) = p_muons.P(2); 
end

% constant attenuation lengths
% consts.L_mu_nmc = 1500;   % g/cm²
% consts.L_mu_fm  = 4320;

end