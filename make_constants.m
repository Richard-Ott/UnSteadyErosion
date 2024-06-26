function [consts,Nmu] = make_constants()
% get Cronus v3 constants and add some missing constants. Also, load
% pre-calulated muon production rates

consts = load("consts_v3.mat"); consts = consts.consts;
consts.density = 2.65;   % density g/cm3
consts.L_sp    = 160;    % spallation attenuation length g/cm2
load('..\..\..\CRONUS cosmo calculation\cronus 3.0\online-calculators-v3\data\Nmu.mat');    % load precaulcated muon rates

end