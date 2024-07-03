function N = Nforward_wrapper(model,sp,consts,Nmu,scenario,nsteps)
% This is a wrapper function for Nforward_discretized for easy use with the
% MCMC algorithm.
% Richard Ott, 2024
%
% Input:
%       - model: inversion model as column vector
%       - sp: sample specific parameters (Pspal & pressure)
%       - consts: the consts_v3.mat file from Cronus v3
%       - precalculated muon production rates from Cronus v3 (Nmu.mat)
%       - type: type of erosion scenario. Either: 
%               - 'step': step changes in erosion rates with individual
%               rates per sample
%               - 'samestep': step changes in erosion with common change
%               factor
%               - 'spike': spikes of soil loss
%       - nsteps:  the number of step changes of spike erosion losses
%
% Output:
%       - N = [N10,N14]: concentrations at/g of 10Be and 14C
%
% Richard Ott 2024

nSamp = length(sp.P10spal);

T = [model(1:nsteps); 0];
if ~strcmp(scenario,'step')
    E = model(nsteps+1:end-nsteps);
    changevar = model(end-nsteps+1:(end));
end
%% reshape input
switch scenario
    case 'step'
        E = model(nsteps+1:end);
        E = reshape(E,[nSamp, length(T)]);

    case 'samestep'
        E = model(nsteps+1:end-nsteps);
        changevar = model(end-nsteps+1:(end));

    case 'samebackground_step'
        E = model(nsteps+1);
        changevar = model(nsteps+2:end);
        changevar = reshape(changevar,[nSamp, length(T)-1]);

    case 'samebackground_samestep'
        E = model(nsteps+1);
        changevar = model(nsteps+2:end);

    case 'spike'
        E = model(nsteps+1: nsteps+nSamp);
        changevar = model(nsteps+nSamp+1:end);
        changevar = reshape(changevar,[nSamp, length(T)-1]);

    case 'samespike'
        E = model(nsteps+1: nsteps+nSamp);
        changevar = model(nsteps+nSamp+1:end);

    case 'samebackground_spike'
        E = model(nsteps+1);
        changevar = model(nsteps+2:end);
        changevar = reshape(changevar,[nSamp, length(T)-1]);

    case 'samebackground_samespike'
        E = model(nsteps+1);
        changevar = model(nsteps+2:end);

end

%% run forward model

if strcmp(scenario,'step')
    N = Nforward_discretized(E,T,sp,consts,Nmu,scenario);
else
    N = Nforward_discretized(E,T,sp,consts,Nmu,scenario,changevar);
end

end
