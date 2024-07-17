clear
clc
close all

addpath('.\online-calculators-v3\')
addpath('.\Matlab MCMC ensemble sampler\')

nWalks = 4;       % how many MCMC chains?
export = 0;        % do you want to export the data and plots?
filetag = 'test';  % filetag for export

% choose one of these erosion scenarios: 
% 'step',  'samestep',  'samebackground_step', 'samebackground_samestep'
% 'spike', 'samespike', 'samebackground_spike','samebackground_samespike'
scenario = 'samebackground_spike'; 

%% Test data. Use this to see if inversion can recover input
n = 7;   % number of samples
tdata = make_test_data(scenario,n);
Nlogical = [true(n,2) false(n,1)];  % only 10Be and 14C

%% Priors -----------------------------------------------------------------
T =  [1,10e3];      % time of step change OR spike in yrs [min,max]
E =  [10,5e2];      % range of expected erosion rates in mm/ka  [min,max]
LOSS = [0,200];     % loss of soil in cm [min,max], can be commented if no spike model
CHG  = [0, 50];     % change factor of erosion rate, can be commented if no samestep model

[prior_range,var_names] = make_prior_and_varnames(scenario,T,E,LOSS,CHG,n,tdata.steps);

%% Constants

[consts,Nmu] = make_constants();

%% Production rates

sp = Cronus_v3_spallation(tdata.lat,tdata.lon,tdata.altitude,consts);   % get sample parameters (surface procution, pressure)

%% initial guess

mini  = initialmodel_flatprior(prior_range,nWalks,2);

%% Forward model (model parameters: time, erosion, CHG/Loss)

forward_model = @(m) Nforward_wrapper(m,sp,consts,Nmu,scenario,tdata.steps,Nlogical);

%% generate test data to see if inversion can succesfully identify these data

mtest = [tdata.t'; tdata.e'; tdata.changeVariable'];
testObs = forward_model(mtest);


%% Posterior sampling

walkers = Egholm_MCMC(nWalks,testObs,testObs*0.08,mini,prior_range,forward_model);

%% best model
best_obs_err = inf;
for i = 1:nWalks
    [best_obs_err_walker, best_ind] = min(walkers{i}.restot);
    best_model_walker = walkers{i}.u(:,best_ind);
    if best_obs_err_walker < best_obs_err
        best_model = best_model_walker;
    end
end

%% Chain plots

h1 = chainplot(walkers,var_names,prior_range,mtest);

%% Corner plot

h2 = ecornerplot(walkers,'ks',true,'color',[.3 .3 .3],'name',var_names,'bestmodel',best_model,"truevals",mtest);
