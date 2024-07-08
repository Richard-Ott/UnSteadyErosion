clear
clc
close all

addpath('.\online-calculators-v3\')
addpath('.\Matlab MCMC ensemble sampler\')



% choose one of these erosion scenarios: 
% 'step',  'samestep',  'samebackground_step', 'samebackground_samestep'
% 'spike', 'samespike', 'samebackground_spike','samebackground_samespike'

for i = 5
nWalks = 30;       % how many MCMC chains?
export = 1;        % do you want to export the data and plots?
filetag = 'test';  % filetag for export
scenarios = {'step', 'samestep', 'samebackground_step', 'samebackground_samestep',...
    'spike', 'samespike', 'samebackground_spike', 'samebackground_samespike'}; 
%% Test data. Use this to see if inversion can recover input
n = 7;   % number of samples
tdata = make_test_data(scenarios{i},n);

%% Priors -----------------------------------------------------------------
T =  [1,10e3];      % time of step change OR spike in yrs [min,max]
if i == 1; E = [10,5e3]; else; E = [10,5e2]; end    % range of expected erosion rates in mm/ka  [min,max]
LOSS = [0,200];     % loss of soil in cm [min,max], can be commented if no spike model
CHG  = [0, 50];     % change factor of erosion rate, can be commented if no samestep model

[prior_range,var_names] = make_prior_and_varnames(scenarios{i},T,E,LOSS,CHG,n,tdata.steps);

%% Constants

[consts,Nmu] = make_constants();

%% Production rates

sp = Cronus_v3_spallation(tdata.lat,tdata.lon,tdata.altitude,consts);   % get sample parameters (surface procution, pressure)

%% initial guess

mini  = initialmodel_flatprior(prior_range,nWalks,2);

%% Forward model (model parameters: time, erosion, CHG/Loss)

forward_model = @(m) Nforward_wrapper(m,sp,consts,Nmu,scenarios{i},tdata.steps);

%% generate test data to see if inversion can succesfully identify these data

mtest = [tdata.t'; tdata.e'; tdata.changeVariable'];
testObs = forward_model(mtest);

%% log likelihood functions
% First we define a helper function equivalent to calling log(normpdf(x,mu,sigma))
% but has higher precision because it avoids truncation errors associated with calling
% log(exp(xxx)).
lognormpdf = @(x,mu,sigma)-0.5*((x-mu)./sigma).^2  -log(sqrt(2*pi).*sigma);
logLike    = @(m) sum(lognormpdf(testObs, forward_model(m), testObs*0.08));

logical_prior = @(m) sum(and(m > prior_range(:,1), m < prior_range(:,2))) == size(prior_range,1);

%% Posterior sampling
tic
[models, logLike] = gwmcmc(mini,{logical_prior logLike},1e8,'ThinChain',5,'burnin',.2,'ProgressBar',false);
toc
models = single(models); logLike = single(logLike); % save some memory

%% Best-fit model

posterior_like = squeeze(logLike(2,:,:));
[best_walker_like, best_walker_index] = max(posterior_like,[],2);
[best_model_like, best_index] = max(best_walker_like);
best_model = models(:,best_index,best_walker_index(best_index));
best_pred = forward_model(best_model);

%% Autocorrelation

h1 = autocorrelationplot(models);

%% Chain plots

h2 = chainplot(models,var_names,prior_range,mtest);

%% Corner plot of parameters

h3 = ecornerplot(models,'ks',true,'color',[.3 .3 .3],'name',var_names,'bestmodel',best_model,"truevals",mtest);

%% Barplot of parameters

h4 = barplot_parameters(models,var_names,'bestmodel',best_model,'truevals',mtest);

%% Comparison best model and data

h5 = conc_modelledVSobserved(best_pred,testObs(1:n),testObs(1:n).*0.08,testObs(n+1:end),testObs(n+1:end)*0.08);

%% Export

if export
    exportgraphics(h1,['./output/' filetag '_' scenarios{i} '_autocorrelation.png'],'Resolution',300)
    exportgraphics(h2,['./output/' filetag '_' scenarios{i} '_chains.png'],'Resolution',300)
    exportgraphics(h3,['./output/' filetag '_' scenarios{i} '_cornerplot.png'],'Resolution',300)
    exportgraphics(h4,['./output/' filetag '_' scenarios{i} '_barplot.png'],'Resolution',300)
    exportgraphics(h5,['./output/' filetag '_' scenarios{i} '_datafit.png'],'Resolution',300)
end
i
clear
close all
end
