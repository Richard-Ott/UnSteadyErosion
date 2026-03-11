clear
clc
close all


%filePath = matlab.desktop.editor.getActiveFilename;
%ix = strfind(filePath,filesep);
%basepath = filePath(1:ix(end));
%cd(basepath)
%addpath(genpath(basepath))

nWalks = 5;         % number of independent HMC chains
export = 1;         % do you want to export the data and plots?
filetag = 'test';   % filetag for export

useTuning = false;  % numerical-gradient HMC tuning is very expensive for this model

numTune = 20;       % tuneSampler iterations
numBurnin = 5000;    % discarded warmup draws per chain
numSamples = 200000;    % retained draws per chain
thinSize = 2;       % thinning of retained chain

stepSize = 0.01;    % manual HMC step size for quick tests
numLeapfrogSteps = 10;  % manual leapfrog steps for quick tests

% choose one of these erosion scenarios:
% 'step',  'samestep',  'samebackground_step', 'samebackground_samestep'
% 'spike', 'samespike', 'samebackground_spike','samebackground_samespike'
% 'curve'
scenario = 'step';

%% Test data. Use this to see if inversion can recover input
n = 7;   % number of samples
tdata = make_test_data(scenario,n);
Nlogical = [true(n,2) false(n,1)];  % only 10Be and 14C

%% Priors -----------------------------------------------------------------
T =  [1e2,6e3];     % time of step change OR spike in yrs [min,max]
E =  [0,4.5e3];     % expected erosion rates in mm/ka [min,max]
LOSS = [0,200];     % loss of soil in cm [min,max]
CHG  = [0,50];      % change factor of erosion rate

[prior_range,var_names] = make_prior_and_varnames(scenario,T,E,LOSS,CHG,n,tdata.steps);

%% Constants

consts = make_constants();

%% Production rates

sp = sample_parameters(tdata.lat,tdata.lon,tdata.altitude,consts);

%% Initial guesses inside prior

mini = initialmodel_flatprior(prior_range,nWalks,tdata.steps);

%% Forward model (model parameters: time, erosion, CHG/Loss)
if strcmp(scenario,'curve')
    sp.curvechange = tdata.curvechange;
    sp.t = tdata.t;
end

forward_model = @(m) Nforward_wrapper(m,sp,consts,scenario,tdata.steps,Nlogical);

%% Generate test data to see if inversion can identify these data

mtest = [tdata.t'; tdata.e'; tdata.changeVariable'];
testObs = forward_model(mtest);
sigmaObs = testObs * 0.08;

%% Log likelihood in constrained parameter space
lognormpdf = @(x,mu,sigma)-0.5*((x-mu)./sigma).^2 - log(sqrt(2*pi).*sigma);
logLike = @(m) sum(lognormpdf(testObs, forward_model(m), sigmaObs));

%% Unconstrained transform for HMC
lower_bounds = prior_range(:,1);
upper_bounds = prior_range(:,2);

models = nan(size(prior_range,1), nWalks, numSamples);
logLikeStore = nan(2, nWalks, numSamples);
accRatio = nan(nWalks,1);
endPoints = nan(size(prior_range,1), nWalks);
hmcSamplers = cell(nWalks,1);
tuningInfo = cell(nWalks,1);

tic
for wix = 1:nWalks
    wix
    startModel = mini(:,wix);
    startPoint = bounded_to_unbounded(startModel, lower_bounds, upper_bounds);
    logpdf = @(u) transformed_logpdf(u, lower_bounds, upper_bounds, logLike);

    hmc = hmcSampler(logpdf, startPoint, ...
        'UseNumericalGradient', true, ...
        'CheckGradient', false, ...
        'VariableNames', var_names, ...
        'StepSize', stepSize, ...
        'NumSteps', numLeapfrogSteps, ...
        'StepSizeTuningMethod', 'dual-averaging', ...
        'MassVectorTuningMethod', 'hessian');

    if useTuning
        [hmcSamplers{wix}, tuningInfo{wix}] = tuneSampler(hmc, ...
            'NumStepSizeTuningIterations', numTune, ...
            'VerbosityLevel', 0);
    else
        hmcSamplers{wix} = hmc;
        tuningInfo{wix} = struct();
    end

    [chainU, endpointU, accRatio(wix)] = drawSamples(hmcSamplers{wix}, ...
        'Burnin', numBurnin, ...
        'NumSamples', numSamples, ...
        'ThinSize', thinSize, ...
        'VerbosityLevel', 0);

    chainM = bounded_from_unbounded(chainU', lower_bounds, upper_bounds);
    models(:,wix,:) = reshape(chainM, size(chainM,1), 1, size(chainM,2));
    endPoints(:,wix) = bounded_from_unbounded(endpointU, lower_bounds, upper_bounds);

    for six = 1:numSamples
        logLikeStore(1,wix,six) = 0;
        logLikeStore(2,wix,six) = logLike(models(:,wix,six));
    end
end
toc

models = single(models);
logLikeStore = single(logLikeStore);

%% Best-fit model

posterior_like = reshape(logLikeStore(2,:,:), nWalks, numSamples);
[best_model_like, linear_index] = max(posterior_like(:));
[best_index, best_sample_index] = ind2sub(size(posterior_like), linear_index);
best_model = models(:,best_index,best_sample_index);
best_pred = forward_model(best_model);

fprintf('Mean HMC acceptance ratio: %.3f\n', mean(accRatio));
fprintf('Best model log-likelihood: %.3f\n', best_model_like);

%% Autocorrelation

h1 = autocorrelationplot(models);

%% Chain plots

h2 = chainplot(models,var_names,prior_range,mtest);

%% Corner plot of parameters

h3 = ecornerplot(models,'ks',true,'color',[.3 .3 .3],'name',var_names,'bestmodel',best_model,'truevals',mtest);

%% Barplot of parameters

h4 = barplot_parameters(models,var_names,prior_range,'bestmodel',best_model,'truevals',mtest);

%% Comparison best model and data

h5 = conc_modelledVSobserved(best_pred,testObs(1:n),sigmaObs(1:n),testObs(n+1:end),sigmaObs(n+1:end));

%% Export

if export
    exportgraphics(h1,['./output/' filetag '_' scenario '_hmc_autocorrelation.png'],'Resolution',300)
    exportgraphics(h2,['./output/' filetag '_' scenario '_hmc_chains.png'],'Resolution',300)
    exportgraphics(h3,['./output/' filetag '_' scenario '_hmc_cornerplot.png'],'Resolution',300)
    exportgraphics(h4,['./output/' filetag '_' scenario '_hmc_barplot.png'],'Resolution',300)
    exportgraphics(h5,['./output/' filetag '_' scenario '_hmc_datafit.png'],'Resolution',300)

    save(['./output/' filetag '_' scenario '_hmc_workspace.mat'], ...
        'models', 'logLikeStore', 'best_model', 'best_pred', 'mtest', ...
        'testObs', 'sigmaObs', 'prior_range', 'var_names', 'accRatio', ...
        'hmcSamplers', 'tuningInfo', 'endPoints', 'scenario');
end


function u = bounded_to_unbounded(m, lower_bounds, upper_bounds)
span = upper_bounds - lower_bounds;
z = (m - lower_bounds) ./ span;
z = min(max(z, 1e-10), 1 - 1e-10);
u = log(z ./ (1 - z));
end


function m = bounded_from_unbounded(u, lower_bounds, upper_bounds)
span = upper_bounds - lower_bounds;
s = 1 ./ (1 + exp(-u));
m = lower_bounds + span .* s;
end


function lpdf = transformed_logpdf(u, lower_bounds, upper_bounds, logLike)
m = bounded_from_unbounded(u, lower_bounds, upper_bounds);
span = upper_bounds - lower_bounds;

log_sigmoid = -softplus(-u);
log_one_minus_sigmoid = -softplus(u);
log_jacobian = sum(log(span) + log_sigmoid + log_one_minus_sigmoid);

lpdf = logLike(m) + log_jacobian;
end


function y = softplus(x)
y = log1p(exp(-abs(x))) + max(x,0);
end