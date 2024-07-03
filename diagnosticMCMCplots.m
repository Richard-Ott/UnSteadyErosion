function [h1,h2,h3,h4,h5] = diagnosticMCMCplots(models,logLike,forward_model, var_names,prior_range,varargin)
% This function produces several diagnostic and analysis plots for the MCMC
% varargin -  test data with known true values can be passed on this way.
% Richard Ott, 2024

if nargin > 3
    mtrue = varargin{1};
end

%% Autocorrelation

h1 = autocorrelationplot(models);

%% Chain plots

h2 = chainplot(models,var_names,prior_range,mtrue);

%% Corner plot of parameters

h3 = ecornerplot(models,'ks',true,'color',[.3 .3 .3],'name',var_names,"truevals",mtrue);

%% Barplot of parameters

h4 = barplot_parameters(models,var_names,mtrue);

%% Best-fit model
posterior_like = squeeze(logLike(2,:,:));

[best_walker_like, best_walker_index] = max(posterior_like,[],2);
[best_model_like, best_index] = max(best_walker_like);

best_model = models(:,best_index,best_walker_index(best_index));

%% Comparison best model and data

best_pred = forward_model(best_model);
difference = testObs - best_pred

h5 = conc_modelledVSobserved(best_pred,testObs(1:n),testObs(1:n).*0.08,testObs(n+1:end),testObs(n+1:end)*0.08);


end