clear
clc
close all

addpath('..\..\..\CRONUS cosmo calculation\cronus 3.0\online-calculators-v3\')
addpath('..\..\..\..\data\geochronology\cosmo\in_situ_14C\')
addpath('..\..\..\Richards functions\Matlab MCMC ensemble sampler\')

data = readtable('14Crete_Richard.xlsx'); % AMS data
data([5,12],:) = [];                      % remove unwanted data
nWalks = 50;                              % how many chains per sample?

%% Observations

Nobs = [data.N10 ; data.N14_atG_1_];
dNobs= [data.dN10; data.x_N14_atG_1_];
 
%% Priors -----------------------------------------------------------------
T =  [1,10e3];      % time of step change in yrs [min,max]
E1 = [10,3e2];      % old erosion rate in mm/ka  [min,max]
E2 = [100,10e3];    % new erosion rates in mm/ka [min,max]

% calculate prior ranges
prior_range = repmat([E1; E2], height(data), 1);
prior_range = [T; prior_range];
erosion_variables = arrayfun(@(i) {sprintf('E1_samp%d', i), sprintf('E2_samp%d', i)}, 1:length(data.N10), 'UniformOutput', false);
var_names = ['T1' , horzcat(erosion_variables{:})];

%% Constants

[consts,Nmu] = make_constants();

%% Production rates

lat= [30,30];
lon= [0,0];
altitude=[500,500];

sp = Cronus_v3_spallation(lat,lon,altitude,consts);   % get sample parameters (surface procution, pressure)

%% initial guess

mini  = initialmodel_flatprior(prior_range,nWalks);

%% Forward model

forward_model = @(m) Nforward_E_discretized_multisample(m(3:end),[m(1:2); 0],sp,consts,Nmu);

%% generate test data to see if inversion can succesfully identify these data
t1 = [5000,100];
e1 = [200,200,200];
e2 = [200,200,200];
% mtest = [t1; repmat([e1;e2], 10,1) ];
mtest = [t1'; e1'; e2'];
testObs = forward_model(mtest);

%% log likelihood function
% First we define a helper function equivalent to calling log(normpdf(x,mu,sigma))
% but has higher precision because it avoids truncation errors associated with calling
% log(exp(xxx)).
lognormpdf = @(x,mu,sigma)-0.5*((x-mu)./sigma).^2  -log(sqrt(2*pi).*sigma);
logLike    = @(m) sum(lognormpdf(testObs, forward_model(m), testObs*0.08));

logical_prior = @(m) sum(and(m > prior_range(:,1), m < prior_range(:,2))) == size(prior_range,1);

%% Posterior sampling
tic
[models, logLike] = gwmcmc(mini,{logical_prior logLike},1e8,'ThinChain',2,'burnin',.2,'StepSize',2);
toc

%% Autocorrelation

figure
[C,lags,ESS]=eacorr(models);
plot(lags,C,'.-',lags([1 end]),[0 0],'k');
grid on
xlabel('lags')
ylabel('autocorrelation');
text(lags(end),0,sprintf('Effective Sample Size (ESS): %.0f_ ',ceil(mean(ESS))),'verticalalignment','bottom','horizontalalignment','right')
title('Markov Chain Auto Correlation')

%% remove "bad" chains

meanC = mean(C,1); % chains that got stuck with high autocorrelation

thres = 0.2; % mean autocorrelation of chains that is allowed. The bad chains have values around 0.4-0.5

remove_inds = meanC > thres;

models(:,remove_inds,:) = [];   % remove bad chains
logLike(:,remove_inds,:)= [];


%% Chain plots

chainplot(models,var_names,prior_range)

%% Corner plot of parameters

% figure
% ecornerplot(out,'ks',true,'color',[.3 .3 .3])

%% Barplot of parameters

barplot_parameters(models,var_names)

%% Concentration observed vs modelled



%% Best-fit model
posterior_like = squeeze(logLike(2,:,:));

[best_walker_like, best_walker_index] = max(posterior_like,[],2);
[best_model_like, best_index] = max(best_walker_like);

best_model = models(:,best_index,best_walker_index(best_index));


%% Comparison best model and data

best_pred = forward_model(best_model);
difference = testObs - best_pred

