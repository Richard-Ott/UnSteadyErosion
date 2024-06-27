clear
clc
close all

addpath('.\online-calculators-v3\')
addpath('.\Matlab MCMC ensemble sampler\')
addpath('.\CosmoTools\')

data = readtable('data\WCdata_RFO.xlsx'); % AMS data
nWalks = 50;                              % how many chains per sample?

%% outline basins for binning (currently CosmoTools isnt properly integrated, should be improved)

SAMS = cosmosampleread('data\WCdata_RFO.xlsx');
DEM  = GRIDobj('.\data\crete_clipped_utm.tif');
SAMS = cosmowatersheds(SAMS,DEM);

% get median lat, lon, altitude for production calculation
lat = arrayfun(@(x) median(x.WSLat), SAMS);
lon = arrayfun(@(x) median(x.WSLon), SAMS);
alt = arrayfun(@(x) median(x.WSDEM.Z(:),'omitnan'), SAMS);

%% Observations

Nobs = [data.N10 ; data.N14];
dNobs= [data.N10sigma; data.N14sigma];
 
%% Priors -----------------------------------------------------------------
T =  [1,10e3];      % time of step change in yrs [min,max]
E1 = [10,3e2];      % old erosion rate in mm/ka  [min,max]
CHG = [0.1 200];    % increase [ ] 

% calculate prior ranges
n = length(Nobs)/2; % number of samples
prior_range = [T; repmat(E1, n, 1); CHG];
var_names = ['T1', arrayfun(@(x) sprintf('E1_sample%d', x), 1:n, 'UniformOutput', false),'ChangeFactor'];

%% Constants

[consts,Nmu] = make_constants();

%% Production rates

sp = Cronus_v3_spallation(lat,lon,alt,consts);   % get sample parameters (surface procution, pressure)

%% initial guess

mini  = initialmodel_flatprior(prior_range,nWalks);

%% Forward model

forward_model = @(m) Nforward_E_change_multisample(m(2:end-1),[m(1); 0],m(end),sp,consts,Nmu);

%% log likelihood function
% First we define a helper function equivalent to calling log(normpdf(x,mu,sigma))
% but has higher precision because it avoids truncation errors associated with calling
% log(exp(xxx)).
lognormpdf = @(x,mu,sigma)-0.5*((x-mu)./sigma).^2  -log(sqrt(2*pi).*sigma);
logLike    = @(m) sum(lognormpdf(Nobs, forward_model(m), dNobs));

logical_prior = @(m) sum(and(m > prior_range(:,1), m < prior_range(:,2))) == size(prior_range,1);

%% Posterior sampling
tic
[models, logLike] = gwmcmc(mini,{logical_prior logLike},1e8,'ThinChain',10,'burnin',.2,'StepSize',5);
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

h2 = chainplot(models,var_names,prior_range);

%% Corner plot of parameters

% figure
% ecornerplot(out,'ks',true,'color',[.3 .3 .3])

%% Barplot of parameters

h3 = barplot_parameters(models,var_names);

%% Best-fit model
posterior_like = squeeze(logLike(2,:,:));

[best_walker_like, best_walker_index] = max(posterior_like,[],2);
[best_model_like, best_index] = max(best_walker_like);

best_model = models(:,best_index,best_walker_index(best_index));


%% Comparison best model and data

best_pred = forward_model(best_model);
difference = Nobs - best_pred

%%
exportgraphics(h2,'testMCMC_chains.png','Resolution',300)
exportgraphics(h3,'testMCMC_barplot.png','Resolution',300)

