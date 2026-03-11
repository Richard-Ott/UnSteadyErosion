clear
clc
close all

addpath('.\online-calculators-v3\')
addpath('.\Matlab MCMC ensemble sampler\')
addpath('.\CosmoTools\')

data = readtable('data\WCdata_RFO.xlsx'); % AMS data
nWalks = 5;                              % how many chains per sample?

scenarios = {'step', 'spike'}; 

nsteps = 1;
export = 1;

%% outline basins for binning (currently CosmoTools isnt properly integrated, should be improved)

SAMS = cosmosampleread('data\WCdata_RFO.xlsx');
DEM  = GRIDobj('.\data\crete_clipped_utm.tif');
SAMS = cosmowatersheds(SAMS,DEM);

% get median lat, lon, altitude for production calculation
lat = arrayfun(@(x) median(x.WSLat), SAMS);
lon = arrayfun(@(x) median(x.WSLon), SAMS);
alt = arrayfun(@(x) median(x.WSDEM.Z(:),'omitnan'), SAMS);

%% Observations

Nobs = [data.N10 ; data.N14; data.N26];
dNobs= [data.N10sigma; data.N14sigma; data.N26sigma];

% cannot accept solution because the uncerts on these two samples are so
% low
dNobs(15) = dNobs(15)*3;
dNobs(16) = dNobs(16)*2;

Nlogical = [~isnan(data.N10) ~isnan(data.N14) ~isnan(data.N26)];

for i = 1

if i ~= 1; nWalks = 20;end   % long chains for step  
%% Priors -----------------------------------------------------------------
T   = [1,6e3];      % time of step change in yrs [min,max]
if i == 1; E = [10,5e3]; else; E = [10,3e2]; end    % range of expected erosion rates in mm/ka  [min,max]
CHG = [0.1 100];     % increase [ ] 
LOSS = [0,250];     % loss of soil in cm [min,max], can be commented if no spike model

% calculate prior ranges
[prior_range,var_names] = make_prior_and_varnames(scenarios{i},T,E,LOSS,CHG,length(data.N10),nsteps);
prior_range(2:11,2) = 300;
% make initial erosion rates lower
%% Constants

consts = make_constants();

%% Production rates

sp = sample_parameters(lat,lon,alt,consts);   % get sample parameters (surface procution, pressure)

%% initial guess

mini  = initialmodel_flatprior(prior_range,nWalks);

%% Forward model

forward_model = @(m) Nforward_wrapper(m,sp,consts,scenarios{i},nsteps,Nlogical);

%% log likelihood function
% First we define a helper function equivalent to calling log(normpdf(x,mu,sigma))
% but has higher precision because it avoids truncation errors associated with calling
% log(exp(xxx)).
lognormpdf = @(x,mu,sigma)-0.5*((x-mu)./sigma).^2  -log(sqrt(2*pi).*sigma);
logLike    = @(m) sum(lognormpdf(Nobs(Nlogical), forward_model(m), dNobs(Nlogical)));

logical_prior = @(m) sum(and(m > prior_range(:,1), m < prior_range(:,2))) == size(prior_range,1);

%% Posterior sampling
tic
[models, logLike] = gwmcmc(mini,{logical_prior logLike},1.5e9,'ThinChain',5,'burnin',.3,'StepSize',3);
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

h2 = chainplot(models,var_names,prior_range);

%% Corner plot of parameters

h3 = ecornerplot(models,'ks',true,'color',[.3 .3 .3],'name',var_names,'bestmodel',best_model);

%% Barplot of parameters

h4 = barplot_parameters(models,var_names,prior_range,'bestmodel',best_model);

%% Comparison best model and data

h5 = conc_modelledVSobserved(best_pred,data.N10,data.N10sigma,data.N14,data.N14sigma);

%% Export

if export
    exportgraphics(h1,['./output/WC/WC_' scenarios{i} '_longchain_autocorrelation.png'],'Resolution',300)
    exportgraphics(h2,['./output/WC/WC_' scenarios{i} '_longchain_chains.png'],'Resolution',300)
    exportgraphics(h3,['./output/WC/WC_' scenarios{i} '_longchain_cornerplot.png'],'Resolution',300)
    exportgraphics(h4,['./output/WC/WC_' scenarios{i} '_longchain_barplot.png'],'Resolution',300)
    exportgraphics(h5,['./output/WC/WC_' scenarios{i} '_longchain_datafit.png'],'Resolution',300)

    exportgraphics(h1,['./output/WC/WC_' scenarios{i} '_longchain_autocorrelation.png'],'ContentType','vector')
    exportgraphics(h2,['./output/WC/WC_' scenarios{i} '_longchain_chains.png'],'ContentType','vector')
    exportgraphics(h3,['./output/WC/WC_' scenarios{i} '_longchain_cornerplot.png'],'ContentType','vector')
    exportgraphics(h4,['./output/WC/WC_' scenarios{i} '_longchain_barplot.png'],'ContentType','vector')
    exportgraphics(h5,['./output/WC/WC_' scenarios{i} '_longchain_datafit.png'],'ContentType','vector')

    clear h1 h2 h3 h4 h5
    save(['./output/WC/WC_workspace_' scenarios{i} '.mat'],  '-v7.3')
end
disp([scenarios{i} ' completed :)'])
clear
close all
clear 'consts' 'sp' 'mini' 'models' 'logLike' 'var_names' 'prior_range' 'h1' 'h2' 'h3' 'h4' 'h5' 'C' 'lags' 'ESS' 'best_model' 'best_pred' 'posterior_like'
end