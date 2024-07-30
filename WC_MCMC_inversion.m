clear
clc
close all

addpath('.\online-calculators-v3\')
addpath('.\Matlab MCMC ensemble sampler\')
addpath('.\CosmoTools\')

data = readtable('data\WCdata_RFO.xlsx'); % AMS data
nWalks = 30;                              % how many chains per sample?

scenario = 'samestep'; 

nsteps = 1;

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

Nlogical = [~isnan(data.N10) ~isnan(data.N14) ~isnan(data.N26)];

%% Priors -----------------------------------------------------------------
T   = [1,10e3];      % time of step change in yrs [min,max]
E1  = [10,3e2];      % old erosion rate in mm/ka  [min,max]
CHG = [0.1 100];     % increase [ ] 
LOSS = [0,200];     % loss of soil in cm [min,max], can be commented if no spike model

% calculate prior ranges
[prior_range,var_names] = make_prior_and_varnames(scenario,T,E1,LOSS,CHG,length(data.N10),nsteps);

%% Constants

[consts,Nmu] = make_constants();

%% Production rates

sp = Cronus_v3_spallation(lat,lon,alt,consts);   % get sample parameters (surface procution, pressure)

%% initial guess

mini  = initialmodel_flatprior(prior_range,nWalks);

%% Forward model

forward_model = @(m) Nforward_wrapper(m,sp,consts,Nmu,scenario,nsteps,Nlogical);

%% log likelihood function
% First we define a helper function equivalent to calling log(normpdf(x,mu,sigma))
% but has higher precision because it avoids truncation errors associated with calling
% log(exp(xxx)).
lognormpdf = @(x,mu,sigma)-0.5*((x-mu)./sigma).^2  -log(sqrt(2*pi).*sigma);
logLike    = @(m) sum(lognormpdf(Nobs(Nlogical), forward_model(m), dNobs(Nlogical)));

logical_prior = @(m) sum(and(m > prior_range(:,1), m < prior_range(:,2))) == size(prior_range,1);

%% Posterior sampling
%% Posterior sampling

walkers = Egholm_MCMC(nWalks,Nobs,dNobs,mini,prior_range,forward_model,Nlogical);

%% best model
best_obs_err = inf;
for j = 1:nWalks
    [best_obs_err_walker, best_ind] = min(walkers{j}.restot);
    best_model_walker = walkers{j}.u(:,best_ind);
    if best_obs_err_walker < best_obs_err
        best_model = best_model_walker;
    end
end

%% Chain plots

h1 = chainplot(walkers,var_names,prior_range);

%% Corner plot

h2 = ecornerplot(walkers,'ks',true,'color',[.3 .3 .3],'name',var_names,'bestmodel',best_model,'support',prior_range');

%% Barplot of parameters

h3 = barplot_parameters(walkers,var_names,'bestmodel',best_model);

%% Comparison best model and data

h4 = conc_modelledVSobserved(walkers,data.N10,data.N10sigma,data.N14,data.N14sigma);

%% Export

if export
    exportgraphics(h1,['./output/WC_' scenarios{i} '_chains.png'],'Resolution',300)
    exportgraphics(h2,['./output/WC_' scenarios{i} '_cornerplot.png'],'Resolution',300)
    exportgraphics(h3,['./output/WC_' scenarios{i} '_barplot.png'],'Resolution',300)
    exportgraphics(h4,['./output/WC_' scenarios{i} '_datafit.png'],'Resolution',300)
end
