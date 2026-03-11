clear
clc
close all

addpath('.\online-calculators-v3\')
addpath('.\Matlab MCMC ensemble sampler\')
addpath('.\CosmoTools\')

nWalks = 10;              % number of independent HMC chains
scenarios = {'step', 'spike'};
nsteps = 1;
export = 1;

useTuning = strcmp(getenv('WC_HMC_TUNE'), '1');  % set env var WC_HMC_TUNE=1 to enable tuning
numTune = 20;            % tuneSampler iterations when enabled
numBurnin = 5000;        % warmup draws per chain
numSamples = 200000;      % retained draws per chain
thinSize = 2;            % thinning for retained chain

stepSize = 0.01;         % manual HMC step size
numLeapfrogSteps = 10;   % manual number of leapfrog steps

%% Outline basins for binning
SAMS = cosmosampleread('data\WCdata_RFO.xlsx');
DEM = GRIDobj('.\data\crete_clipped_utm.tif');
SAMS = cosmowatersheds(SAMS, DEM);

% Median basin coordinates and elevation for production calculations.
lat = arrayfun(@(x) median(x.WSLat), SAMS);
lon = arrayfun(@(x) median(x.WSLon), SAMS);
alt = arrayfun(@(x) median(x.WSDEM.Z(:), 'omitnan'), SAMS);

for i = 1:numel(scenarios)
    data = readtable('data\WCdata_RFO.xlsx');

    %% Observations
    Nobs = [data.N10; data.N14; data.N26];
    dNobs = [data.N10sigma; data.N14sigma; data.N26sigma];

    % Inflate two very small analytical uncertainties to keep residuals balanced.
    dNobs(15) = dNobs(15) * 3;
    dNobs(16) = dNobs(16) * 2;

    Nlogical = [~isnan(data.N10) ~isnan(data.N14) ~isnan(data.N26)];

    %% Priors
    T = [1, 6e3];
    if i == 1
        E = [10, 5e3];
    else
        E = [10, 3e2];
    end
    CHG = [0.1, 100];
    LOSS = [0, 250];

    [prior_range, var_names] = make_prior_and_varnames(scenarios{i}, T, E, LOSS, CHG, length(data.N10), nsteps);
    
    % make initial erosion rate low for step scenario
    if i == 1
        prior_range(2:11, 2) = 300;
    end

    %% Constants and production rates
    consts = make_constants();
    sp = sample_parameters(lat, lon, alt, consts);

    %% Initial guess inside prior support
    mini = initialmodel_flatprior(prior_range, nWalks);

    %% Forward model and log likelihood
    forward_model = @(m) Nforward_wrapper(m, sp, consts, scenarios{i}, nsteps, Nlogical);

    lognormpdf = @(x, mu, sigma) -0.5 * ((x - mu) ./ sigma).^2 - log(sqrt(2 * pi) .* sigma);
    logLike = @(m) sum(lognormpdf(Nobs(Nlogical), forward_model(m), dNobs(Nlogical)));

    %% Unconstrained transform for HMC
    lower_bounds = prior_range(:, 1);
    upper_bounds = prior_range(:, 2);

    nParams = size(prior_range, 1);
    models = nan(nParams, nWalks, numSamples);
    logLikeStore = nan(2, nWalks, numSamples);
    accRatio = nan(nWalks, 1);
    endPoints = nan(nParams, nWalks);
    hmcSamplers = cell(nWalks, 1);
    tuningInfo = cell(nWalks, 1);

    tic
    for wix = 1:nWalks
        wix
        startModel = mini(:, wix);
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
        models(:, wix, :) = reshape(chainM, nParams, 1, size(chainM, 2));
        endPoints(:, wix) = bounded_from_unbounded(endpointU, lower_bounds, upper_bounds);

        for six = 1:numSamples
            logLikeStore(1, wix, six) = 0;
            logLikeStore(2, wix, six) = logLike(models(:, wix, six));
        end
    end
    toc

    models = single(models);
    logLikeStore = single(logLikeStore);

    %% Best-fit model
    posterior_like = reshape(logLikeStore(2, :, :), nWalks, numSamples);
    [best_model_like, linear_index] = max(posterior_like(:));
    [best_index, best_sample_index] = ind2sub(size(posterior_like), linear_index);
    best_model = models(:, best_index, best_sample_index);
    best_pred = forward_model(double(best_model));

    fprintf('%s mean HMC acceptance ratio: %.3f\n', scenarios{i}, mean(accRatio));
    fprintf('%s best model log-likelihood: %.3f\n', scenarios{i}, best_model_like);

    %% Plots
    h1 = autocorrelationplot(models);
    h2 = chainplot(models, var_names, prior_range);
    h3 = ecornerplot(models, 'ks', true, 'color', [.3 .3 .3], 'name', var_names, 'bestmodel', best_model);
    h4 = barplot_parameters(models, var_names, prior_range, 'bestmodel', best_model);
    h5 = conc_modelledVSobserved(best_pred, data.N10, data.N10sigma, data.N14, data.N14sigma);

    %% Export
    if export
        if ~exist('./output/WC', 'dir')
            mkdir('./output/WC');
        end

        exportgraphics(h1, ['./output/WC/WC_' scenarios{i} '_hmc_autocorrelation.png'], 'Resolution', 300)
        exportgraphics(h2, ['./output/WC/WC_' scenarios{i} '_hmc_chains.png'], 'Resolution', 300)
        exportgraphics(h3, ['./output/WC/WC_' scenarios{i} '_hmc_cornerplot.png'], 'Resolution', 300)
        exportgraphics(h4, ['./output/WC/WC_' scenarios{i} '_hmc_barplot.png'], 'Resolution', 300)
        exportgraphics(h5, ['./output/WC/WC_' scenarios{i} '_hmc_datafit.png'], 'Resolution', 300)

        exportgraphics(h1, ['./output/WC/WC_' scenarios{i} '_hmc_autocorrelation.pdf'], 'ContentType', 'vector')
        exportgraphics(h2, ['./output/WC/WC_' scenarios{i} '_hmc_chains.pdf'], 'ContentType', 'vector')
        exportgraphics(h3, ['./output/WC/WC_' scenarios{i} '_hmc_cornerplot.pdf'], 'ContentType', 'vector')
        exportgraphics(h4, ['./output/WC/WC_' scenarios{i} '_hmc_barplot.pdf'], 'ContentType', 'vector')
        exportgraphics(h5, ['./output/WC/WC_' scenarios{i} '_hmc_datafit.pdf'], 'ContentType', 'vector')

        save(['./output/WC/WC_workspace_hmc_' scenarios{i} '.mat'], ...
            'models', 'logLikeStore', 'best_model', 'best_pred', 'prior_range', ...
            'var_names', 'accRatio', 'hmcSamplers', 'tuningInfo', ...
            'endPoints', 'scenarios', 'i', '-v7.3');

        clear h1 h2 h3 h4 h5
    end

    disp([scenarios{i} ' completed :)'])
    clear consts sp mini models logLikeStore var_names prior_range h1 h2 h3 h4 h5 best_model best_pred posterior_like
    close all
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
y = log1p(exp(-abs(x))) + max(x, 0);
end
