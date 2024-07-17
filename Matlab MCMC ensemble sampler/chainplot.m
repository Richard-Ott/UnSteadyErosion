function H = chainplot(models,var_names,prior_range,varargin)
% plot chains of MCMC sampling. In case you're running an inversion test
% you can supply the true parameter values as an array in varargin.
% Richard Ott, 2024

if nargin > 3
    true_vals = varargin{1};
end

nWalks = length(models);
nparameters = length(var_names);

H=figure('Units', 'normalized', 'Position', [0.08, 0.08, 0.85, 0.85]);
cc = lines(nWalks);
col1 = 0.9*[1,1,1];
col2 = 0.6*[1,1,1];


tiledlayout('flow');
nexttile
for i = 1:nparameters
	
    for j = 1:nWalks
        hold on
        I = find(models{j}.status == 0);
        plot(I,models{j}.up(i,I),'.','color',col1);
        I = find(models{j}.status == -1);
        plot(I,models{j}.up(i,I),'.','color',col2);
        I = find(models{j}.status == 1);
        plot(I,models{j}.up(i,I),'.','color',cc(j,:));
    end

    if nargin > 3; yline(true_vals(i),'LineWidth',2); end  % plot true values if available

    ylim([ prior_range(i,1), prior_range(i,2)] )
    ylabel(var_names(i), 'Interpreter', 'none')
    if i < nparameters; nexttile; end
end

end