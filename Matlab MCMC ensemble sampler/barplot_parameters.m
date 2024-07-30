function H=barplot_parameters(m,varnames,varargin)
% makes a barplot of of MCMC chain parameters. Additionally you can plot
% the best model as varargin{1}, and the true values if this is a test case
% as varargin{2}.
% Richard Ott, 2024
p = inputParser;
addParameter(p, 'bestmodel',[]);
addParameter(p, 'truevals',[]);
parse(p,varargin{:});


models= [];
for i = 1:length(m)
    models = [ models, m{i}.u(:,m{i}.status == 1)];
end

H=figure();
boxchart(models','Orientation','horizontal')
yticklabels(varnames)
legendlabels{1} = 'posterior';

if ~isempty(p.Results.truevals)
    hold on 
    plot(p.Results.truevals, 1:length(p.Results.truevals), 'p', 'MarkerSize', 10, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r')
    legendlabels = [legendlabels ; {'true value'}];
end

if ~isempty(p.Results.bestmodel)
    hold on 
    plot(p.Results.bestmodel, 1:length(p.Results.bestmodel), 'd', 'MarkerSize', 8, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'none')
    legendlabels = [legendlabels ; {'best fit'}];
end


legend(legendlabels)

end
