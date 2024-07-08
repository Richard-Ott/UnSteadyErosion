function H=barplot_parameters(m,varnames,varargin)
% makes a barplot of of MCMC chain parameters. Additionally you can plot
% the best model as varargin{1}, and the true values if this is a test case
% as varargin{2}.
% Richard Ott, 2024
p = inputParser;
addParameter(p, 'bestmodel',[]);
addParameter(p, 'truevals',[]);
parse(p,varargin{:});

m=m(:,:);   % reshape and collapse all the individual chains together

H=figure();
boxchart(m','Orientation','horizontal')
yticklabels(varnames)

if ~isempty(p.Results.truevals)
    hold on 
    plot(p.Results.truevals, 1:length(p.Results.truevals), 'p', 'MarkerSize', 10, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r')
end

if ~isempty(p.Results.bestmodel)
    hold on 
    plot(p.Results.bestmodel, 1:length(p.Results.bestmodel), 'd', 'MarkerSize', 8, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'none')
end

end
