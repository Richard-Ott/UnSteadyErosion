function H=barplot_parameters(m,varnames,varargin)
% makes a barplot of of MCMC chain parameters. Additionally you can supple
% an array of the true values (if this is a test run) as varargin.
% Richard Ott, 2024

if nargin > 2
    truevals = varargin{1};
end

m=m(:,:);   % reshape and collapse all the individual chains together

H=figure();
boxchart(m','Orientation','horizontal')
yticklabels(varnames)

if nargin > 2
    hold on 
    plot(truevals, 1:length(truevals), 'p', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r')
end
end
