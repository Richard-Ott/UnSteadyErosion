function [mean,stdev] = ewmean(values,errors);

% [Mean, Stdev] = ewmean(values,errors) returns the error-weighted
% 	mean of the data. values and errors have to be the same size. 
%
%	The error-weighted mean is Sigma(Wi * Xi) / Sigma(Wi)
%   where Xi are the values and Wi are the weights, such that
%   Wi = 1/(errori ^2).
%
%	The uncertainty is given by 1/sqrt(Sigma(Wi)).
%

weights = 1./(errors.^2);

mean = (sum(sum(weights.*values)))./(sum(sum(weights)));

stdev = 1./sqrt(sum(sum(weights)));

