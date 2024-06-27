function out = chisquarecdf_rt(x,v)

% Returns right-tailed probability from chi-squared CDF
% uses gammainc.m, which is in normal MATLAB, to avoid use of stats toolbox
% 
% syntax: out = chisquarecdf_rt(x,DOF)
% where x is the chi-squared (not the reduced chi-squared) value, and DOF
% is the degrees of freedom. 

out = 1-gammacdf(x,v/2,2); % defn of chi-squared cdf, 1- makes it right-tailed

function p = gammacdf(x,a,b)

z = x ./ b;
p = gammainc(z,a);

p(z == Inf) = 1;