function out = WF_retention(in)

% This is from Wolf and Farley paper
% 
% Syntax: out = WF_retention(in);
%
% where:
%  in.TC temp (degC)
%  in.r grain radius (cm)
%  in.Ea = activation energy kJ/mol
%  in.lnD0 ln(1/s)
%  in.t_yr = holding time (yr)
%
% out is retention - in Wolf and Farley this is age -- apparent age/true
% age...this is the same as amount of diffusant present/amount produced. 


T = in.TC + 273.15; % deg K

Ea = in.Ea.*1e3; % J/mol; 
lnD0 = in.lnD0; % ln(cm2/s); 

R = 8.314; % J/mol K

a = in.r; % model grain radius cm

d0 = exp(lnD0); % cm2/s
D = d0.*exp(-Ea./(R.*T)); % cm2/s
Da2 = D./(a.^2); % 1/s

sPerYear = 60*60*24*365.25;

t = in.t_yr.*sPerYear;
  
maxn = 2048; % summation length
n = 1:maxn;

toSum = (6./((pi.^4).*(n.^4))) .* exp(-(n.^2).*(pi.^2).*Da2.*t);

tprime_s = (1./Da2).*(1/15 - sum(toSum));

tprime_yr = tprime_s./sPerYear;

out = tprime_yr./in.t_yr;

