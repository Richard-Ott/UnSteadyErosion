function [sp] = Cronus_v3_spallation(lat,lon,altitude,consts)
% get top spallation production from Cronus v3 for 10Be, 14C, 26Al
% multi sample input possible (provide variables as vectors)
% Richard Ott, 2024

addpath('C:\Users\rott\OneDrive - UvA\Richard\PhD_ETH\matlab\CRONUS cosmo calculation\cronus 3.0\online-calculators-v3')
addpath('C:\Users\rott\OneDrive - UvA\Richard\PhD_ETH\matlab\CRONUS cosmo calculation\cronus 3.0\online-calculators-v3\data')

sp.pressure = ERA40atm(lat, lon, altitude);  % atmospheric pressure

if isrow(sp.pressure)   % make sure its a column vector
    sp.pressure = sp.pressure';
end

P10_SLHL = consts.refP_St(4);                % Stone scaling ref spallation production
P14_SLHL = consts.refP_St(5);
P26_SLHL = consts.refP_St(7);

sf_St = (stone2000(lat,sp.pressure,1))';        % scaling factor

sp.P10spal = sf_St .* P10_SLHL;              % site specific spallation production
sp.P14spal = sf_St .* P14_SLHL;
sp.P26spal = sf_St .* P26_SLHL;

end