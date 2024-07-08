% This code computes the isolines of erosion ratios for a step
% change in erosion rate, given 10Be-14C nuclice concentration measurements.
% Richard Ott, 2024
clc
clear
close all
addpath('./online-calculators-v3/')
addpath('./data/')
addpath('./CosmoTools/')

%% load data ------------------------------------------------------------ %

SAMS = cosmosampleread('data\WCdata_RFO.xlsx');
DEM  = GRIDobj('.\data\crete_clipped_utm.tif');
SAMS(6:end) = []; % dont overload the plot

t = 1e1:1e2:1e4;        % time range of step change 

%% calculate isolines
[E1,E2,E1up,E1low,E2up,E2low] = calc_isoline(SAMS,DEM,t,'step');

%% Plot isolines

figure()
cc = viridis(numel(SAMS));
counter = 1;
for i = 1:numel(SAMS)
    plot(t(2:end),E2(i,2:end)./E1(i,2:end),'-','Color',cc(counter,:),'LineWidth',1.5)
    hold on
    indsUp = ~isnan(E2up(i,2:end)./E1up(i,2:end)) & ~isinf(E2up(i,2:end)./E1up(i,2:end));
    indsLow = ~isnan(E2low(i,2:end)./E1low(i,2:end)) & ~isinf(E2low(i,2:end)./E1low(i,2:end));
    f = fill([t([false indsUp]),fliplr(t([false indsLow]))],[E2up(i,[false indsUp])./E1up(i,[false indsUp]), fliplr(E2low(i,[false indsLow]))...
        ./fliplr(E1low(i,[false indsLow] ))] ,cc(counter,:),'FaceAlpha',0.3,'EdgeColor','none');
    set(get(get(f,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    hold on
    counter = counter+1;
end
xlabel('years since step change')
ylabel('E2/E1')
title('Erosion rate ratios')
legend(SAMS.ID)
ylim([0,500])   % for some times the ratios will get quiet crazy, so good to use limits
