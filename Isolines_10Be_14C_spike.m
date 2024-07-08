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
[E1,loss,E1up,E1low,lossup,losslow] = calc_isoline(SAMS,DEM,t,'spike');

%% Plot isolines 
% the search boundaries I'm setting and numerical instabilities lead to
% weird values typically at small and large times. Currently, I just rmeove
% these areas with xlim and ylim. You can also choose wieder search
% boundaries in the calc_isoline function (will make it slower).
figure()
cc = viridis(numel(SAMS));
for i = 1:numel(SAMS)
    subplot(1,2,1)
    plot(t,loss(i,:),'-','Color',cc(i,:),'LineWidth',1.5)
    ylim([0,200])
    xlim([0,5000])
    hold on
    subplot(1,2,2)
    plot(t,E1(i,:),'-','Color',cc(i,:),'LineWidth',1.5)
    hold on
    ylim([0,100])
    xlim([0,5000])
end
legend(SAMS.ID)

