function out = billipse(xin,delxin,yin,delyin,plotSelect,plotString);

% Plots 68% and 95% confidence ellipses for a bivariate normal
% distribution.
% 
% Syntax: out = billipse(xin,dxin,yin,dyin,plotSelect,plotString);
% 
% where xin and dxin, and dxin and dyin, are the means and standard errors
% for the two measurements. 
%
% NOTE: assumes uncorrelated measurements. This should be fixed later. 
%
% plotSelect is 1 for 68% ellipse and 2 for both 68% and 95% ellipses.
% Defaults to 1.
% enter 0 to return x and y vectors -- in this case you get out.x and out.y
% 
% plotstring is optional; allows color specification, e.g., 'b'.
%
% Returns the function handles
%
% Greg Balco -- UW Cosmogenic Nuclide Lab -- February, 2007
%
% Accompanies G. Balco and C. Rovey, 'An isochron method for cosmogenic-
% nuclide dating of buried soils and sediments,' for publication in the 
% American Journal of Science.
% 
% Developed in part with funding from the National Science Foundation.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License, version 2,
% as published by the Free Software Foundation (www.fsf.org).


% check args


if nargin < 5
	plotSelect = 1;
end;

% mesh

[x,y] = meshgrid((xin-4*delxin):(0.1*delxin):(xin+4*delxin),...
		(yin-4*delyin):(0.1*delyin):(yin+4*delyin));

% calculate PDF...note assumes uncorrelated

Prob = exp(-0.5 .* (((x-xin)./delxin).^2 + ((y-yin)./delyin).^2));
		
% So we now know the probability for each cell in the grid. 
% Now find the 68% probability contour. 

% Normalize to volume = 1

normP = Prob ./ sum(sum(Prob));

% Now we need to figure out cumulative probabilities. 

% multiply by 10000 to achieve manageable number of values, round to get integers:

normP = normP * 10000;
intP = round(normP);

for a = 1:max(max(intP));
	cumprob(a) = sum(intP(find(intP >= a)))/10000;
end;

probs = 1:max(max(intP));

sigma1 = find(abs(cumprob - 0.68) == min(abs(cumprob - 0.68)));
sigma2 = find(abs(cumprob - 0.95) == min(abs(cumprob - 0.95)));

% weed out cases where adjacent probs are same -- rounding error

if length(sigma1) ~= 1;
	sigma1 = min(sigma1);
end;

if length(sigma2) ~= 1;
	sigma2 = min(sigma2);
end;

% Now draw the contours.

cmat = contourc(x(1,:),y(:,1),normP,[sigma1 sigma2]);

cl1 = cmat(2,1);
s = size(cmat);

% Sometimes contourc returns a bunch of contours for one level -- grid size issue?
% This is spurious, so plot only the major one. 

contourStarts = find(cmat(1,:) == sigma1);
contourSizes = cmat(2,contourStarts);
contourToPlot = find(contourSizes == max(contourSizes));

x1 = cmat(1,(contourStarts(contourToPlot)+1):(contourStarts(contourToPlot) + contourSizes(contourToPlot)));
y1 = cmat(2,(contourStarts(contourToPlot)+1):(contourStarts(contourToPlot) + contourSizes(contourToPlot)));

if plotSelect == 0;
    out.x = x1;
    out.y = y1;
    return;
end;

if nargin < 6;
	plotString = 'k';
end;

sd1 = plot(x1,y1,plotString);

if plotSelect == 2;
	
	hold on;
	
	contourStarts = find(cmat(1,:) == sigma2);
	contourSizes = cmat(2,contourStarts);
	contourToPlot = find(contourSizes == max(contourSizes));

	x2 = cmat(1,(contourStarts(contourToPlot)+1):(contourStarts(contourToPlot) + contourSizes(contourToPlot)));
	y2 = cmat(2,(contourStarts(contourToPlot)+1):(contourStarts(contourToPlot) + contourSizes(contourToPlot)));

	sd2 = plot(x2,y2,plotString);
end;