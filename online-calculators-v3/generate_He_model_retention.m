function out = generate_He_model_retention()

% This generates model retention curves for cosmogenic He-3 in quartz using
% the Wolf and Farley formula. Returns a t vector and a matrix of retention
% vectors for various temperatures. Note assumes grain size and diffusion
% kinetics. 

%% First, figure diffusion kinetics

% Extract diffusion kinetics from Shuster and Farley irradiation

% Experimental data for He-3 in quartz

Ea = 84.5; % kJ/mol; 
lnD0a2 = 11.1; % ln(1/s); 

% Want ln(D0), i.e. non-radius-normalized units

expa = 0.0215; % radius of S&F analysed grain (cm)
D0a2 = exp(lnD0a2); %1/s 
D0 = D0a2.*(expa.^2); % cm2/s
lnD0 = log(D0); % ln(cm2/s)  


%% Calculate retention as function of hold time for various temperatures

d3.Ea = Ea; d3.lnD0 = lnD0;
% Specify grain size
rmm = 0.25; % That is the radius in mm
d3.r = rmm.*0.1; % cm

% time vector
tv = [1 logspace(2,7,40)];
% temps 
TCs = [-40 -30 -20 -10 0];

% Loop 
for row = 1:length(TCs);
    for col = 1:length(tv);
        d3.t_yr = tv(col);
        d3.TC = TCs(row);
        R(row,col) = WF_retention(d3);
    end;
end;

% Return

out.tv = tv; out.R = R; out.TCs = TCs;


