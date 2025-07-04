function N = Nforward_discretized(E,T,sp,consts,Nmu,scenario,Nlogical,varargin)
% This function calculates concentrations N10 and N14 for mutiple step 
% changes in erosion and for multiple samples at once.
%
% Input:
%       - E: erosion rates in mm/ka from old to current for n samples for 
%           m step changes. Input format:
%               - n x m for type 'step'
%               - 1 x n for type 'samestep'
%               - scalar for type 'spike' (background erosion rate)
%        
%       - T: 1 x m vector of step change timings (yrs BP)
%           from old to young
%       - sp: sample specific parameters (Pspal & pressure)
%       - consts: the consts_v3.mat file from Cronus v3
%       - precalculated muon production rates from Cronus v3 (Nmu.mat)
%       - type: type of erosion scenario. Either: 
%               - 'step': step changes in erosion rates with individual
%               rates per sample
%               - 'samestep': step changes in erosion with common change
%               factor
%               - 'spike': spikes of soil loss
%       - varargin: empty for 'step'
%
% Output:
%       - [N10,N14]: concentrations at/g of 10Be and 14C
%
% Richard Ott & Dirk Scherler, 2024

% Parse optional inputs
p = inputParser;
addParameter(p, 'change_variable', []);
addParameter(p, 'mixing_depth', []);
parse(p, varargin{:});
CHG = p.Results.change_variable;
mix_depth = p.Results.mixing_depth;

nSamp = length(sp.P10spal);

E = E./1e4;   % convert to cm/a
T_time_spans   = [inf, diff(T')*(-1)];        % time span of every time interval between erosion changes

%% reshape erosion rate array and calculate depths at steps for different scenarios

switch scenario
    case 'step'
        segment_depths = E.*T_time_spans;            % the exhumation occuring during every erosion time interval in cm

    case 'samestep'
        % make matrix of erosion rates for different time steps by multiplying E
        % with the change factors
        E = [E, repmat(E,1,length(CHG))];
        E(:,2:end) = E(:,2:end) * diag(CHG);
        segment_depths = E.*T_time_spans;            % the exhumation occuring during every erosion time interval in cm

    case 'samebackground_step'
        E = repmat(E,nSamp,length(T));
        E(:,2:end) = E(:,2:end) .* CHG;
        segment_depths = E.*T_time_spans;            % the exhumation occuring during every erosion time interval in cm

    case 'samebackground_samestep'
        E = repmat(E,nSamp,length(T));
        E(:,2:end) = E(:,2:end) * diag(CHG);
        segment_depths = E.*T_time_spans;            % the exhumation occuring during every erosion time interval in cm

    case 'spike'
        E = repmat(E,1, length(T));                      % make matrix of erosion rates for easy calling in concentration loop
        segment_depths = E.*T_time_spans+[inf(nSamp,1) CHG];     % calcuate the exhumation occuring during every erosion time interval in cm

    case 'samespike'
        E = repmat(E,1, length(T));
        segment_depths = E.*T_time_spans+[inf(nSamp,1) repmat(CHG',nSamp,1)];     % calcuate the exhumation occuring during every erosion time interval in cm

    case 'samebackground_spike'
        E = repmat(E,nSamp,length(T));                      % make matrix of erosion rates for easy calling in concentration loop
        segment_depths = E.*T_time_spans+[inf(nSamp,1) CHG];

    case 'samebackground_samespike'
        E = repmat(E,nSamp,length(T));                      % make matrix of erosion rates for easy calling in concentration loop
        segment_depths = E.*T_time_spans+[inf(nSamp,1) repmat(CHG',nSamp,1)];
    case 'curve'
        scaleFactor = CHG;

        absolute_erosion_change = sp.curvechange .* scaleFactor; % calculate erosion rate change from relative differences from curve and scaling factor
        % make matrix of erosion rates for different time steps by multiplying E with the change factors
        E = repmat(E,1,length(T));
        E(:,2:end) = E(:,2:end)+ E(:,2:end) .* absolute_erosion_change;
        segment_depths = E.*T_time_spans;            % the exhumation occuring during every erosion time interval in cm
end

t_depths       = [fliplr(cumsum(fliplr(segment_depths(:,2:end)),2)), zeros(nSamp,1)];     % depths at every step change T in cm

% add depths for soil mixing
interval_spacing = 10;   % depth intervals for which to compute concentration
mixing_depths = 0:interval_spacing:mix_depth;

% Add mix_depth offsets to create new slices
if isempty(mix_depth) % no mixing
    nmix = 1;
    t_depths_3D = t_depths;

else % soil mixing depths
    nmix = length(mixing_depths);
    for i = 1:nmix
    t_depths_3D(:,:,i) = t_depths + mixing_depths(i);
    end
end



%% start cosmo calulcation

% assign variables --------------------------------------------------------
rho = consts.density;  

% attentuation lengths
att_l_10(1) = consts.L_sp;          % spallation
att_l_10(2) = consts.L_mu_atm(4);   % muon 
att_l_14(1) = consts.L_sp;          % spallation
att_l_14(2) = consts.L_mu_atm(5);   % muon 
att_l_26(1) = consts.L_sp;          % spallation
att_l_26(2) = consts.L_mu_atm(7);   % muon 

% production rates 
P10(:,1) = sp.P10spal; 
P14(:,1) = sp.P14spal; 
P26(:,1) = sp.P26spal; 

% do we need to calculate Alumimium (currently I run 26Al as a separate,
% assuming that it will be the least measured nuclide, so instead of
% modyfying the current loop, I just add a new loop with an if statement)
% taking Al into the main loop would be cleaner, but I'm trying to optimize
% for speed in the inversion.
Al = any(Nlogical(:,3));

% calculate concentrations --------------------------------------------
N10 = zeros(nSamp, nmix);
N14 = zeros(nSamp, nmix); 

% loop through production pathways (first spallation, then muons)
for i = 1:2 
    
    % segments of production profile
    N10i = zeros(nSamp, nmix);
    N14i = zeros(nSamp,nmix);
    
    for j = 1:length(T)  % loop through erosion segments

        beta10 = rho .* E(:,j) ./ att_l_10(i) + consts.l10;
        beta14 = rho .* E(:,j) ./ att_l_14(i) + consts.l14;
        
        N10_segment = zeros(nSamp, nmix);
        N14_segment = zeros(nSamp, nmix); 

        for d = 1:nmix  % loop through soil mixing depths
            if i==1
                N10_segment(:,d) = P10(:,1)./beta10 .* ...                    % production
                    exp(-(rho .* t_depths_3D(:,j,d) ./ att_l_10(i))) .* ...    % depth cut-off
                    (1 - exp(-beta10.*T_time_spans(j)));              % time to develop the concentration profile
                N14_segment(:,d) = P14(:,1)./beta14 .* ...
                    exp(-(rho .* t_depths_3D(:,j,d) ./ att_l_14(i))) .* ...
                    (1 - exp(-beta14.*T_time_spans(j)));
            else 
                % muon production
                P10(:,2) = intNmu(E(:,j).*rho,sp.pressure,Nmu.pp,Nmu.logee,Nmu.N10quartz); 
                P14(:,2) = intNmu(E(:,j).*rho,sp.pressure,Nmu.pp,Nmu.logee,Nmu.N14quartz); 
    
                if any(isnan(P10(:,2)))   % this is necessary for very fast erosion rates that surpass the pre-calculated rates in Cronus
                    P10(:,2) = zeros(nSamp,1);
                    P14(:,2) = zeros(nSamp,1);
                end
    
                N10_segment(:,d) = P10(:,2) .* ...                           % production, here now beta is needed because Cronus outputs total production
                    exp(-(rho .* t_depths_3D(:,j,d) ./ att_l_10(2))) .* ...    % depth cut-off
                    (1 - exp(-beta10.*T_time_spans(:,j)));              % time to develop the concentration profile
                N14_segment(:,d) = P14(:,2) .* ...
                    exp(-(rho .* t_depths_3D(:,j,d) ./ att_l_14(2))) .* ...
                    (1 - exp(-beta14.*T_time_spans(:,j)));
            end
        end

        N10i = N10i .* exp(-consts.l10.*T_time_spans(j)) + N10_segment; % add segments and don't forget decay
        N14i = N14i .* exp(-consts.l14.*T_time_spans(j)) + N14_segment; % add segments and don't forget decay       
    end

    N10 = N10 + N10i;
    N14 = N14 + N14i;
end
 

if Al       % calculate Al separetely assuming that it will be used the least and we run code faster if we dont calculate this
N26=zeros(nmix,1);
    for i = 1:2 
    N26i = zeros(nSamp, nmix);  
    for j = 1:length(T)  % loop through erosion segments

        beta26 = rho .* E(:,j) ./ att_l_26(i) + consts.l26;
        
        N26_segment = zeros(nSamp, nmix);
        for d = 1:nmix
            if i==1
                N26_segment(:,d) = P26(:,1)./beta26 .* ...                    % production
                    exp(-(rho .* t_depths_3D(:,j,d) ./ att_l_26(i))) .* ...    % depth cut-off
                    (1 - exp(-beta26.*T_time_spans(j)));              % time to develop the concentration profile
            else 
                % muon production
                P26(:,2) = intNmu(E(:,j).*rho,sp.pressure,Nmu.pp,Nmu.logee,Nmu.N26quartz); 
                if any(isnan(P26(:,2)))   % this is necessary for very fast erosion rates that surpass the pre-calculated rates in Cronus
                    P26(:,2) = zeros(nSamp,1);
                end
                N26_segment(:,d) = P26(:,2) .* ...                           % production, here now beta is needed because Cronus outputs total production
                    exp(-(rho .* t_depths_3D(:,j,d) ./ att_l_26(2))) .* ...    % depth cut-off
                    (1 - exp(-beta26.*T_time_spans(:,j)));              % time to develop the concentration profile
            end
        end

        
        N26i = N26i .* exp(-consts.l26.*T_time_spans(j)) + N26_segment; % add segments and don't forget decay
    end
    N26 = N26 + N26i;
    end
end

%% do soil mixing ---------------------------------------------------------
% simply average concentrations acros depth intervals to get concentration
% under soil mixing
N10_mixed = mean(N10,2);
N14_mixed = mean(N14,2);
if Al
    N26_mixed = mean(N26,2);
end
%% output

if Al
    Ncalculated = [N10_mixed;N14_mixed;N26_mixed];
else
    Ncalculated = [N10_mixed;N14_mixed; zeros(nSamp,1)];
end
N = Ncalculated(Nlogical);  % only pass on nuclide values that were measured
end
