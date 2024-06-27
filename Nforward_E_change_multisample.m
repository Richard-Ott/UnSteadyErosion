function N = Nforward_E_change_multisample(E,T,CHG,sp,consts,Nmu)
% This function calculates concentrations N10 and N14 for mutiple step 
% changes in erosion and for multiple samples at once.
%
% Input:
%       - E: n x 1 vector with erosion rates in mm/ka for n samples
%       - T: 1 x m vector of erosion rate step change timings (yrs BP)
%       from old to young
%       - CHG: 1 x m change factors of erosion rates at times T
%       - sp: sample specific parameters (Pspal & pressure)
%       - consts: the consts_v3.mat file from Cronus v3
%       - precalculated muon production rates from Cronus v3 (Nmu.mat)
% Output:
%       - [N10,N14]: concentrations at/g of 10Be and 14C
%
% Richard Ott & Dirk Scherler, 2024

nSamp = length(sp.P10spal);

% make matrix of erosion rates for different time steps by multiplying E
% with the change factors
E = [E, repmat(E,1,length(CHG))];
E(:,2:end) = E(:,2:end) * diag(CHG);

E = E./1e4;   % convert to cm/a

T_time_spans   = [inf, diff(T')*(-1)];        % time span of every time interval between erosion changes
segment_depths = E.*T_time_spans;            % the exhumation occuring during every erosion time interval in cm
t_depths       = [fliplr(cumsum(fliplr(segment_depths(:,2:end)),2)), zeros(nSamp,1)];     % depths at every step change T in cm

% assign variables --------------------------------------------------------
rho = consts.density;  

% attentuation lengths
att_l_10(1) = consts.L_sp;          % spallation
att_l_10(2) = consts.L_mu_atm(4);   % muon 
att_l_14(1) = consts.L_sp;          % spallation
att_l_14(2) = consts.L_mu_atm(5);   % muon 

% production rates 
P10(:,1) = sp.P10spal; 
P14(:,1) = sp.P14spal; 


% calculate concentrations --------------------------------------------
N10 = 0;
N14 = 0; 

% loop through production pathways (first spallation, then muons)
for i = 1:2 
    
    % segments of production profile
    N10i = 0;
    N14i = 0;
    
    for j = 1:size(E,2)  % loop through erosion segments

        beta10 = rho .* E(:,j) ./ att_l_10(i) + consts.l10;
        beta14 = rho .* E(:,j) ./ att_l_14(i) + consts.l14;
    
        if i==1
            N10_segment = P10(:,1)./beta10 .* ...                    % production
                exp(-(rho .* t_depths(:,j) ./ att_l_10(i))) .* ...    % depth cut-off
                (1 - exp(-beta10.*T_time_spans(j)));              % time to develop the concentration profile
            N14_segment = P14(1)./beta14 .* ...
                exp(-(rho .* t_depths(:,j) ./ att_l_14(i))) .* ...
                (1 - exp(-beta14.*T_time_spans(j)));
        else 
            % muon production
            P10(:,2) = intNmu(E(:,j).*rho,sp.pressure',Nmu.pp,Nmu.logee,Nmu.N10quartz); 
            P14(:,2) = intNmu(E(:,j).*rho,sp.pressure',Nmu.pp,Nmu.logee,Nmu.N14quartz); 
            
            if any(isnan(P10(:,2)))   % this is necessary for very fast erosion rates that surpass the pre-calculated rates in Cronus
                P10(:,2) = zeros(nSamp,1);
                P14(:,2) = zeros(nSamp,1);
            end

            N10_segment = P10(:,2) .* ...                           % production, here now beta is needed because Cronus outputs total production
                exp(-(rho .* t_depths(:,j) ./ att_l_10(2))) .* ...    % depth cut-off
                (1 - exp(-beta10.*T_time_spans(:,j)));              % time to develop the concentration profile
            N14_segment = P14(:,2) .* ...
                exp(-(rho .* t_depths(:,j) ./ att_l_14(2))) .* ...
                (1 - exp(-beta14.*T_time_spans(:,j)));
        end

        N10i = N10i .* exp(-consts.l10.*T_time_spans(j)) + N10_segment; % add segments and don't forget decay
        N14i = N14i .* exp(-consts.l14.*T_time_spans(j)) + N14_segment; % add segments and don't forget decay       
    end

    N10 = N10 + N10i;
    N14 = N14 + N14i;
end
 

N = [N10;N14];
end
