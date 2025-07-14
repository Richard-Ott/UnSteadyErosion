function N = Nforward_discretized_FEM(E,T,sp,consts,scenario,Nlogical,varargin)
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
parse(p, varargin{:});
CHG = p.Results.change_variable;

nSamp = length(sp.P10spal);

E = E./1e6;   % convert to m/a
T_time_spans   = [inf, diff(T')*(-1)];        % time span of every time interval between erosion changes

%% reshape erosion rate array and calculate depths at steps for different scenarios

switch scenario
    case 'samestep'
        % make matrix of erosion rates for different time steps by multiplying E
        % with the change factors
        E = [E, repmat(E,1,length(CHG))];
        E(:,2:end) = E(:,2:end) * diag(CHG);
    case 'samebackground_step'
        E = repmat(E,nSamp,length(T));
        E(:,2:end) = E(:,2:end) .* CHG;
    case 'samebackground_samestep'
        E = repmat(E,nSamp,length(T));
        E(:,2:end) = E(:,2:end) * diag(CHG);
    case 'spike'
        E = repmat(E,1, length(T));                      % make matrix of erosion rates for easy calling in concentration loop
    case 'samespike'
        E = repmat(E,1, length(T));
    case 'samebackground_spike'
        E = repmat(E,nSamp,length(T));                      % make matrix of erosion rates for easy calling in concentration loop
    case 'samebackground_samespike'
        E = repmat(E,nSamp,length(T));                      % make matrix of erosion rates for easy calling in concentration loop
    case 'curve'
        scaleFactor = CHG;
        absolute_erosion_change = sp.curvechange .* scaleFactor; % calculate erosion rate change from relative differences from curve and scaling factor
        % make matrix of erosion rates for different time steps by multiplying E with the change factors
        E = repmat(E,1,length(T));
        E(:,2:end) = E(:,2:end)+ E(:,2:end) .* absolute_erosion_change;
end

%% start cosmo calulcation

% Geometry of profile -----------------------------------------------------

D = 25;   % depth of profile (m)
Nn = 100; % number of nodes
z = linspace(0,1,Nn);
z = D*z.^2; % depth node spacing
ze = .5*(z(1:(Nn-1))+z(2:Nn));


TOPO = [1:Nn-1;2:Nn]';  % finite element topogology
dl = diff(z);       

Ne = Nn - 1;            % number of elements

% useful matrices for FEM
M1 = [1/3,1/6;1/6,1/3]; %phi*phi
M2 = [-1,1;-1,1]/2;     %phi*dphi
M3 = [1,-1;-1,1];       %dphi*dphi

% assign variables --------------------------------------------------------
rho = consts.density;  
dmix = sp.mixing*exp(-ze(:)/sp.dzmix);

% decay constants
l10 = consts.l10;
l14 = consts.l14;
l26 = consts.l26;

% attentuation lengths
att_l_10(1) = consts.L_sp;          % spallation
att_l_10(2) = consts.L_mu_atm(4);   % muon 
att_l_14(1) = consts.L_sp;          % spallation
att_l_14(2) = consts.L_mu_atm(5);   % muon 
att_l_26(1) = consts.L_sp;          % spallation
att_l_26(2) = consts.L_mu_atm(7);   % muon 

% spallation surface production rates 
P10(:,1) = sp.P10spal; 
P14(:,1) = sp.P14spal; 
P26(:,1) = sp.P26spal; 

% muon surfce production rates (no erosion)
P10(:,2) = consts.Pmu0(4).* exp((1013.25-sp.pressure)./att_l_10(2)'); 
P14(:,2) = consts.Pmu0(5).* exp((1013.25-sp.pressure)./att_l_14(2)');

% production profiles
P10profile = P10(:,1)*exp(-z(:)*100*rho/att_l_10(1))' + P10(:,2).* exp(-z(:)*100*rho/att_l_10(2))'; % multiply by hundred to convert depth to cm
P14profile = P14(:,1)*exp(-z(:)*100*rho/att_l_14(1))' + P14(:,2).* exp(-z(:)*100*rho/att_l_14(2))';
P10profile = P10profile';  P14profile = P14profile';

%% debugging test compare to Knudsen
E = ones(size(E)).*5e-6;  % debugging
% rho = 1.9;  % debugging test
% P10(:,1) = 4;
% P10(:,2) = P10(:,1).*0.015;
% P10(:,3) = P10(:,1).*0.005;
% 
% att_l_10(1) = 150;          % spallation
% att_l_10(2) = 1500;   % muon 
% att_l_10(3) = 4320;   % muon 
% 
% P10 = P10(:,1)*exp(-z(:)*100*rho/att_l_10(1))' + P10(:,2).* exp(-z(:)*100*rho/att_l_10(2))'; % multiply by hundred to convert depth to cm
% P10 = P10(:);
%% calculate steady state starting profile %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N10 = cell(nSamp,1);
N14 = cell(nSamp,1);

for s = 1:nSamp
    %initiate matrices
    K10 = sparse(Nn,Nn);
    f10 = sparse(Nn,1);
    K14 = sparse(Nn,Nn);
    f14 = sparse(Nn,1);
    
    %loop elements
    for i=1:Ne
    
      %add to stiffness matrix
      Ke10 = l10*M1*dl(i)-E(s,1) *M2+dmix(i)*M3/dl(i);
      K10(TOPO(i,:),TOPO(i,:)) = K10(TOPO(i,:),TOPO(i,:)) + Ke10;
      Ke14 = l14*M1*dl(i)-E(s,1) *M2+dmix(i)*M3/dl(i);
      K14(TOPO(i,:),TOPO(i,:)) = K14(TOPO(i,:),TOPO(i,:)) + Ke14;
      
      %add to Load vector
      fe10 = M1*dl(i)*P10profile(TOPO(i,:)',s);
      fe14 = M1*dl(i)*P14profile(TOPO(i,:)',s);
      
      f10(TOPO(i,:)) = f10(TOPO(i,:)) + fe10;
      f14(TOPO(i,:)) = f14(TOPO(i,:)) + fe14;
    
    end
    
    %Lower boundary condition
    [K10,f10] = ebc(K10,f10,Nn,0);
    [K14,f14] = ebc(K14,f14,Nn,0);
    
    %solve system
    N10{s} = K10\f10;
    N14{s} = K14\f14;
end

%% calculate transient solution concentrations ----------------------------

% loop through erosion history segments
for s = 1:nSamp
    for t = 1:length(T)-1  

        % duration of this erosion period
        dt = T_time_spans(t+1);

        %initiate matrices
        K10 = sparse(Nn,Nn);
        f10 = sparse(Nn,1);
        K14 = sparse(Nn,Nn);
        f14 = sparse(Nn,1);

        % loop elements
        for i=1:Ne

          %add to stiffness matrix
          Ke10 = (1+dt*l10)*M1*dl(i) - E(s,t)*dt*M2+dmix(i) * dt*M3/dl(i);
          K10(TOPO(i,:),TOPO(i,:)) = K10(TOPO(i,:),TOPO(i,:)) + Ke10;
          Ke14 = (1+dt*l14)*M1*dl(i) - E(s,t)*dt*M2+dmix(i)*dt*M3/dl(i);
          K14(TOPO(i,:),TOPO(i,:)) = K14(TOPO(i,:),TOPO(i,:)) + Ke14;

          %add to Load vector
          fe10 = M1*N10{s}(TOPO(i,:)')*dl(i)+dt*M1*dl(i)*P10profile(TOPO(i,:)',s);
          fe14 = M1*N14{s}(TOPO(i,:)')*dl(i)+dt*M1*dl(i)*P14profile(TOPO(i,:)',s);

          f10(TOPO(i,:)) = f10(TOPO(i,:)) + fe10;
          f14(TOPO(i,:)) = f14(TOPO(i,:)) + fe14;

        end

        %Lower boundary condition
        [K10,f10] = ebc(K10,f10,Nn,0);
        [K14,f14] = ebc(K14,f14,Nn,0);

        %update profiles
        N10{s} = full(K10\f10);
        N14{s} = full(K14\f14);
    end
end

% take surface concentration
N10_surface = cellfun(@(x) x(1), N10);
N14_surface = cellfun(@(x) x(1), N14);


% do we need to calculate Alumimium (currently I run 26Al as a separate,
% assuming that it will be the least measured nuclide, so instead of
% modyfying the current loop, I just add a new loop with an if statement)
% taking Al into the main loop would be cleaner, but I'm trying to optimize
% for speed in the inversion.
Al = any(Nlogical(:,3));

if Al       % calculate Al separetely assuming that it will be used the least and we run code faster if we dont calculate this
    N26=zeros(nmix,1);
end


% ADD Al CALUCLATION HERE!!!!!!!!!!!!!!!!

%% output

if Al
    Ncalculated = [N10_surface;N14_surface;N26_surface];
else
    Ncalculated = [N10_surface;N14_surface; zeros(nSamp,1)];
end
N = Ncalculated(Nlogical);  % only pass on nuclide values that were measured
end
