function out = get_erates_v3(in,control,consts)

% This computes a block of erosion rates from a
% validated data structure. 
% 
% input args in,control,consts
%
% in is a structure that is the output of validate_v3 or validate_v2. 
% this has s,n,c sub-structures for sample, nuclide, calibration data. 
% 
% control has .resultType = 'short' or 'long' to control NTD or NTD + TD
% If not supplied, default is 'short.' 
%
% consts is consts structure from make_consts_v3; optional
%
% Returns output arg out which is just the input arg in with a lot of 
% new fields added to it. 
%
% Greg Balco
% Berkeley Geochronology Center
% March, 2017
% 
% Development. Not licensed for use or distribution. 

% Should pass consts to minimize file loads, but it is not a very big file
% so, whatever. Passing consts is needed if one wants to use non-default
% production rates or other parameters. 

if nargin < 3;
    load consts_v3;
end;

% Also load precalculated muon data
% load([consts.datadir 'Nmu.mat']);
load('C:\Users\rott\OneDrive - UvA\Richard\PhD_ETH\matlab\CRONUS cosmo calculation\cronus 3.0\online-calculators-v3\data\Nmu.mat')

% Defaults
if nargin < 2;
    control.resultType = 'short'; % Just compute NTD results
end;

% Version
in.version_erates = '3.0';
in.version_muons = Nmu.version;

% Initialize

in.flags = [];

% ---------------- COMPUTE SAMPLE-SPECIFIC PARAMETERS ------------------

% 1. Compute sample-specific parameters: pressure, various SFs

% Calculate pressures with minimal number of calls to ERA40atm.m.
% In this section is where we would pull an elevation history from an 
% isostatic rebound model. 

calc_ps_ERA40 = find(strcmp(in.s.aa,'std')); % Find indices to calculate
calc_ps_ant = find(strcmp(in.s.aa,'ant'));

% Main point here is to call ERA40.m only once
in.s.pressure(calc_ps_ERA40) = ERA40atm(in.s.lat(calc_ps_ERA40),in.s.long(calc_ps_ERA40),in.s.elv(calc_ps_ERA40));
in.s.pressure(calc_ps_ant) = antatm(in.s.elv(calc_ps_ant));

% Note if atmospheric pressure were made time-dependent, then one would 
% need a cell array rather than a vector, which will incur lots of down-
% stream code changes. 

% ------------ COMPUTE NUCLIDE-SPECIFIC PARAMETERS  --------------------

% Get L; this is placeholder for now but may be spatially variable, so it 
% is a vector. It belongs to n not s so that it can also vary by nuclide if
% needed. 

in.n.Lsp = zeros(size(in.n.index)) + Lsp();

% Thickness correction; by extension this is also a property of n
in.n.sf_thick = (thickness(in.s.thick(in.n.index),in.n.Lsp,in.s.rho(in.n.index)))';    
% Transpose is because thickness.m returns wrong shape

% St SF. Then on the other hand this one is technically a nuclide-
% independent property of the sample, but we put it in the nuclide 
% properties for convenience. 

in.n.sf_St = (stone2000(in.s.lat(in.n.index),in.s.pressure(in.n.index),1))'; 
% Note transpose because stone2000.m returns wrong shape.

% Need to obtain index identifying what nuclide each entry in n is, for 
% use in extracting nuclide-specific properties from the consts structure. 
% Really this should probably be done upstream, but whatever. 
in.n.nindex = zeros(size(in.n.index));
for a = 1:length(in.n.index);
    % Loop over all nuclide measurements
    % Determine what index to use in consts file
    in.n.nindex(a) = find(strcmp(in.n.nuclide{a},consts.nuclides));
    % This isn't robust against errors in the nuclide specification, but 
    % shouldn't need to be if validate_v3_input is doing its job correctly. 
end;

% Decay constants

in.n.l = consts.l(in.n.nindex)'; % Transpose needed for some reason.
in.n.dell = consts.dell(in.n.nindex)'; % Likewise. 


% ----------- COMPUTE NTD EROSION RATES ----------------------------------

% 3. If unknowns rather than calibration data, compute St exposure ages.

% Total production rate on St
in.n.Psp_St = in.n.sf_St.*in.s.othercorr(in.n.index).*in.n.sf_thick.*consts.refP_St(in.n.nindex)';
in.n.delPsp_St = in.n.sf_St.*in.s.othercorr(in.n.index).*in.n.sf_thick.*consts.delrefP_St(in.n.nindex)';

% Normalized nuclide concentrations (for plotting).
% This requires the surface production rate due to muons. Use fast
% approximation. That's also needed for saturation checks. 
in.n.P_mu = consts.Pmu0(in.n.nindex)'.* exp((1013.25-in.s.pressure(in.n.index))./consts.L_mu_atm(in.n.nindex)');
in.n.Ptot_St = in.n.Psp_St + in.n.P_mu; % No thick/topo scaling for muon production
in.n.Nnorm_St = in.n.N./in.n.Ptot_St; 
in.n.delNnorm_St = in.n.delN./in.n.Ptot_St;
% Also include production rate uncert for potential use in plotting
in.n.delNnorm_ext_St = sqrt((in.n.delN./in.n.Ptot_St).^2 + (in.n.delPsp_St.*in.n.N./(in.n.Ptot_St.^2)).^2);

% Non-time-dependent age calculations (St)

% Allocate
in.n.E_St = zeros(size(in.n.index));
in.n.delE_ext_St = in.n.E_St; % Used much later
in.n.delE_int_St = in.n.E_St; % Likewise

% Set opts. Not sure what correct tolerance is for Octave. Work on this. 
opts = optimset('TolX',1e-4);

% Compute the erosion rates

for a = 1:length(in.n.N);        
    % Saturation check
    if in.n.N(a) >= in.n.Ptot_St(a)./in.n.l(a);
        % Case saturated
        in.n.E_St(a) = 0;
        flag = ['Sample ' char(in.s.sample_name{in.n.index(a)}) ' -- ' char(in.n.nuclide{a}) ' appears to be saturated WRT Stone(2000).'];
        in.flags = [in.flags ' <br> ' sprintf('%s\n','') flag]; 
    else
        % Case below saturation; can calculate erosion rate
        % rootfinding scheme with precalculated muon inventories
        % compute limits
        % upper bound on erosion rate is given by long attenuation length
        xmax = (4000./in.n.N(a)).*(in.n.Ptot_St(a) - in.n.N(a).*in.n.l(a));
        % lower bound on erosion rate is given by short attenuation length
        xmin = (100./in.n.N(a)).*(in.n.Ptot_St(a) - in.n.N(a).*in.n.l(a));
        
        % Define which thing to use
        eval(['this_Nmu = Nmu.' in.n.nuclide{a} ';'])
        
        % Do rootfinding
        [this_x,this_fval,this_exitflag] = fzero(@(x) (in.n.N(a) - (in.n.Psp_St(a)./(in.n.l(a) + x./in.n.Lsp(a))) - intNmu(x,in.s.pressure(in.n.index(a)),Nmu.pp,Nmu.logee,this_Nmu)),[xmin xmax]);
        if this_exitflag ~= 1;
            % Rootfinder failed
            in.n.E_St(a) = 0;
            flag = ['Sample ' char(in.s.sample_name{in.n.index(a)}) ' -- ' char(in.n.nuclide{a}) ' fzero failed for St scaling.'];
            in.flags = [in.flags ' <br> ' sprintf('%s\n','') flag];
        end;  
        % Assign
        in.n.E_St(a) = this_x;       
    end;
end;
    
% ---------- DONE WITH NTD EROSION RATES ----------------------

% ---------- NOW GENERATE TIME-DEPENDENT SF ------------------------------

% Default case only NTD scaling
sfa = {'St'};

% 5. If long-form results called for, 
if strcmp(control.resultType,'long');
    % Case NTD scaling schemes
    sfa = {'St' 'Lm' 'LSDn'}; % Add here when adding code
    
    % 5A. Generate cutoff rigidity time series for all samples. 
    % Assemble input data
    rin.lat = in.s.lat(in.n.index);
    rin.long = in.s.long(in.n.index);
    rin.yr = in.s.yr(in.n.index);
    rin.t = zeros(size(rin.yr)) + 2005000; % Set to end of pmag record
    % because, of course, production rate is not time varying before that.
    % This is a bit of a waste of time for C-14 data, but whatever. 
    
    % Get Rc vector. This builds a new structure called sfdata. 
    % Currently using dipole approximation in get_DipRc.m. 
      
    sfdata = get_DipRc(rin,consts.RCfname,consts.Sfname);
   
    % Now the structure sfdata contains a bunch of Rc vectors, one for each
    % analysis. To compute scaling factors, we'll also need to know sample
    % elevations as well as which nuclide we are dealing with. 
    
    sfdata.nuclide = in.n.nuclide;
    sfdata.pressure = in.s.pressure(in.n.index);
    
    % Now sfdata is a structure with:
    %   cell array sfdata.t contains time vectors
    %   cell array sfdata.Rc contains Rc at corresponding times
    %   cell array sfdata.S contains Lifton solar parameter 
    %   cell array nuclide is nuclide identifier string
    %   vector pressure is atmospheric pressure. Potentially this could be
    %   a cell array if there is an isostatic rebound model upstream. 
       
    % 5B. Get SF time series from precalculated results.
    % Relies on get_LmSF.m, get_LSDnSF.m
    for a = 1:length(sfa);
        if ~strcmp(sfa{a},'St');
            eval(['sfdata = get_' sfa{a} 'SF(sfdata,consts.datadir);'])
        end;
    end;   
    
    % sfdata now has:
    %   cell array LSDn
    %   cell array Lm
    % etc. These have vectors of "scaling factor". 
    
    % Add sfdata to output structure so it is available later if needed

    in.f = sfdata;
    
    % 5C. Now calculate erosion rates using scaling data. 
    for a = 2:length(sfa); % Assume sfa(1) = 'St'
        % Initialize
        eval(['in.n.E_' sfa{a} '= zeros(size(in.n.N));']);
        eval(['in.n.delE_int_' sfa{a} ' = zeros(size(in.n.N));']);
        eval(['in.n.delE_ext_' sfa{a} '= zeros(size(in.n.N));']);
        eval(['in.n.max_N_' sfa{a} ' = zeros(size(in.n.N));']);
        for b = 1:length(in.n.N);
            % Generate P0 with site scaling adjustments
            % This is a scalar and includes thickness/shielding corrections. 
            eval(['thisP0 = consts.refP_' sfa{a} '(in.n.nindex(b)).*in.s.othercorr(in.n.index(b)).*in.n.sf_thick(b);']);

            % Obtain spallogenic surface production rate in each time step. 
            % This is a vector. 
            eval(['pp = thisP0.*sfdata.' sfa{a} '{b};']);

            % Define which thing to use
            eval(['this_Nmu = Nmu.' in.n.nuclide{b} ';'])
            
            % Do a saturation check
            % Spallation integral
            if in.n.l(b) > 0;
                max_N = sum(-(pp./in.n.l(b)).*(exp(-sfdata.tmax{b}.*in.n.l(b)) - exp(-sfdata.tmin{b}.*in.n.l(b)))) ...
                    + pp(end).*exp(-in.n.l(b).*rin.t(b))./(in.n.l(b)) + in.n.P_mu(b)./in.n.l(b);
                eval(['in.n.max_N_' sfa{a} '(b) = max_N;']);
            else
                max_N = Inf; % Can't saturate for stable nuclide
            end;
            
            if in.n.N(b) > max_N;
               % Case saturated
               eval(['in.n.E_' sfa{a} '(b) = 0;']);
               flag = ['Sample ' char(in.s.sample_name(in.n.index(b))) ' -- ' char(in.n.nuclide{b}) ' appears to be saturated WRT ' char(sfa{a}) ' scaling.'];
               in.flags = [in.flags ' <br> ' sprintf('%s\n','') flag]; 
            else  
                % Define function of which to find zeros
                % This has the measured nuclide concentration, less:
                % 1. the muon integral, 
                % 2. the steady-state concentration at the beginning of
                % the time-dependent part, decayed until the present time,
                % and 3. the concentration developed during
                % the time-dependent part. 
                % (2) is probably very small. 
                efun = @(x) (in.n.N(b) - intNmu(x,in.s.pressure(in.n.index(b)),Nmu.pp,Nmu.logee,this_Nmu) ...
                    - exp(-rin.t(b).*in.n.l(b)).*pp(end).*exp(-rin.t(b).*x./in.n.Lsp(b))./(in.n.l(b)+x/in.n.Lsp(b)) ...
                    - sum(-(pp./(in.n.l(b)+x/in.n.Lsp(b))).*(exp(-sfdata.tmax{b}.*(in.n.l(b)+x/in.n.Lsp(b))) - exp(-sfdata.tmin{b}.*(in.n.l(b)+x/in.n.Lsp(b))))) );

                % compute limits
                % upper bound on erosion rate is given by long attenuation length
                % and max production rate
                xmax = (4000./in.n.N(b)).*(max(pp) - in.n.N(b).*in.n.l(b));
                % lower bound on erosion rate is given by short attenuation length
                % and min production rate
                xmin = (100./in.n.N(b)).*(min(pp) - in.n.N(b).*in.n.l(b));
                
                % Catch pathological situations near saturation
                if xmin < 0; xmin = 0; end;
                if xmax < 0; xmax = 1e-2; end; % Unclear what this should be. Making it too wide probably just slows things down. 

                % Do rootfinding
                try
                    [this_x,this_fval,this_exitflag] = fzero(@(x) efun(x),[xmin xmax]);
                catch
                    % Rootfinder failed
                    eval(['in.n.E_' sfa{a} '(b) = 0;']);
                    flag = ['Sample ' char(in.s.sample_name{in.n.index(a)}) ' -- ' char(in.n.nuclide{a}) ' fzero crashed for ' sfa{a} ' scaling.'];
                    in.flags = [in.flags ' <br> ' sprintf('%s\n','') flag];
                end;
                if this_exitflag ~= 1;
                    % Rootfinder failed
                    eval(['in.n.E_' sfa{a} '(b) = 0;']);
                    flag = ['Sample ' char(in.s.sample_name{in.n.index(a)}) ' -- ' char(in.n.nuclide{a}) ' fzero returned failure for ' sfa{a} ' scaling.'];
                    in.flags = [in.flags ' <br> ' sprintf('%s\n','') flag];
                end;  
                % Assign
                eval(['in.n.E_' sfa{a} '(b) = this_x;']);   
            end;
        end;
    end;   
end; % End loop executed when result type is 'long'

% --------- END ALL TIME-DEPENDENT EROSION RATES -------------
    
% 6. Compute linearized uncertainties. 

for b = 1:length(sfa);
        sfs = sfa{b}; % String for scaling
        % Extract fractional production rate uncertainty for current scaling
        eval(['fdelP = (consts.delrefP_' sfs '(in.n.nindex)./consts.refP_' sfs '(in.n.nindex))'';']);
        
        % extract erosion rate for current scaling
        eval(['this_ee = in.n.E_' sfs ';']); 
        
        % Flag saturated
        ok = find(this_ee > 0); notok = find(this_ee == 0); % flag saturated ones
        
        % Allocate
        FP = zeros(size(in.n.index));
        delFP = FP; dedN = FP; dedP = FP; dedl = FP;
        
        % Compute effective production rate
        % This works for both stable and radionuclides.  
        % This is fine for error propagation, but it's pathological for
        % computing the normalized nuclide concentrations. That aspect
        % needs work and is not implemented at present. 
        FP = in.n.N.*(in.n.l + this_ee./in.n.Lsp);
            
        % Continue with error propagation. 
        % If saturated, no error propagation. 
        
        delFP(ok) = fdelP(ok) .* FP(ok);
        
        % Partials
        dedN(ok) = -in.n.Lsp(ok).*FP(ok)./(in.n.N(ok).^2);
        dedP(ok) = in.n.Lsp(ok)./in.n.N(ok);
        dedl(ok) = -in.n.Lsp(ok);

        % make respective delt's
        eval(['in.n.delE_ext_' sfa{b} '(ok) = sqrt((dedN(ok).*in.n.delN(ok)).^2 + (dedP(ok).*delFP(ok)).^2 + (dedl(ok).*in.n.dell(ok)).^2);']);
        eval(['in.n.delE_int_' sfa{b} '(ok) = sqrt((dedN(ok).* in.n.delN(ok)).^2);']);
        
        clear ttt fdelP ok case1 case2 FP delFP dedN dedP temp_sat_FP
    end;

        
% ---------- DONE WITH ERROR PROPAGATION ----------------------

% 6. Return output data, added to input data

out = in;