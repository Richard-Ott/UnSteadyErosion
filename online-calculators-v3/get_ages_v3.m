function out = get_ages_v3(in,control,consts)

% This computes a block of exposure ages or production rates from a
% validated data structure. Overall aim is to compute all exposure ages as
% quickly as possible whilst minimizing file loads. 
% 
% input args in,control,consts
%
% in is a structure that is the output of validate_v3 or validate_v2. 
% this has s,n,c sub-structures for sample, nuclide, calibration data. 
% 
% control has .cal (0 = compute exposure ages; 1 = compute production
% rates);
% and .resultType = 'short' or 'long' to control NTD or NTD + TD
% If not supplied, defaults are 0 and 'short.' 
%
% consts is consts structure from make_consts_v3; optional
%
% Returns output arg out which is just the input arg in with a lot of 
% new fields added to it. 
%
% This will generate production rates for any nuclide measurements that are
% supplied in the input structure. So if you want to work with only one
% nuclide, others should be stripped out of the input structure upstream. 
%
% Greg Balco
% Berkeley Geochronology Center
% February, 2016
% 
% November 2016: Switched to simple muon approximation (3.0.2). 
% Also switched to piecewise-constant production rate for time-dependent
% SF. Checked entire script pretty thoroughly and fixed several bugs in the
% error propagation. 


% Should pass consts to minimize file loads, but it is not a very big file
% so, whatever. Passing consts is needed if one wants to use non-default
% production rates or other parameters. 

if nargin < 3;
    load consts_v3;
end;

% Defaults
if nargin < 2;
    control.cal = 0; % Not a calibration data set, calculate ages
    control.resultType = 'short'; % Just compute NTD results
end;

% Version
in.version_ages = '3.0.2';

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

% Muon production rate calculations use fast approximation scheme. 
% Variable by nuclide; belongs to n

in.n.P_mu = consts.Pmu0(in.n.nindex)'.* exp((1013.25-in.s.pressure(in.n.index))./consts.L_mu_atm(in.n.nindex)');
in.n.Lmu = 1./(consts.Lmu.a + consts.Lmu.b.*in.s.pressure(in.n.index));
in.version_muons = consts.version_muons;

% Calculate A's for simplicity later on. A = l + E/L. 
% Note exponential approximation for muon production in this code, so need
% A_sp and A_mu
% Note if stable nuclide, no erosion, these evaluate to 0
in.n.A_sp = (in.n.l + in.s.E(in.n.index).*in.s.rho(in.n.index)./in.n.Lsp);
in.n.A_mu = (in.n.l + in.s.E(in.n.index).*in.s.rho(in.n.index)./in.n.Lmu); 
% Note A_mu will not be zero for He-3, but it should always be multiplied
% by zero later. 


% ----------- COMPUTE NTD EXPOSURE AGES ----------------------------------

if control.cal == 0;
    % 3. If unknowns rather than calibration data, compute St exposure ages.
    % This is mostly vectorizable. 
        
    % Total production rate on St
    in.n.Psp_St = in.n.sf_St.*in.s.othercorr(in.n.index).*in.n.sf_thick.*consts.refP_St(in.n.nindex)';
    in.n.delPsp_St = in.n.sf_St.*in.s.othercorr(in.n.index).*in.n.sf_thick.*consts.delrefP_St(in.n.nindex)';
    in.n.P_St = in.n.Psp_St + in.n.P_mu; % No thick/topo scaling for muon production
    
    % Normalized nuclide concentrations (for plotting); 
    in.n.Nnorm_St = in.n.N./in.n.P_St; 
    in.n.delNnorm_St = in.n.delN./in.n.P_St;
    % Also include production rate uncert for potential use in plotting
    in.n.delNnorm_ext_St = sqrt((in.n.delN./in.n.P_St).^2 + (in.n.delPsp_St.*in.n.N./(in.n.P_St.^2)).^2);
    
    % Non-time-dependent age calculations (St)
    
    % Allocate
    in.n.t_St = zeros(size(in.n.index));
    in.n.delt_ext_St = in.n.t_St; % Used much later
    in.n.delt_int_St = in.n.t_St; % Likewise
    
    % Set opts. Not sure what correct tolerance is for Octave. Work on this. 
   	opts = optimset('TolX',1e-4);
    
    % Two cases: A = 0 (stable nuclide, no erosion) or not. 
    is_simple = (in.n.A_sp == 0);
    case1 = find(is_simple); case2 = find(~is_simple);
    
    % Case 1 is vectorizable; also saturation is impossible
    in.n.t_St(case1) = in.n.N(case1)./in.n.P_St(case1);
    
    % Case 2 requires looping and potentially use of fzero 
    for a = 1:length(case2);        
        % Case radioactive decay, erosion, or both
        % Two-term simple age equation
        Pv = [in.n.Psp_St(case2(a)) in.n.P_mu(case2(a))];
        Av = [in.n.A_sp(case2(a)) in.n.A_mu(case2(a))];
        % Need to check for saturation first
        if in.n.N(case2(a)) >= sum(Pv./Av);
            % Case saturated
            in.n.t_St(case2(a)) = 0;
            flag = ['Sample ' char(in.s.sample_name{in.n.index(case2(a))}) ' -- ' char(in.n.nuclide{case2(a)}) ' appears to be saturated WRT Stone(2000) SF at this erosion rate.'];
            in.flags = [in.flags ' <br> ' sprintf('%s\n','') flag]; 
        else
            % Case below saturation; can calculate age
            if in.s.E(in.n.index(case2(a))) == 0;
                % Case zero erosion; analytical solution feasible because 
                % A_sp = A_mu = lambda
                in.n.t_St(case2(a)) = (-1./in.n.A_sp(case2(a))).*log(1-(in.n.N(case2(a)).*in.n.A_sp(case2(a))./(in.n.P_St(case2(a)))));
            else
                % Case nonzero erosion; implicit solution 
                % Need bounds for fzero. Basically, if we make the attenuation
                % length longer, then the surface nuclide concentration is 
                % bigger for the same age. Thus, for given nuclide concentration, 
                % including muon production, which increases the attenuation length,
                % results in a lower age. So the simple age with L = Lsp is 
                % an upper bound on
                % the actual age, and the lower bound is the simple age 
                % calculated with L = Lmu.
                % Figure upper bound by one-term S.A.E., but also have to do
                % saturation check
                if (in.n.N(case2(a)) >= ((in.n.P_St(case2(a)))./in.n.A_sp(case2(a))));
                    % if above nominal saturation, simple age equation
                    % would evaluate to NaN. Avoid this; set upper bound
                    % Inf. Note this may cause trouble. 
                    ub = Inf;
                else
                    % Not saturated for 1-term S.A.E., calculate it
                    ub = (-1./in.n.A_sp(case2(a))).*log(1-(in.n.N(case2(a)).*in.n.A_sp(case2(a))./(in.n.P_St(case2(a)))));
                end;
                % Figure lower bound with S.A.E. 
                lb = (-1./in.n.A_mu(case2(a))).*log(1-(in.n.N(case2(a)).*in.n.A_mu(case2(a))./(in.n.P_St(case2(a)))));
                % Now use fzero
                % two-term S.A.E. 
                fun = @(t) (in.n.N(case2(a)) - sum((Pv./Av).*(1 - exp(-t.*Av))));
                % Solve
                in.n.t_St(case2(a)) = fzero(fun,[lb.*0.99 ub.*1.01],opts);
            end;
        end;
    end;
    
    % At this point, we have computed all the NTD exposure ages, but not
    % yet the uncertainties. 
        
% ---------------- COMPUTE NTD PRODUCTION RATES --------------------------
    
elseif control.cal == 1;
    % 4. If calibration site, compute St production rate. 
    % For saturated sites, e.g. for C-14, it should not be necessary to
    % treat them separately if age is specified appropriately. That is,
    % one should ehter saturated sites with truet old enough to ensure
    % saturation, i.e. 8 or 10x the half-life.
    
    % Allocate a lot of stuff for max/min/true age constraints
    in.n.calc_P_St = zeros(size(in.n.index));
    in.n.calc_delP_St = in.n.calc_P_St;
    in.n.calc_minP_St = in.n.calc_P_St;
    in.n.calc_delminP_St = in.n.calc_P_St;
    in.n.calc_maxP_St = in.n.calc_P_St;
    in.n.calc_delmaxP_St = in.n.calc_P_St;
    in.n.Nmu = in.n.calc_P_St;
    in.n.Nsp = in.n.Nmu;
    in.n.max_Nmu = in.n.Nmu;
    in.n.max_Nsp = in.n.Nsp;
    in.n.min_Nmu = in.n.Nmu;
    in.n.min_Nsp = in.n.Nsp;
    % Determine when min/max calcs are needed
    istruet = find(in.c.truet(in.n.index) > 0);
    ismint = find(in.c.mint(in.n.index) > 0);
    ismaxt = find(in.c.maxt(in.n.index) < Inf);
    
    % Correct independent ages (supplied in years before 1950) to years before
    % sample collection. This gives a true exposure age. 
    % Note that structures c and s have same number of elements in each
    % field. So,
    % (age before sample collection) = (age before 1950 + (date of sample
    % collection - 1950)). 
    % Note this permanently changes the data in c. 
    in.c.truet(in.c.truet > 0) = in.c.truet(in.c.truet > 0) + in.s.yr(in.c.truet > 0) - 1950;
    in.c.mint(in.c.mint > 0) = in.c.mint(in.c.mint > 0) + in.s.yr(in.c.mint > 0) - 1950;
    in.c.maxt(in.c.maxt < Inf) = in.c.maxt(in.c.maxt < Inf) + in.s.yr(in.c.maxt < Inf) - 1950;   
    
    % Two cases: A = 0 (stable nuclide, no erosion) or not. 
    is_simple = (in.n.A_sp == 0);
    
    % Break up to deal with exact-max-min cases, stay vectorized. 
    % This probably could be shortened with main equations moved to a
    % subroutine, but it is not that much code. Do be careful to avoid
    % errors in redundant equations. 
    
    if ~isempty(istruet);
        % First, sites with exact age
        ok = (in.c.truet(in.n.index) > 0); 
        % Case 1 (A = 0)/ Case 2
        ok1 = ok & is_simple; ok2 = ok & (~is_simple);
        % Estimate muon-produced inventory
        % Case A = 0
        in.n.Nmu(ok1) = in.n.P_mu(ok1).*in.c.truet(in.n.index(ok1));
        % Case A ~= 0
        in.n.Nmu(ok2) = (in.n.P_mu(ok2)./in.n.A_mu(ok2)).*(1-exp(-in.n.A_mu(ok2).*in.c.truet(in.n.index(ok2)))); 
        % Remove it; 
        in.n.Nsp(ok) = (in.n.N(ok) - in.n.Nmu(ok));
        % Compute implied sample-specific Psp
        % Case A = 0, Psp just is Nsp/t
        in.n.calc_P_St(ok1) = in.n.Nsp(ok1)./(in.c.truet(in.n.index(ok1)));
        % Case A ~= 0, Psp is Nsp*Asp./(1-exp(-Asp*t))
        in.n.calc_P_St(ok2) = (in.n.Nsp(ok2).*in.n.A_sp(ok2))./(1 - exp(-in.n.A_sp(ok2).*in.c.truet(in.n.index(ok2))));
        % Finally, correct both cases to a reference production rate
        in.n.calc_P_St(ok) = in.n.calc_P_St(ok)./(in.n.sf_St(ok).*in.s.othercorr(in.n.index(ok)).*in.n.sf_thick(ok));
    end;
    
    if ~isempty(ismint);
        % Second, sites with minumum age
        % Of course, this yields a maximum production rate, so in var
        % structure mint leads to maxp. 
        ok = (in.c.mint(in.n.index) > 0);
        ok1 = ok & is_simple; ok2 = ok & (~is_simple);
        % Estimate muon-produced inventory
        % Case A = 0
        in.n.min_Nmu(ok1) = in.n.P_mu(ok1).*in.c.mint(in.n.index(ok1));
        % Case A ~= 0
        in.n.min_Nmu(ok2) = (in.n.P_mu(ok2)./in.n.A_mu(ok2)).*(1-exp(-in.n.A_mu(ok2).*in.c.mint(in.n.index(ok2))));
        % Remove it; 
        in.n.min_Nsp(ok) = (in.n.N(ok) - in.n.min_Nmu(ok));
        % Compute implied sample-specific Psp
        % Case A = 0, Psp just is Nsp/t
        in.n.calc_maxP_St(ok1) = in.n.min_Nsp(ok1)./(in.c.mint(in.n.index(ok1)));
        % Case A ~= 0
        in.n.calc_maxP_St(ok2) = (in.n.min_Nsp(ok2).*in.n.A_sp(ok2))./(1 - exp(-in.n.A_sp(ok2).*in.c.mint(in.n.index(ok2))));
        % Finally, correct both cases to a reference production rate
        in.n.calc_maxP_St(ok) = in.n.calc_maxP_St(ok)./(in.n.sf_St(ok).*in.s.othercorr(in.n.index(ok)).*in.n.sf_thick(ok));
    end;

    if ~isempty(ismaxt);
        % Third, sites with maximum age
        % Of course, this yields a minimum production rate, so in var 
        % structure maxt leads to minp.
        ok = (in.c.maxt(in.n.index) < Inf); % Exists a max age
        ok1 = ok & is_simple;
        ok2 = ok & (~is_simple);
        % Estimate muon-produced inventory
        % Case A = 0
        in.n.max_Nmu(ok1) = in.n.P_mu(ok1).*in.c.maxt(in.n.index(ok1));
        % Case A ~= 0
        in.n.max_Nmu(ok2) = (in.n.P_mu(ok2)./in.n.A_mu(ok2)).*(1-exp(-in.n.A_mu(ok2).*in.c.maxt(in.n.index(ok2))));
        % Remove it; 
        in.n.max_Nsp(ok) = (in.n.N(ok) - in.n.max_Nmu(ok));
        % Compute implied sample-specific Psp
        % Case A = 0, Psp just is Nsp/t
        in.n.calc_minP_St(ok1) = in.n.max_Nsp(ok1)./(in.c.maxt(in.n.index(ok1)));
        % Case A ~= 0
        in.n.calc_minP_St(ok2) = (in.n.max_Nsp(ok2).*in.n.A_sp(ok2))./(1 - exp(-in.n.A_sp(ok2).*in.c.maxt(in.n.index(ok2))));
        % Finally, correct both cases to a reference production rate
        in.n.calc_minP_St(ok) = in.n.calc_minP_St(ok)./(in.n.sf_St(ok).*in.s.othercorr(in.n.index(ok)).*in.n.sf_thick(ok));
    end;
       
end;

% ---------- DONE WITH NTD AGES OR PRODUCTION RATES ----------------------

% ---------- NOW GENERATE TIME-DEPENDENT SF ------------------------------

% Default case only NTD scaling
sfa = {'St'};

% 5. If long-form results called for, 
if strcmp(control.resultType,'long');
    % Case NTD scaling schemes
    sfa = {'St' 'Lm' 'LSDn'}; % Add here when adding code
    
    % 5A. Generate cutoff rigidity time series for all samples. 
    % Note that at this stage we know either the St age or the true 
    % calibration age, so we can pull only the length needed. 
    % Assemble input data
    rin.lat = in.s.lat(in.n.index);
    rin.long = in.s.long(in.n.index);
    rin.yr = in.s.yr(in.n.index);
    if control.cal == 0;
        % Case exposure age; add buffer
        rin.t = in.n.t_St.*1.5; % That should be enough. 
        % If in.n.t_St is zero, then was saturated. To try again to see if
        % saturated using other SF, use time limit defined by 8x the
        % effective half-life (the effective half-life is -ln(0.5)/Asp).  
        is_sat_St = find(in.n.t_St == 0);
        rin.t(is_sat_St) = -8.*log(0.5)./in.n.A_sp(is_sat_St);
        if any (rin.t(is_sat_St) == 0);
            error('get_ages_v3.m: saturated WRT St but A = 0?');
        end;
        % Also, for carbon-14, standardize at 100 ka to make sure nothing
        % goes awry for near-saturated situations. 
        is_C14 = find(in.n.nindex == 5);
        rin.t(is_C14) = 1e5 + zeros(size(is_C14));
    else
        % Case calibration. Idea is to do calculation for either truet (if
        % known) or longest of maxt and mint, depending on which are known.
        rin.t = in.c.truet(in.n.index); % Normal case
        % If maxt exists (i.e., not set to Inf), use that
        usemaxt = ~isinf(in.c.maxt(in.n.index));
        rin.t(usemaxt) = in.c.maxt(in.n.index(usemaxt));
        % If maxt doesn't exist but mint does, use mint
        usemint = isinf(in.c.maxt(in.n.index)) & (in.c.mint(in.n.index) > 0);
        rin.t(usemint) = in.c.mint(in.n.index(usemint));
        if any(rin.t <= 0); error('get_ages_v3.m: zero in rin.t?'); end;
        if any(isinf(rin.t)); error('get_ages_v3.m: Inf in rin.t?'); end; 
        % Note that truet and maxt should not be set at the same time. That
        % would cause trouble. 
        % For C-14, we assume the user put in enough time for saturated
        % samples. 
    end;
    
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
    
    % 5C. Now calculate ages/production rates using scaling data. 
    % Both exposure age and production rate cases, compute integrated N
    % over time for each analysis. 

    % Generate integrated N including muon production

    for a = 2:length(sfa); % Always assume sfa{1} is 'St'. 
        % Do things for each other scaling scheme
        for b = 1:length(in.n.index);
            % Loop over analyses
            % Determine vector of P for this scheme
            % Do spallation integral by hybrid scheme - eqs 39-40 in
            % 2008 appendix
            
            % Generate P0 with site scaling adjustments
            % This is a scalar. 
            eval(['thisP0 = consts.refP_' sfa{a} '(in.n.nindex(b)).*in.s.othercorr(in.n.index(b)).*in.n.sf_thick(b);']);
    
            % Obtain spallogenic surface production rate in each time step. 
            % This is a vector. 
            % This is corrected for thickness/shielding. 
            eval(['pp = thisP0.*sfdata.' sfa{a} '{b};']);

            % Do spallation integral
            % Need to act differently if A == 0, i.e. stable nuclide no
            % erosion. 
            if in.n.A_sp(b) == 0; 
                % Case stable nuclide, no erosion
                int_Nsp = pp.*(sfdata.tmax{b}-sfdata.tmin{b}); % Total production in each time step (atoms/g)
                int_Nmu = in.n.P_mu(b).*(sfdata.tmax{b}-sfdata.tmin{b}); % Same thing
            else
                % Case A > 0
                % Spallation integral
                int_Nsp = -(pp./in.n.A_sp(b)).*(exp(-sfdata.tmax{b}.*in.n.A_sp(b)) - exp(-sfdata.tmin{b}.*in.n.A_sp(b)));
                % Muon integral
                int_Nmu = -(in.n.P_mu(b)./in.n.A_mu(b)).* (exp(-sfdata.tmax{b}.*in.n.A_mu(b)) - exp(-sfdata.tmin{b}.*in.n.A_mu(b)));
                % That gives total production in each time step, corrected
                % for decay to the present time. 
            end;
            % Cumulate. Note recording spallation and muons separately for
            % use in production rate calibration case. 
            eval(['sfdata.cum_N_mu_' sfa{a} '{b} = cumsum(int_Nmu);']);
            eval(['sfdata.cum_N_sp_' sfa{a} '{b} = cumsum(int_Nsp);']);
            eval(['sfdata.cum_N_' sfa{a} '{b} = sfdata.cum_N_mu_' sfa{a} '{b} + sfdata.cum_N_sp_' sfa{a} '{b};']);
        end;
    end;
      
    % sfdata now has cumulated N for each scaling scheme. 
        
    if control.cal == 0;
        % Case determining exposure ages
        for a = 2:length(sfa); % Always assume sfa{1} is 'St'
            % Initialize
            eval(['in.n.t_' sfa{a} ' = zeros(size(in.n.t_St));']);
            eval(['in.n.delt_int_' sfa{a} ' = in.n.t_' sfa{a} ';']);
            eval(['in.n.delt_ext_' sfa{a} ' = in.n.t_' sfa{a} ';']);
            for b = 1:length(in.n.index);
                % Loop over analyses
                % Perform saturation check
                eval(['tempmax = sfdata.cum_N_' sfa{a} '{b}(end);']);
                if in.n.N(b) > tempmax;
                    % Case saturated
                    eval(['in.n.t_' sfa{a} '(b) = 0;']);
                    flag = ['Sample ' char(in.s.sample_name(in.n.index(b))) ' -- ' char(in.n.nuclide{b}) ' appears to be saturated WRT ' char(sfa{a}) ' scaling at this erosion rate.'];
                    in.flags = [in.flags ' <br> ' sprintf('%s\n','') flag]; 
                else
                    % Case not saturated   
                    eval(['in.n.t_' sfa{a} '(b) = interp1(sfdata.cum_N_' sfa{a} '{b},sfdata.tmax{b},in.n.N(b));']);
                end;
            end;
        end;
                    
        
    
    elseif control.cal == 1;
        % Case determining production rates
        % Generate integrated Nsp
        % Determine P for min/max/exact; 
        for a = 2:length(sfa); % Always assume sfa{1} is 'St'.
            % Initialize
            eval(['in.n.calc_P_' sfa{a} ' = zeros(size(in.n.calc_P_St));']);
            eval(['in.n.calc_delP_' sfa{a} ' = in.n.calc_P_' sfa{a} ';']);
            eval(['in.n.calc_minP_' sfa{a} ' = zeros(size(in.n.calc_P_St));']);
            eval(['in.n.calc_delminP_' sfa{a} ' = in.n.calc_minP_' sfa{a} ';']);
            eval(['in.n.calc_maxP_' sfa{a} ' = zeros(size(in.n.calc_P_St));']);
            eval(['in.n.calc_delmaxP_' sfa{a} ' = in.n.calc_maxP_' sfa{a} ';']);
            for b = 1:length(in.n.index);
                % Loop over analyses
                % We are only considering spallogenic production here.
                % Thus, we assume calculated muon production is correct and
                % subtract it from the measured nuclide concentration. 
                % Get reference production rate out of cumulated N.
                % To do this, we need cumulated N_sp and N_mu separately,
                % for the scaling scheme of interest. Get them. 
                eval(['this_cum_N_mu = sfdata.cum_N_mu_' sfa{a} '{b};']); 
                % That's an atoms/g concentration. 
                % Note the below still contains the thickness and topo
                % shielding factors. The only thing divided out is the
                % reference production rate. 
                eval(['this_cum_SF_sp = sfdata.cum_N_sp_' sfa{a} '{b}./consts.refP_' sfa{a} '(in.n.nindex(b));']);
                % That's a nondimensional integrated scaling factor. 
                
                % Now calculate production rates
                if any(b == istruet);
                    % Case exact t; time vector is correct length. 
                    % Obtain target spallogenic N. 
                    this_target_Nsp = in.n.N(b) - this_cum_N_mu(end);  
                    % Note this would fail if ismaxt and is truet, because the time vector would 
                    % have been calculated out to maxt, but that
                    % shouldn't occur. 
                    % Now calculate implied reference production rate. 
                    eval(['in.n.calc_P_' sfa{a} '(b) = this_target_Nsp./this_cum_SF_sp(end);']);
                end;
                if any(b == ismint);
                    % Sites with minimum age give maximum production rate
                    % Obtain target spallogenic N - interpolate Nmu to mint
                    this_target_Nsp = in.n.N(b) - interp1(sfdata.tmax{b},this_cum_N_mu,in.c.mint(in.n.index(b)));
                    % Now calculate implied reference production rate
                    eval(['in.n.calc_maxP_' sfa{a} '(b) = this_target_Nsp./( interp1(sfdata.tmax{b},this_cum_SF_sp,in.c.mint(in.n.index(b))) );']);
                end;
                if any(b == ismaxt);
                    % In this case the time vector should have been
                    % calculated out to maxt. 
                    this_target_Nsp = in.n.N(b) - this_cum_N_mu(end);
                    % Now calculate implied reference production rate
                    eval(['in.n.calc_minP_' sfa{a} '(b) = this_target_Nsp./this_cum_SF_sp(end);']);
                end;

            end;
        end;
    
    end;
    
end; % End loop executed when result type is 'long'

% --------- END ALL TIME-DEPENDENT AGES OR PRODUCTION RATES -------------
    
% 6. Compute linearized uncertainties. Need separate cases for cal and not
% cal. 

if control.cal == 0;
    % Case working with exposure ages. 
    for b = 1:length(sfa);
        sfs = sfa{b}; % String for SF
        
        % Extract fractional external production rate uncertainty for
        % current SS
        eval(['fdelP = (consts.delrefP_' sfs '(in.n.nindex)./consts.refP_' sfs '(in.n.nindex))'';']);
        
        % extract t for current SS
        eval(['tt = in.n.t_' sfs ';']); 
        % Flag saturated
        ok = find(tt > 0); notok = find(tt == 0); % flag saturated ones
        case1 = find(tt > 0 & in.n.A_sp > 0); % Case is an age and A > 0
        case2 = find(tt > 0 & in.n.A_sp == 0); % Case age exists, A = 0
        
        % Allocate
        FP = zeros(size(in.n.index));
        delFP = FP; dtdN = FP; dtdP = FP;
        
        FP(case1) = (in.n.N(case1).*in.n.A_sp(case1))./(1 - exp(-in.n.A_sp(case1).*tt(case1)));
        FP(case2) = in.n.N(case2)./tt(case2);
        
        % Also, use these to compute effective normalized nuclide concentrations for
        % time-dependent scaling schemes, for messing around with in
        % banana plots. 
        
        if ~strcmp(sfs,'St');
            % initialize
            eval(['in.n.Nnorm_' sfs ' = zeros(size(in.n.index));']);
            eval(['in.n.delNnorm_' sfs ' = in.n.Nnorm_Lm;']);
            eval(['in.n.delNnorm_ext_' sfs ' = in.n.Nnorm_Lm;']);

            
            % do easy ones that aren't saturated
            eval(['in.n.Nnorm_' sfs '(ok) = in.n.N(ok)./FP(ok);']);
            eval(['in.n.delNnorm_' sfs '(ok) = in.n.delN(ok)./FP(ok);']);
            eval(['in.n.delNnorm_ext_' sfs '(ok) = sqrt( (in.n.delN(ok)./FP(ok)).^2 + (fdelP(ok).*FP(ok).*in.n.N(ok)./(FP(ok).^2)).^2 );']);

            % If saturated, the production rate is the max N achieved in the
            % time dependent calculation times lambda. 
            for c = 1:length(notok);
                eval(['temp_sat_FP = sfdata.cum_N_' sfs '{notok(c)}(end).*in.n.l(notok(c));']);
                eval(['in.n.Nnorm_' sfs '(notok(c)) = in.n.N(notok(c))./temp_sat_FP;']);
                eval(['in.n.delNnorm_' sfs '(notok(c)) = in.n.delN(notok(c))./temp_sat_FP;']);
                eval(['in.n.delNnorm_ext_' sfs '(notok(c)) = in.n.delNnorm_' sfs '(notok(c));']); % External uncert doesn't depend on P here. 
            end;
        end;
            
        % Continue with error propagation. 
        % If saturated, no error propagation. 
        
        delFP(ok) = fdelP(ok) .* FP(ok);
       
        dtdN(ok) = 1./(FP(ok) - in.n.N(ok).*in.n.A_sp(ok));  
        dtdP(ok) = -in.n.N(ok)./(FP(ok).*FP(ok) - in.n.N(ok).*in.n.A_sp(ok).*FP(ok));
        
        % Note that uncertainty in the decay consts is not considered here.
        % This is just because it is basically never relevant in exposure
        % dating. To do this, however, one could observe that delA = dell 
        % and use dtdA. 

        % make respective delt's
        eval(['in.n.delt_ext_' sfa{b} '(ok) = sqrt( dtdN(ok).^2 .* in.n.delN(ok).^2 + dtdP(ok).^2 .* delFP(ok).^2);']);
        eval(['in.n.delt_int_' sfa{b} '(ok) = sqrt(dtdN(ok).^2 .* in.n.delN(ok).^2);']);
        
        clear ttt fdelP ok case1 case2 FP delFP dtdN dtdP temp_sat_FP
    end;
else
    % Case working with production rates
    for b = 1:length(sfa);
        sfs = sfa{b}; % String for SF
        
        % Again broken into blocks for exact/min/max; should fix that.  
        % Inherit min/max/exact flag variables from above
        % This really needs to be shortened for less redundancy. 
        if ~isempty(istruet);
            ok = istruet;
            % Split cases for A = 0 and not
            case1 = find(in.n.A_sp(ok) > 0);
            case2 = find(in.n.A_sp(ok) == 0);
            % Reset vectors to make sure correct shape
            SFeff = zeros(size(in.n.index)); dPdNsp = SFeff; dPdt = SFeff;
            
            % Determine effective SF and partials
            
            % First deal with case 1
            eval(['SFeff(ok(case1)) = (in.n.Nsp(ok(case1)).*in.n.A_sp(ok(case1)))./(in.n.calc_P_' sfs '(ok(case1)).*in.s.othercorr(in.n.index(ok(case1))).*in.n.sf_thick(ok(case1)).*(1-exp(-in.n.A_sp(ok(case1)).*in.c.truet(in.n.index(ok(case1))))));']);
            dPdNsp(ok(case1)) = in.n.A_sp(ok(case1))./(SFeff(ok(case1)).*in.s.othercorr(in.n.index(ok(case1))).*in.n.sf_thick(ok(case1)).*(1-exp(-in.n.A_sp(ok(case1)).*in.c.truet(in.n.index(ok(case1))))));
            dPdt(ok(case1)) =(in.n.Nsp(ok(case1))).*in.n.A_sp(ok(case1))./(SFeff(ok(case1)).*in.n.sf_thick(ok(case1)).*in.s.othercorr(in.n.index(ok(case1)))).*in.n.A_sp(ok(case1)).*exp(-in.n.A_sp(ok(case1)).*in.c.truet(in.n.index(ok(case1)))).*((1-exp(-in.n.A_sp(ok(case1)).*in.c.truet(in.n.index(ok(case1))))).^-2);
            
            % Now case 2
            eval(['SFeff(ok(case2)) = (in.n.Nsp(ok(case2)).*in.c.truet(in.n.index(ok(case2))))./(in.n.calc_P_' sfs '(ok(case2)).*in.s.othercorr(in.n.index(ok(case2))).*in.n.sf_thick(ok(case2)));']);
            dPdNsp(ok(case2)) = in.c.truet(in.n.index(ok(case2)))./(SFeff(ok(case2)).*in.s.othercorr(in.n.index(ok(case2))).*in.n.sf_thick(ok(case2)));
            dPdt(ok(case2)) = in.n.Nsp(ok(case2))./(SFeff(ok(case2)).*in.s.othercorr(in.n.index(ok(case2))).*in.n.sf_thick(ok(case2)));
            
            % Uncert
            eval(['in.n.calc_delP_' sfs '(ok) = sqrt( (in.n.delN(ok).*dPdNsp(ok)).^2 + (in.c.dtruet(in.n.index(ok)).*dPdt(ok)).^2 );']);    
        end;
        
        if ~isempty(ismint);
            % mint corresponds to maxp, remember
            ok = ismint;
            % Split cases for A = 0 and not
            case1 = find(in.n.A_sp(ok) > 0);
            case2 = find(in.n.A_sp(ok) == 0);
            % Reset vectors to make sure correct shape
            SFeff = zeros(size(in.n.index)); dPdNsp = SFeff; dPdt = SFeff;
            
            % Determine effective SF and partials
            
            % First deal with case 1
            eval(['SFeff(ok(case1)) = (in.n.min_Nsp(ok(case1)).*in.n.A_sp(ok(case1)))./(in.n.calc_maxP_' sfs '(ok(case1)).*in.s.othercorr(in.n.index(ok(case1))).*in.n.sf_thick(ok(case1)).*(1-exp(-in.n.A_sp(ok(case1)).*in.c.mint(in.n.index(ok(case1))))));']);
            dPdNsp(ok(case1)) = in.n.A_sp(ok(case1))./(SFeff(ok(case1)).*in.s.othercorr(in.n.index(ok(case1))).*in.n.sf_thick(ok(case1)).*(1-exp(-in.n.A_sp(ok(case1)).*in.c.mint(in.n.index(ok(case1))))));
            dPdt(ok(case1)) =(in.n.min_Nsp(ok(case1))).*in.n.A_sp(ok(case1))./(SFeff(ok(case1)).*in.n.sf_thick(ok(case1)).*in.s.othercorr(in.n.index(ok(case1)))).*in.n.A_sp(ok(case1)).*exp(-in.n.A_sp(ok(case1)).*in.c.mint(in.n.index(ok(case1)))).*((1-exp(-in.n.A_sp(ok(case1)).*in.c.mint(in.n.index(ok(case1))))).^-2);
            
            % Now case 2
            eval(['SFeff(ok(case2)) = (in.n.min_Nsp(ok(case2)).*in.c.mint(in.n.index(ok(case2))))./(in.n.calc_maxP_' sfs '(ok(case2)).*in.s.othercorr(in.n.index(ok(case2))).*in.n.sf_thick(ok(case2)));']);
            dPdNsp(ok(case2)) = in.c.mint(in.n.index(ok(case2)))./(SFeff(ok(case2)).*in.s.othercorr(in.n.index(ok(case2))).*in.n.sf_thick(ok(case2)));
            dPdt(ok(case2)) = in.n.min_Nsp(ok(case2))./(SFeff(ok(case2)).*in.s.othercorr(in.n.index(ok(case2))).*in.n.sf_thick(ok(case2)));
            
            % Uncert
            eval(['in.n.calc_delmaxP_' sfs '(ok) = sqrt( (in.n.delN(ok).*dPdNsp(ok)).^2 + (in.c.dmint(in.n.index(ok)).*dPdt(ok)).^2 );']);    
        end;
        
        if ~isempty(ismaxt);
            % maxt corresponds to minp, remember
            ok = ismaxt;
            % Split cases for A = 0 and not
            case1 = find(in.n.A_sp(ok) > 0);
            case2 = find(in.n.A_sp(ok) == 0);
            % Reset vectors to make sure correct shape
            SFeff = zeros(size(in.n.index)); dPdNsp = SFeff; dPdt = SFeff;
            
            % Determine effective SF and partials
            
            % First deal with case 1
            eval(['SFeff(ok(case1)) = (in.n.max_Nsp(ok(case1)).*in.n.A_sp(ok(case1)))./(in.n.calc_minP_' sfs '(ok(case1)).*in.s.othercorr(in.n.index(ok(case1))).*in.n.sf_thick(ok(case1)).*(1-exp(-in.n.A_sp(ok(case1)).*in.c.maxt(in.n.index(ok(case1))))));']);
            dPdNsp(ok(case1)) = in.n.A_sp(ok(case1))./(SFeff(ok(case1)).*in.s.othercorr(in.n.index(ok(case1))).*in.n.sf_thick(ok(case1)).*(1-exp(-in.n.A_sp(ok(case1)).*in.c.maxt(in.n.index(ok(case1))))));
            dPdt(ok(case1)) =(in.n.max_Nsp(ok(case1))).*in.n.A_sp(ok(case1))./(SFeff(ok(case1)).*in.n.sf_thick(ok(case1)).*in.s.othercorr(in.n.index(ok(case1)))).*in.n.A_sp(ok(case1)).*exp(-in.n.A_sp(ok(case1)).*in.c.maxt(in.n.index(ok(case1)))).*((1-exp(-in.n.A_sp(ok(case1)).*in.c.maxt(in.n.index(ok(case1))))).^-2);   
            
            % Now case 2
            eval(['SFeff(ok(case2)) = (in.n.max_Nsp(ok(case2)).*in.c.maxt(in.n.index(ok(case2))))./(in.n.calc_minP_' sfs '(ok(case2)).*in.s.othercorr(in.n.index(ok(case2))).*in.n.sf_thick(ok(case2)));']);
            dPdNsp(ok(case2)) = in.c.maxt(in.n.index(ok(case2)))./(SFeff(ok(case2)).*in.s.othercorr(in.n.index(ok(case2))).*in.n.sf_thick(ok(case2)));
            dPdt(ok(case2)) = in.n.max_Nsp(ok(case2))./(SFeff(ok(case2)).*in.s.othercorr(in.n.index(ok(case2))).*in.n.sf_thick(ok(case2)));
            
            % Uncert
            eval(['in.n.calc_delminP_' sfs '(ok) = sqrt( (in.n.delN(ok).*dPdNsp(ok)).^2 + (in.c.dmaxt(in.n.index(ok)).*dPdt(ok)).^2 );']);    
        end;
        
        
    end;
end;
        
% ---------- DONE WITH ERROR PROPAGATION ----------------------

% 6. Return output data, added to input data

out = in;