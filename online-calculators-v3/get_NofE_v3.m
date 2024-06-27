function out = get_NofE_v3(in,consts)

% This computes a block of nuclide concentrations as a function of erosion
% rate from a validated data structure. 
%
% The purpose of this script is to service some external piece of code that
% wants to invert something for an erosion rate, but doesn't want to do the
% production rate calculations. 
% 
% input args in,consts
%
% in is a structure that is the output of validate_v3. 
% this has s,n sub-structures for sample, nuclide data. 
% 
% consts is consts structure from make_consts_v3; optional
%
% Returns output arg out which is just the input arg in with a lot of 
% new fields added to it. 
%
% Greg Balco
% Berkeley Geochronology Center
% May, 2022
% 
% Development. Not licensed for use or distribution. 
%
% This is simplified from get_erates_v3.m. 

% Should pass consts to minimize file loads, but it is not a very big file
% so, whatever. Passing consts is needed if one wants to use non-default
% production rates or other parameters. 

if nargin < 2
    load consts_v3;
end

% Also load precalculated muon data needed for forward integration
load([consts.datadir 'Nmu.mat']);

% Version
in.version_NofE = '3.0-dev';
in.version_muons = Nmu.version;

% Initialize

in.flags = [];

% ---------------- COMPUTE SAMPLE-SPECIFIC PARAMETERS ------------------

% 1. Compute sample-specific parameters: pressure, various SFs

% Calculate pressures with minimal number of calls to ERA40atm.m.
% In this section is where we would pull an elevation history from an 
% isostatic rebound model if such existed. 

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
for a = 1:length(in.n.index)
    % Loop over all nuclide measurements
    % Determine what index to use in consts file
    in.n.nindex(a) = find(strcmp(in.n.nuclide{a},consts.nuclides));
    % This isn't robust against errors in the nuclide specification, but 
    % shouldn't need to be if validate_v3_input is doing its job correctly. 
end

% Decay constants

in.n.l = consts.l(in.n.nindex)'; % Transpose needed for some reason.
in.n.dell = consts.dell(in.n.nindex)'; % Likewise. 

% Compute input erosion rates

in.s.Egcm2 = in.s.E .* in.s.rho; % now in g/cm2

% ----------- COMPUTE NTD PREDICTED N ---------------------------------- 

% Total production rate on St
in.n.Psp_St = in.n.sf_St.*in.s.othercorr(in.n.index).*in.n.sf_thick.*consts.refP_St(in.n.nindex)';


% Allocate
in.n.Npred_St = zeros(size(in.n.index));

% Compute the predicted nuclide concentrations

for a = 1:length(in.n.N)               
    % Define which muon thing to use
    eval(['this_Nmu = Nmu.' in.n.nuclide{a} ';'])   
    in.n.Npred_St(a) =  in.n.Psp_St(a)./(in.n.l(a) + in.s.Egcm2(in.n.index(a))./in.n.Lsp(a)) + intNmu(in.s.Egcm2(in.n.index(a)),in.s.pressure(in.n.index(a)),Nmu.pp,Nmu.logee,this_Nmu);    
end
    
% ---------- DONE WITH NTD EROSION RATES ----------------------

% ---------- NOW GENERATE TIME-DEPENDENT SF ------------------------------

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

% 5B. Get SF time series from precalculated results. Loop for various SF. 
% Relies on get_LmSF.m, get_LSDnSF.m
for a = 1:length(sfa)
    if ~strcmp(sfa{a},'St')
        eval(['sfdata = get_' sfa{a} 'SF(sfdata,consts.datadir);'])
    end
end

% sfdata now has:
%   cell array LSDn
%   cell array Lm
% etc. These have vectors of "scaling factor". 

% Add sfdata to output structure so it is available later if needed

in.f = sfdata;

% 5C. Now calculate predicted nuclide concentrations using scaling data. 
for a = 2:length(sfa) % Assume sfa(1) = 'St', which we already did
    % Initialize
    eval(['in.n.Npred_' sfa{a} '= zeros(size(in.n.N));']);
    for b = 1:length(in.n.N)
        % This loop goes through each nuclide concentration measurement
        % Generate P0 with site scaling adjustments
        % This is a scalar and includes thickness/shielding corrections. 
        eval(['thisP0 = consts.refP_' sfa{a} '(in.n.nindex(b)).*in.s.othercorr(in.n.index(b)).*in.n.sf_thick(b);']);

        % Obtain spallogenic surface production rate in each time step. 
        % This is a vector. 
        eval(['pp = thisP0.*sfdata.' sfa{a} '{b};']);

        % Define which thing to use
        eval(['this_Nmu = Nmu.' in.n.nuclide{b} ';'])
       
        % Define function 
        % This has 
        % 1. the muon integral, 
        % 2. the steady-state concentration at the beginning of
        % the time-dependent part, decayed until the present time,
        % and 3. the concentration developed during
        % the time-dependent part. 
        % (2) is probably very small. 
        efun = @(x) (intNmu(x,in.s.pressure(in.n.index(b)),Nmu.pp,Nmu.logee,this_Nmu) ...
            + exp(-rin.t(b).*in.n.l(b)).*pp(end).*exp(-rin.t(b).*x./in.n.Lsp(b))./(in.n.l(b)+x/in.n.Lsp(b)) ...
            + sum(-(pp./(in.n.l(b)+x/in.n.Lsp(b))).*(exp(-sfdata.tmax{b}.*(in.n.l(b)+x/in.n.Lsp(b))) - exp(-sfdata.tmin{b}.*(in.n.l(b)+x/in.n.Lsp(b))))) );
             
        eval(['in.n.Npred_' sfa{a} '(b) = efun(in.s.Egcm2(in.n.index(b)));']); 
    end
end 


% --------- END ALL TIME-DEPENDENT N CALCS -------------
    
% Finaly, return output data, added to input data

out = in;