function out = make_consts_v3()

% This function creates and saves a structure with relevant constants and
% external data for the exposure age and erosion rate calculators.  
%
% This is for v3 code that ingests multiple-nuclide data. 
%
% Syntax: make_consts_v3
% (no arguments)
%
%
% Greg Balco - Berkeley Geochronology Center - February, 2015.
% Development, not licensed for use or distribution. 
%
% Updated Dec 2016 to switch to simplified muon scheme. 
% Updated many times subsequently. 

consts.version = datestr(now,'yyyy-mm-dd');
consts.prepdate = fix(clock);

% Determine what machine we are running on. 
% Note: this means you can't just upload the mat-file containing the 
% constants; you need to upload this m-file and execute it. 
% Also, you need to edit a couple of places depending on whether you are on
% hess or something else. 
%
% Also set machine-specific data file locations.
% 
% Note: hard machine urls are needed in the following places:
%   1. location of image files 
%   2. form link in calibration return input page
% 


consts.isLocal = 1;
consts.scratchdir = '.\scratch\';
consts.RCfname = '.\data\PavonDip.mat';
consts.Sfname = '.\data\LiftonS2015.mat';
consts.datadir = '.\data\';
% consts.plotURL = ['https://' consts.machine_name '/scratch/'];   

% Determine if running Octave

% consts.is_octave = 1;
% temp = ver('Octave');
% if isempty(temp); consts.is_octave = 0; end

% Be-10 standardization info.

% Standards comparison/conversion lookup tables

% For Be-10. A zero placeholder is allowed, but suppresses 
% calculations even if corresponding concentrations are entered. 

consts.be_stds_names = {'07KNSTD','KNSTD','NIST_Certified','LLNL31000',...
    'LLNL10000','LLNL3000','LLNL1000','LLNL300','NIST_30000','NIST_30200',...
    'NIST_30300','NIST_30600','NIST_27900','S555','S2007','BEST433',...
    'BEST433N','S555N','S2007N','STD11','0','NIST_30500','SMDBe12'};
consts.be_stds_cfs = [1.0000 0.9042 1.0425 0.8761 0.9042 0.8644 0.9313 ...
    0.8562 0.9313 0.9251 0.9221 0.9130 1 0.9124 0.9124 0.9124 1 1 1 1 0 0.9124 1];

% Same for Al-26. A zero placeholder is also allowed, same idea.  

consts.al_stds_names = {'KNSTD','ZAL94','SMAL11','0','ZAL94N','ASTER','Z92-0222','SMDAl11'};
consts.al_stds_cfs = [1.0000 0.9134 1.021 0 1 1.021 1 1.018];

% Handled slightly differently for He/Ne, here we define reference standard
% materials and a reference concentration for each that we normalize to
% based on a lab's results for the standard. "NONE" is also allowed, in
% which case no restandardization happens. Thus, in contrast to Al-26 and
% Be-10, you can calculate He-3 and Ne-21 exposure ages without
% standardization information. For now. 

consts.he_stds_names = {'CRONUS-P','NONE'};
consts.he_stds_ref = [5.02e9 0];

consts.ne_stds_names = {'CRONUS-A','CREU-1','NONE'};
consts.ne_stds_ref = [320e6 348e6 0];

% Also intended to be handled this way for C-14, but not used at present. 

consts.c_stds_names = {'CRONUS-A','NONE'};
consts.c_stds_ref = [690000 0]; % What to use for this? 

% Define minerals allowed for various nuclides. Obviously, there needs to
% be code in get_ages_v3 for each one of these, the lookup tables below 
% need to be in agreement, and the wrapper code needs to agree as well. 
% You can't just randomly add things here. Note that this is kind of
% redundant with the nuclide/mineral pair names below...could be improved. 

consts.ok_mins_10 = {'quartz','pyroxene'};
consts.ok_mins_26 = {'quartz'};
consts.ok_mins_3 = {'quartz','pyroxene','olivine'};
consts.ok_mins_21 = {'quartz'};
consts.ok_mins_14 = {'quartz'};

% Define scaling schemes; mostly this defines loops in XML and HTML dumps,
% it is redundantly defined in get_ages_v3. Could be improved. 

consts.sschemes = {'St','Lm','LSDn','Ag'}; 
consts.sspropernames = {'"St": Lal(1991)/Stone(2000)','"Lm": Time-variable Lal/Stone',...
    '"LSDn": Fully loaded Lifton/Sato','"Ag": Some implementation of Argento (2014)'};

% Other odd bits
consts.default_yr = 2010; % Default year sample was collected

% Decay constants

% Be-10
consts.l10 = -log(0.5)./1.387e6; % Chmeleff/Korschinek value
dldt = -log(0.5).*(1.387e6^-2);
consts.dell10 = sqrt((dldt.*0.012e6)^2); 

% Al-26 -- value compatible with Nishiizumi standards
% lambda = 9.83e-7 --> t(1/2) = 0.705e6 yr
% See Nishiizumi (2004) for details.

consts.l26 = 9.83e-7;
consts.dell26 = 2.5e-8;

% C-14 decay constant

consts.l14 = -log(0.5)./5730; 
consts.dell14 = consts.l14.*0.005; % ????

% Format nuclide-related data as lookup table

% Define nuclide/mineral pairs
consts.nuclides = {'N3quartz','N3olivine','N3pyroxene','N10quartz','N14quartz','N21quartz','N26quartz','N10pyroxene'};
consts.properName = {'He-3 (qtz)','He-3 (ol)','He-3 (px)','Be-10 (qtz)','C-14 (qtz)','Ne-21 (qtz)','Al-26 (qtz)','Be-10 (px)'};

% Production ratios.  
% consts.R2610q = 6.75; % Spallogenic 26/10 ratio in quartz
% I don't think the above is needed any more.
% consts.R2110q = 3.89; % Spallogenic 21/10 ratio in quartz, normalized to CRONUS-A = 320e6
% I believe the above are now obsolete. 

% Production rates for use with Lal-Stone scheme. 

consts.P10q_St = 4.086; % From recalibrating with CRONUS primary 20161204
consts.delP10q_St = 4.086.*0.079; % bootstrap from secondary cal data set

consts.P26q_St = 28.535; % From recalibrating with CRONUS primary 20161204
consts.delP26q_St = 28.535.*0.104;

consts.P3q_St = 115.7; % From Vermeesch data
consts.delP3q_St = 10; % ??

consts.P3op_St = 119.6; % From recalibrating w/CRONUS primary 20160415
consts.delP3op_St = 119.6 * 0.11; % bootstrap from primary data set scatter

consts.P14q_St = 12.1; % From Brent CRONUS-A data 20161204
consts.delP14q_St = 12.1 * 0.05; % ?? 

consts.P21q_St = 16.896; % From SPICE data
consts.delP21q_St = 1.033;

consts.P10px_St = 3.72; % From Bergelin saturation data
consts.delP10px_St = 3.72 * 0.06; % Basically made up to match Be-10 in quartz

% Assemble into 8-element vector. 
% This has reference production rates in order: P3q, P3ol, P3px, P10q,
% P14q, P21q, P26q, P10px. Total 8 elements. 
consts.refP_St = [consts.P3q_St consts.P3op_St consts.P3op_St consts.P10q_St consts.P14q_St ...
    consts.P21q_St consts.P26q_St consts.P10px_St];
consts.delrefP_St = [consts.delP3q_St consts.delP3op_St consts.delP3op_St consts.delP10q_St consts.delP14q_St ...
    consts.delP21q_St consts.delP26q_St consts.delP10px_St];

% Simplified muon scheme, uses model 1A with alpha = 1;
% This is output from approx_1A_with_exp_for_calcs.m
% Revised 20200826
% Be-10 (quartz): P0 = 0.07345; Leff = 299.2
% Al-26: P0 = 0.67644; Leff = 288.0
% C-14: P0 = 3.06690; Leff = 267.8
% Ne-21: P0 = 0.20281; Leff = 482.2

% Added 20240302
% Be-10 in pyroxene: P0 = 0.05983; Leff = 311.0

% Assume no muon production for He-3. This needs to be fixed. 

consts.version_muons = '1A, alpha = 1';
% Also in 8-element vectors corresponding to nuclide/target options
consts.Pmu0 = [0 0 0 0.07345 3.06690 0.20281 0.67644 0.05983]; 
consts.L_mu_atm = [1 1 1 299.2 267.8 482.2 288.0 311.0]; 

consts.Lmu.a = 0.01036; consts.Lmu.b = -9.697e-6; % Effective near-surface attenuation lengths
% Note that these might be inaccurate for Ne-21. That might oughta be fixed. 

% Define order of "half-lives" for axes of 2-nuclide plots.The plotting code 
% can't just figure out which half-life is shorter, because that won't
% work for He-3 in quartz. 

consts.splot_order = [1 5 5 4 2 5 3 4]; % Controls axis orientation in splots

% Decay constants as 8-element array

consts.l = [0 0 0 consts.l10 consts.l14 0 consts.l26 consts.l10];
consts.dell = [0 0 0 consts.dell10 consts.dell14 0 consts.dell26 consts.l10];

% Other scaling methods

% Lm -- these are dimensional production rates
 
% He-3 in quartz is assumed same as for St. Needs work. 
consts.refP_Lm(1) = consts.refP_St(1);
consts.delrefP_Lm(1) = consts.delrefP_St(1);

% He-3 in px/ol is calibrated.
consts.refP_Lm(2) = 119.6; consts.refP_Lm(3) = consts.refP_Lm(2);
consts.delrefP_Lm(2) = 119.6.* 0.11; consts.delrefP_Lm(3) = consts.delrefP_Lm(2);

% Be-10 is calibrated.
consts.refP_Lm(4) = 4.208; % 20161204 calibration from CRONUS primary. 
consts.delrefP_Lm(4) = 4.208.*0.075; % from scatter

% C-14 is calibrated
consts.refP_Lm(5) = 12.1; % Brent CRONUS-A data 20161204
consts.delrefP_Lm(5) = 12.1.*0.05; % ???

% Ne-21 is calibrated from SPICE data
consts.refP_Lm(6) = 16.243;
consts.delrefP_Lm(6) = 0.994;

% Al-26 is calibrated
consts.refP_Lm(7) = 29.416; % CRONUS primary calibration 20161204
consts.delrefP_Lm(7) = 29.416.*0.094; 

% Be-10 in px
consts.refP_Lm(8) = 3.72; % From Bergelin saturation data
consts.delrefP_Lm(8) = 3.72 * 0.06; % Basically fabricated to match Be-10/quartz 

% LSDn -- these are nondimensional correction factors, not dimensional
% production rates. 

% He-3 in quartz is just assumed to be 1 +/- 0.1. Needs work. 
consts.refP_LSDn(1) = 1;
consts.delrefP_LSDn(1) = 0.1;

% He-3 in px/ol is calibrated
consts.refP_LSDn(2) = 1.323; 
consts.delrefP_LSDn(2) = 1.323 .* 0.11; % bootstrap from primary CDS scatter
consts.refP_LSDn(3) = consts.refP_LSDn(2);
consts.delrefP_LSDn(3) = consts.delrefP_LSDn(2);

% Be-10 is calibrated
consts.refP_LSDn(4) = 0.849; % 20161204 Be-10 calibration from CRONUS primary
consts.delrefP_LSDn(4) = 0.849.*0.059;

% C-14 calibrated from Brent's CRONUS-A data 20161204
consts.refP_LSDn(5) = 0.725;
consts.delrefP_LSDn(5) = 0.725.*0.05; % ???

% Ne-21 is calibrated from SPICE data
consts.refP_LSDn(6) = 1.272;
consts.delrefP_LSDn(6) = 0.079;  

% Al-26 is calibrated
consts.refP_LSDn(7) = 0.819; % 20161204 Al-26 calibration from CRONUS primary
consts.delrefP_LSDn(7) = 0.819 .* 0.086; % Bootstrap

% Be-10 in px
consts.refP_LSDn(8) = 0.680; % From Bergelin saturation data
consts.delrefP_LSDn(8) = 0.680 * 0.06; % This is basically made up to match Be-10/quartz

% This is for plotting He retention curves. 

consts.RHe = generate_He_model_retention();

% This is for plotting Antarctic C-14 saturation plots. 

consts.sat14 = generate_C14_isochrons(0:50:4000,consts);

save consts_v3 consts

if consts.isLocal == 1
    % Don't print anything to output if running as part of something else on server
    disp(['Constants version ' consts.version]);
    disp('Saved'); 
end

