function out = generate_C14_isochrons(elv,consts)


% This generates C-14, etc. saturation curves for Antarctica (e.g., Rc =
% 0) at various elevations. Also isochrons. 

% Set up input data

pressure = antatm(elv);

in.lat = -77.88.*ones(size(elv));
in.long = 160.94.*ones(size(elv)); % Actually this is location of CRONUS-A
in.t = 100000.*ones(size(elv));
in.yr = 2010.*ones(size(elv));

sfdata = get_DipRc(in,consts.RCfname,consts.Sfname);

for a = 1:length(elv);
    sfdata.nuclide{a} = 'N14quartz';
end;

sfdata.pressure = antatm(elv);

% Get scaling

sfdata = get_LSDnSF(sfdata,consts.datadir);
P_mu = consts.Pmu0(5).* exp((1013.25-sfdata.pressure)./consts.L_mu_atm(5));


%% Get P
l14 = consts.l(5);
for a = 1:length(elv);
    P14 = sfdata.LSDn{a}; % keep Pref out of it for now... .*consts.refP_LSDn(5);
    % Spallation integral
    int_Nsp = -(P14./l14).*(exp(-sfdata.tmax{a}.*l14) - exp(-sfdata.tmin{a}.*l14));
    % Muon integral
    int_Nmu = -(P_mu(a)./l14).* (exp(-sfdata.tmax{a}.*l14) - exp(-sfdata.tmin{a}.*l14));
    % Cumulate
    sfdata.cum_N_mu{a} = cumsum(int_Nmu);
    sfdata.cum_N_sp{a} = cumsum(int_Nsp);
    sfdata.cum_N{a} = sfdata.cum_N_mu{a} + sfdata.cum_N_sp{a};
end;

%% Extract various age curves; keep spallation and muons separate so 
% can be adjusted later for calibrated P_sp. 

for a = 1:length(elv);
    mu.sat14(a) = sfdata.cum_N_mu{a}(end);
    mu.t2k(a) = interp1(sfdata.tmax{a},sfdata.cum_N_mu{a},2000);
    mu.t4k(a) = interp1(sfdata.tmax{a},sfdata.cum_N_mu{a},4000);
    mu.t6k(a) = interp1(sfdata.tmax{a},sfdata.cum_N_mu{a},6000);
    mu.t8k(a) = interp1(sfdata.tmax{a},sfdata.cum_N_mu{a},8000);
    mu.t10k(a) = interp1(sfdata.tmax{a},sfdata.cum_N_mu{a},10000);
    mu.t15k(a) = interp1(sfdata.tmax{a},sfdata.cum_N_mu{a},15000);
    mu.t20k(a) = interp1(sfdata.tmax{a},sfdata.cum_N_mu{a},20000);  
    sp.sat14(a) = sfdata.cum_N_sp{a}(end);
    sp.t2k(a) = interp1(sfdata.tmax{a},sfdata.cum_N_sp{a},2000);
    sp.t4k(a) = interp1(sfdata.tmax{a},sfdata.cum_N_sp{a},4000);
    sp.t6k(a) = interp1(sfdata.tmax{a},sfdata.cum_N_sp{a},6000);
    sp.t8k(a) = interp1(sfdata.tmax{a},sfdata.cum_N_sp{a},8000);
    sp.t10k(a) = interp1(sfdata.tmax{a},sfdata.cum_N_sp{a},10000);
    sp.t15k(a) = interp1(sfdata.tmax{a},sfdata.cum_N_sp{a},15000);
    sp.t20k(a) = interp1(sfdata.tmax{a},sfdata.cum_N_sp{a},20000);  
end;

out.elv = elv;
out.sp = sp;
out.mu = mu;

% out.sat14 = sat14;
% out.t2k = t2k;
% out.t4k = t4k;
% out.t6k = t6k;
% out.t8k = t8k;
% out.t10k = t10k;
% out.t15k = t15k;
% out.t20k = t20k;
%  





