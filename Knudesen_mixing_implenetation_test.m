% test script to compare our code with Knudesen implementation

scenario = 'step';
tdata = make_test_data(scenario,1);
Nlogical = [true(1,2) false(1,1)];  % only 10Be and 14C

%% Priors -----------------------------------------------------------------
T =  [1.2e3,1.7e3]; % time of step change OR spike in yrs [min,max]
E =  [0,8e2];      % range of expected erosion rates in mm/ka  [min,max]
LOSS = [0,200];     % loss of soil in cm [min,max], can be commented if no spike model
CHG  = [0, 50];     % change factor of erosion rate, can be commented if no samestep model. Also serves as scale factor for curve model
[prior_range,var_names] = make_prior_and_varnames(scenario,T,E,LOSS,CHG,1,tdata.steps);

[consts,Nmu] = make_constants();

sp = Cronus_v3_spallation(tdata.lat,tdata.lon,tdata.altitude,consts);   % get sample parameters (surface procution, pressure)

sp.mix = true;    % make flag
sp.mixing = 3e-5; % mixing diffusivity
sp.dzmix = 0.1;   % mixing depth (m)

mtest = [tdata.t'; tdata.e'; tdata.changeVariable'];

[N1,N14] = Nforward_discretized_FEM(tdata.e(1),T,sp,consts,scenario,Nlogical);

%% Knudsen 

CNprop.PBe = 4; %10Be production due to spallation in atoms/(g*yr)
TBe = 1.387e6; %Half-life of 10Be
TAl = 0.705e6; %Half-life of 26Al
CNprop.lambda_Be = log(2)/TBe;  %Decay constant of 10Be
CNprop.lambda_Al = log(2)/TAl; %Decay constant of 26Al
CNprop.rho = 1.9; %Soil density in g/cm3
CNprop.pratio = 6.97; %Al/Be production ratio
CNprop.pr_fm_Be = 0.005; %Be production due to fast muons (of total p)
CNprop.pr_fm_Al = 0.006; %Al production due to fast muons (of total p)
CNprop.pr_nmc_Be = 0.015; %Be production due to negative muon capture (of total p)
CNprop.pr_nmc_Al = 0.018; %Al production due to negative muon capture (of total p)
CNprop.Lspal = 150; %Attenuation length for fast neutrons
CNprop.Lnmc = 1500; %Attenuation length for negative muon capture
CNprop.Lfm = 4320;  %Attenuation length for fast muon

%other parameters
data.P10spal = 4.0;  %10Be production due to spallation in atoms/(g*yr)
data.P10muon = 0.06+0.02; %10Be production due to muons in atoms/(g*yr)
data.erate_g = 10e-6; %Erosion rate glacial (m/y)
data.erate_ng = 1e-6; %Erosion rate non-glacial (m/y)
data.mixing = 3e-5; %mixing
data.dzmix = 1; %mixing mixing depth
data.age = 1000e3; %time period of model
data.d18O = 4.0; %d18O threshold to seperate glacial/interglacial
data.eini = 5e-6; %Initial steady-state erosion rate

%Geometry
D = 25; %depth of profile
Nn = 100; %nodes in regolith

%node dispribution
z = linspace(0,1,Nn);
z = D*z.^2;
ze = .5*(z(1:(Nn-1))+z(2:Nn));
Ne = Nn - 1;

%model props
dt = 5e3;
theta = 1; %time integration

%make dmix vector
z0 = 1.8;
dz = 0.1;
zb = 0.0;
nomix = zeros(Ne,1);
dmix = data.mixing*exp(-ze(:)/data.dzmix);

%compute steady state profile with mixing
[N2,N26ss,P10K] = steady_state(z,CNprop,data.P10spal,data.P10muon,data.eini,dmix,1);

%%

plot(N1{1},z, 'r-')
hold on
plot(N2,z, 'k-')
plot(N14{1},z,'b-')
legend({ "my code", "Knudsen", "N14"})


%% check N14 and N10 changes with transient erosion

figure()
tiledlayout(1,3)
nexttile
hold on
mixing  = logspace(-7,-2);
cc = viridis(length(mixing));
for i = 1:length(mixing)
    dmix = mixing(i)*exp(-ze(:)/data.dzmix);
    [N2,N26ss,P10K] = steady_state(z,CNprop,data.P10spal,data.P10muon,data.eini,dmix,1);
    plot(N2,z,'Color',cc(i,:))
end

%%

nexttile
hold on
for i = 1:length(mixing)
    dmix = mixing(i)*exp(-ze(:)/data.dzmix);
    sp.mix = true;    % make flag
    sp.mixing = mixing(i); % mixing diffusivity
    sp.dzmix = 1;   % mixing depth (m)
    
    mtest = [tdata.t'; tdata.e'; tdata.changeVariable'];
    
    [N1,N14] = Nforward_discretized_FEM(tdata.e(1),T,sp,consts,scenario,Nlogical);

    plot(N1{1},z,'Color',cc(i,:))
end

%%
nexttile
hold on
for i = 1:length(mixing)
    dmix = mixing(i)*exp(-ze(:)/data.dzmix);
    sp.mix = true;    % make flag
    sp.mixing = mixing(i); % mixing diffusivity
    sp.dzmix = 1;   % mixing depth (m)
    
    mtest = [tdata.t'; tdata.e'; tdata.changeVariable'];
    
    [N1,N14] = Nforward_discretized_FEM(tdata.e(1),T,sp,consts,scenario,Nlogical);

    plot(N14{1},z,'Color',cc(i,:))
end
