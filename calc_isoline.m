function [p1,p2,p1up,p1low,p2up,p2low] = calc_isoline(SAMS,DEM,t,scenario)
% calculates the erosion rates for a range of timings of a step change, or
% the soil loss in a spike scenario
% Input: 
%       - SAMS: Sample structure created by 'cosmosampleread'
%       - DEM: DEM to outline basins (I should add something that makes
%       this available for bedrock point samples).
%        - t: time vector for which to calculate isoline
%        - scenario: 'step'or 'spike'
% Richard Ott, 2024


Nobs = [vertcat(SAMS.N10) ; vertcat(SAMS.N14)];
dNobs= [vertcat(SAMS.N10sigma); vertcat(SAMS.N14sigma)];
Nlogical = [~isnan(vertcat(SAMS.N10)) ~isnan(vertcat(SAMS.N14)) false(length(vertcat(SAMS.N10)),1)];

%% outline basins for binning (currently CosmoTools isnt properly integrated, should be improved)

SAMS = cosmowatersheds(SAMS,DEM);

% get median lat, lon, altitude for production calculation
lat = arrayfun(@(x) median(x.WSLat), SAMS);
lon = arrayfun(@(x) median(x.WSLon), SAMS);
alt = arrayfun(@(x) median(x.WSDEM.Z(:),'omitnan'), SAMS);

%% Constants

[consts,Nmu] = make_constants();

%% Production rates

sp = Cronus_v3_spallation(lat,lon,alt,consts);   % get sample parameters (surface procution, pressure)

%% Calculate nuclide ratios

n = length(vertcat(SAMS.N10));
p1 = nan(n,length(t));  % E1 is step and spike scneario
p2 = nan(n,length(t));  % E2/loss in step and spike scenarios
p1up = nan(n,length(t));    
p2up = nan(n,length(t));
p1low = nan(n,length(t));   
p2low = nan(n,length(t));

switch scenario
    case 'step'        
        x0 = [1e1,1e3];
        UB = [1e3,1e4];   % Upper bound of erosion rates in m/a
    case 'spike'   
        x0 = [50,50];
        UB = [1e3,1e3];   % Upper bound of erosion rates in m/a
end
LB = [0,0];         % Lower bound of erosion rates in m/a   
% Define options for optimization  
options = optimset('MaxIter',5e4,'TolFun',1e3);            % These options may need to be tuned specifically to your problem

for i = 1:n
    samppars.pressure = sp.pressure(i);
    samppars.P10spal = sp.P10spal(i);
    samppars.P14spal = sp.P14spal(i);
    samppars.P26spal = sp.P26spal(i);
    
    iobs = [Nobs(i) Nobs(i+n)];
    idobs = [dNobs(i) dNobs(i+n)];
    for j = 1:length(t)
        

        % average solution ------------------------------------------------
        switch scenario
            case 'step'
                fun = @(x) sum(abs(Nforward_discretized([x(1) x(2)],[t(j) 0],samppars,consts,Nmu,scenario,Nlogical(i,:)) - iobs'));
            case 'spike'
                fun = @(x) sum(abs(Nforward_discretized(x(1),[t(j) 0],samppars,consts,Nmu,scenario,Nlogical(i,:),x(2)) - iobs'));
        end
        % minimize isfit
        [sol,~,exitflag] = fminsearchbnd(fun,x0,LB,UB,options);  

        if exitflag == 0 
            sol = sol*nan;  % set values of failed optimizations to nan
        end

        p1(i,j) = sol(1); p2(i,j) = sol(2);
        
        % lower bound -----------------------------------------------------
        switch scenario
            case 'step'
                fun = @(x) sum(abs(Nforward_discretized([x(1) x(2)],[t(j) 0],samppars,consts,Nmu,'step',Nlogical(i,:)) - iobs'+idobs'));
            case 'spike'
                fun = @(x) sum(abs(Nforward_discretized(x(1),[t(j) 0],samppars,consts,Nmu,scenario,Nlogical(i,:),x(2)) - iobs'+idobs'));
        end
        % minimize isfit
        [sol,~,exitflag] = fminsearchbnd(fun,x0,LB,UB,options);  

        if exitflag == 0 
            sol = sol*nan;  % set values of failed optimizations to nan
        end

        p1low(i,j) = sol(1); p2low(i,j) = sol(2);

        % upper bound -----------------------------------------------------
        switch scenario
            case 'step'
                fun = @(x) sum(abs(Nforward_discretized([x(1) x(2)],[t(j) 0],samppars,consts,Nmu,'step',Nlogical(i,:)) - iobs'-idobs'));
            case 'spike'
                fun = @(x) sum(abs(Nforward_discretized(x(1),[t(j) 0],samppars,consts,Nmu,scenario,Nlogical(i,:),x(2)) - iobs'-idobs'));
        end
        % minimize isfit
        [sol,~,exitflag] = fminsearchbnd(fun,x0,LB,UB,options);  

        if exitflag == 0 
            sol = sol*nan;  % set values of failed optimizations to nan
        end

        p1up(i,j) = sol(1); p2up(i,j) = sol(2);
        
    end
end


% clean up output (necessary because crazy input values can produce crazy
% solutions)
p2(p2 <0) = nan;     % set obviously wrong solutions to nan
p1(p1 <0) = nan;    
p2(isinf(p2)) = nan;     % set obviously wrong solutions to nan
p1(isinf(p1)) = nan;   
p1up(p1up   <0) = nan; 
p2up(p2up   <0) = nan; 
p1low(p1low <0) = nan; 
p1low(p1low <0) = nan; 
p1up(isinf(p1up)) = nan; 
p2up(isinf(p2up)) = nan; 
p1low(isinf(p1low)) = nan; 
p1low(isinf(p1low)) = nan; 

end