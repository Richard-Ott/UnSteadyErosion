function out = get_DipRc(in,Rcfname,Sfname,plotFlag)

% This function extracts cutoff rigidity time series for a set of
% exposure-dating samples. Uses dipole approximation with paleopole
% positions fit to Pavon-Carrasco 13ka reconstruction per Nat Lifton for
% 0-13 ka, then Nat's 'GLOPAD' compilation to 2 Ma, then the mean of that
% for any time prior to 2 Ma. 
% 
% Also spits out solar variability parameter a la Lifton. 
%
% Syntax: out = get_DipRc(in,fname,sfname,plotFlag);
%
% Input arg in must have
%   in.lat sample latitudes (DD)
%   in.long sample longitude (DD); can be in -180-180 or 0-360
%   in.t length of Rc record required for each sample
%   in.yr date of collection of each sample
% These must all be nx1 vectors. 
% 
% Input arg fname is the name of a mat-file with the paleomag data.
%
% Input arg sfname is the name of a mat-file with the solar variability
% data.
%
% plotFlag is optional and defaults to zero (no plotting). Enter 1 for
% diagnostic plots. 
%
% Output has
%   out.tmin and out.tmax cell array (one cell for each sample) of time 
%   vectors with bounds of each time step.
%   out.Rc cell array of Rc for each sample computed from paleopole-
%       dipole approximation
%   out.S solar parameter a la Lifton, on same time vector
% 
% Note that the time vectors in out.tmin, out.tmax have units of years
% before sample
% collection, not absolute ages. So the magnetic field reconstructions
% correspond to the correct year. The aim here is to eventually compute
% exposure durations, not absolute ages. 
% 
% Note assumption is magnetic field is piecewise constant during time steps
%
% Greg Balco
% Berkeley Geochronology Center 
% October, 2015
% Modified Dec 2016 to use time step bounds, not middle of time step

if nargin < 4; plotFlag = 0; end;

% Rectify longitudes
in.long(in.long < 0) = in.long(in.long < 0)+360;
% Convert to radians 
latr = d2r(in.lat);
longr = d2r(in.long);

% OK, what is supplied here is an exposure duration for which we want an Rc
% vector. The exposure duration ends at the year of sample collection. The
% paleomag data are supplied by Nat in years before a zero year which is 
% defined in the var tzero.

% Load data

load(Sfname);
load(Rcfname);

% Define polynomial in cos(lat) from LSD2014, Equation 2
% This generates the geocentric dipole effective cutoff rigidity
pRc = [-448.004 1189.18 -1152.15 522.061 -103.241 6.89901 0];

% Loop
for a = 1:length(in.lat);  
    % Determine geomagnetic latitude at time t for short record
    tempGMLat = abs((pi/2)-angdistr(latr(a),longr(a),m.latr_pp,m.lonr_pp)); % Note in radians
    % Note already took absolute value above
    % Apply polynomial to compute Rc for short record
    tempDipRc = m.MM01.*polyval(pRc,cos(tempGMLat)); 
    
    % Convert real years to sample exposure years
    temptmin = m.t1min + (in.yr(a) - m.tzero);
    temptmax = m.t1max + (in.yr(a) - m.tzero);
    
    % Note must also set youngest tmin in old record accordingly. We don't
    % adjust the others, so given Nat's time spacing, this code will fail
    % in 194 years. Note this resets this value in each loop, so cannot
    % rely on this value in subsequent loops. I think this is OK. 
    
    m.t2min(1) = temptmax(end);
    
    % Adjust short record for before and after-zero-year cases. 
    if in.yr(a) > m.tzero;
        % After-zero-year case assumes unchanging magnetic field between
        % zero year of field reconstructions and year of sampling. 
        % This just involves changing the minimum age of the first time
        % step to zero.
        temptmin(1) = 0;
        tempS = SPhi; % Also assign solar var parameter
    elseif in.yr(a) < m.tzero;
        % Before-zero-year case needs to remove points from front end.
        inpast = find(temptmax > 0);
        temptmin = temptmin(inpast);
        temptmax = temptmax(inpast);
        tempDipRc = tempDipRc(inpast);
        tempS = SPhi(inpast); % Note assumes Sphi has same numel as MM01
        % Now set end of first time step to zero
        temptmin(1) = 0;
    else
        % sampled in zero year; do nothing. 
        % tempDipRc already assigned. 
        tempS = SPhi;
    end;
    
    % Now deal with different length cases. 
    if in.t(a) <= m.t2min(1);
        % Case all time before long record starts, that is, uses only short
        % record. 
        ok = find(temptmin < in.t(a)); 
        % The above gives all time steps whose minimum age is less than the
        % requested t. So the last one of these has a max age
        % greater than the requested t, i.e. it overshoots. 
        out.tmin{a} = temptmin(ok);
        out.tmax{a} = temptmax(ok);
        out.tmax{a}(end) = in.t(a); % Last tmax is equal to requested t
        out.Rc{a} = tempDipRc(ok);
        out.S{a} = tempS(ok);           
    elseif in.t(a) <= m.t2(end);       
        % Case in long record. 
        ok = find(m.t2min < in.t(a));
        % The above gives elts in t2 with min less than requested age
        out.tmin{a} = [temptmin m.t2min(ok)]; % Add to short record
        out.tmax{a} = [temptmax m.t2max(ok)]; % Add to short record
        out.tmax{a}(end) = in.t(a); % Last tmax equals requested t
        out.S{a} = [tempS zeros(size(ok))+SPhiInf]; % Add mean val to short record
        out.Rc{a} = [tempDipRc m.MM02(ok).*polyval(pRc,cos(abs(latr(a))))]; % Long record is GAD       
    else
        % Case longer than long record. Add a single time step onto the end
        % with mean MM0 value for long record. 
        % This needs to be fixed to add a log-spaced series of time steps. 
        out.tmin{a} = [temptmin m.t2min];
        out.tmax{a} = [temptmax m.t2max];
        out.S{a} = [tempS zeros(size(m.MM02))+SPhiInf];
        out.Rc{a} = [tempDipRc m.MM02.*polyval(pRc,cos(abs(latr(a))))]; % Long record is GAD   
        % Now just change oldest time step to have tmax = in.t(a)
        out.tmax{a}(end) = in.t(a); 
    end;
    
    % Deal with rigidities that leak below zero d/t polynomial fit
    out.Rc{a} = abs(out.Rc{a});
    
end;
        
%% Diagnostic plotting

if plotFlag == 1;
    figure;
    subplot(2,1,1);
    for a = 1:length(in.lat);
        thistmin = out.tmin{a};
        thistmax = out.tmax{a};
        thisRc = out.Rc{a};
        for b = 1:length(thistmin);
            xx = [thistmin(b) thistmax(b)];
            yy = [thisRc(b) thisRc(b)];
            plot(xx,yy,'b'); hold on;
        end;
        plot([in.t(a) in.t(a)],[0 16],'k');
    end;
    
    subplot(2,1,2);
    for a = 1:length(in.lat);
        thistmin = out.tmin{a};
        thistmax = out.tmax{a};
        thisS = out.S{a};
        for b = 1:length(thistmin);
            xx = [thistmin(b) thistmax(b)];
            yy = [thisS(b) thisS(b)];
            plot(xx,yy,'b'); hold on;
        end;
        plot([in.t(a) in.t(a)],[0 800],'k');
    end;    
end;


        