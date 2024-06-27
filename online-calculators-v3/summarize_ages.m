function out = summarize_ages(ages,camelplotFlag)

% This script calculates summary statistics on an array of exposure ages
% returned by get_ages_v3.m. Returns an output struct with summary data for
% (1) all nuclides, all samples together, as well as (2) each nuclide
% separately. Calculates (1) and (2) for all scaling methods (St, Lm,
% LSDn). 
%
% out = summarize_ages(ages,camelplotFlag)
% 
% second arg is optional and defaults to 0. If 1, generates camelplot(s)
% and returns filenames. 
% 
% Note that this does not distiguish between ages on different samples and
% multiple ages on same sample. All are treated as separate ages. 
% 
% This includes an outlier-detection scheme that prunes data based on a
% chi-squared criterion. Basically, it will prune outliers until the
% probability that the remanining data belong to a single population is >
% 0.01. Then it will calculate either the mean and SD, or the
% error-weighted mean, of the remaining data depending on what the
% chi-squared of the remanining data actually is. 
%
% Calculation of the external uncertainty on the summary age is a bit
% fudged in that it is linearized. Thus, this is inappropriate for ages
% that are old in relation to the half-life of the nuclide in question. In
% practice this probably should't be much of a problem. 
% 
% Optionally also calls return_camelplot.m.  
%
% Greg Balco
% Berkeley Geochronology Center
% December, 2016.

% Check to make sure not being sent calibration results
if ~isfield(ages.n,'t_St');
    error('summarize_ages.m: calibration results incorrectly passed');
end;

if nargin < 2;
    camelplotFlag = 0;
end;

% Load consts. This is only needed to define possible nuclides and file 
% paths, so there would be no reason to pass it from elsewhere. 
load consts_v3;

% Operate on various arrays of ages

% All ages together
out.all.St = sumstats(ages.n.t_St,ages.n.delt_int_St,ages.n.delt_ext_St);
if isfield(ages.n,'t_Lm');
    out.all.Lm = sumstats(ages.n.t_Lm,ages.n.delt_int_Lm,ages.n.delt_ext_Lm);
else
    out.all.Lm = [];
end;
if isfield(ages.n,'t_LSDn');
    out.all.LSDn = sumstats(ages.n.t_LSDn,ages.n.delt_int_LSDn,ages.n.delt_ext_LSDn);
else
    out.all.LSDn = [];
end;

% Ages by nuclide and scaling scheme
if isfield(ages.n,'nindex'); % This turns off this block for Cl-36 ages. Assumption is for Cl-36 there are only Cl-36 ages. 
    for a = 1:length(consts.nuclides);
        if any(a == ages.n.nindex);
            use = find(ages.n.nindex == a);
            thiss = sumstats(ages.n.t_St(use),ages.n.delt_int_St(use),ages.n.delt_ext_St(use));
            thiss.use = use;
            eval(['out.' consts.nuclides{a} '.St = thiss;']);
            if isfield(ages.n,'t_Lm');
                thiss = sumstats(ages.n.t_Lm(use),ages.n.delt_int_Lm(use),ages.n.delt_ext_Lm(use));
                thiss.use = use;
                eval(['out.' consts.nuclides{a} '.Lm = thiss;']);
            else
                eval(['out.' consts.nuclides{a} '.Lm = [];']);
            end;
            if isfield(ages.n,'t_LSDn');
                thiss = sumstats(ages.n.t_LSDn(use),ages.n.delt_int_LSDn(use),ages.n.delt_ext_LSDn(use));
                thiss.use = use;
                eval(['out.' consts.nuclides{a} '.LSDn = thiss;']);
            else
                eval(['out.' consts.nuclides{a} '.LSDn = [];']);
            end;
        else
            eval(['out.' consts.nuclides{a} '.St = [];']);
            eval(['out.' consts.nuclides{a} '.Lm = [];']);
            eval(['out.' consts.nuclides{a} '.LSDn = [];']);
        end;
    end;  
end

% Now do camelplots.

if camelplotFlag == 1;

    if isfield(ages.n,'nindex')
        % Normal case
        numnuclides = length(unique(ages.n.nindex));
    else
        % Cl-36 case
        numnuclides = 1;
    end

    nloop = {'all'};
    % Enable below to make camelplots for all nuclide/SF pairs
    % if numnuclides > 1;
    %     nloop = [{'all'} consts.nuclides];
    % end;

    out.camels = {};

    for sf = {'St','Lm','LSDn'};
        for n = nloop;
            eval(['thiss = out.' char(n) '.' char(sf) ';']);
            if ~isempty(thiss);
                % Get correct t,dt
                eval(['thist = ages.n.t_' char(sf) ';']);
                eval(['thisdt = ages.n.delt_int_' char(sf) ';']);
                % Append what it is for labeling
                if strcmp(char(n),'all');
                    thiss.label = ['All data: ' char(sf) ' scaling'];
                else
                    thiss.label = [char(n) ': ' char(sf) ' scaling'];
                end;
                out.camels{end+1} = return_camelplot(thist,thisdt,thiss);
            end;
        end;
    end;

    
end;
% ----------------

function ss = sumstats(t,dtint,dtext)

% this subroutine computes various stats from vectors of ages and uncerts

% Remove any zeros from t. This could occur if saturated samples are
% passed, which probably shouldn't happen, but might. 

nonzero = find(t > 0);
t = t(nonzero);
dtint = dtint(nonzero);
dtext = dtext(nonzero);
if isempty(t);
    ss = []; return;
end;

% Start by calculating basic stats 

ss.mean = mean(t);
ss.sd = std(t);
[ss.ewmean,ss.ewse] = ewmean(t,dtint);
ss.chi2 = sum(((t - ss.ewmean)./dtint).^2);
ss.dof = length(t)-1;
ss.pchi2 = chisquarecdf_rt(ss.chi2,ss.dof);

ss.text = ['All data: mean ' sprintf('%0.0f',ss.mean) '; SD ' sprintf('%0.0f',ss.sd) '; chi-squared p-value is ' sprintf('%0.4f',ss.pchi2) '.<br>'];

% Outlier-stripping scheme based on chi-squared probability of belonging to
% a single population. 

ss.strip.ok = 1:length(t);
its = 0;

while 1;
    % Don't reduce data set too much...
    if its > length(t)./2;
        ss.text = [ss.text 'Less than half the data set remaining; stopped pruning. <br>'];
        break;
    end;
    % Determine chi-squared with respect to mean
    ss.strip.mean = mean(t(ss.strip.ok));
    ss.strip.sd = std(t(ss.strip.ok));
    [ss.strip.ewmean,ss.strip.ewse] = ewmean(t(ss.strip.ok),dtint(ss.strip.ok));
    ss.strip.allchi2 = ((t(ss.strip.ok) - ss.strip.mean)./dtint(ss.strip.ok)).^2;
    ss.strip.chi2 = sum(ss.strip.allchi2);
    ss.strip.dof = length(t(ss.strip.ok))-1;
    ss.strip.pchi2 = chisquarecdf_rt(ss.strip.chi2,ss.strip.dof);
    % Check chi-squared
    if ss.strip.pchi2 > 0.01;
        % Good enough; no outliers; stop
        break;
    end;
    % Don't actually do any pruning if you have two or less data.
    if length(ss.strip.ok) < 3;
        ss.text = [ss.text 'Two or less samples left; stopped pruning. <br>'];
        break;
    end;
    
    % If we got this far, chi-squared p is terrible. Remove sample with 
    % worst chi-squared and repeat. 
    worst = ss.strip.ok(find(ss.strip.allchi2 == max(ss.strip.allchi2)));
    ss.strip.ok = ss.strip.ok(ss.strip.ok ~= worst);
    its = its + 1;
end;

if its == 1;
    ss.text = [ss.text 'Pruned ' int2str(its) ' outlier. <br>'];
else
    ss.text = [ss.text 'Pruned ' int2str(its) ' outliers. <br>'];
end;

if ss.strip.pchi2 > 0.05;
    ss.sumval = ss.strip.ewmean;
    ss.sumdel_int = ss.strip.ewse;
    ss.text = [ss.text 'Remaining data have p greater than 0.05; using error-weighted mean. <br>'];
else
    ss.sumval = ss.strip.mean;
    ss.sumdel_int = ss.strip.sd;
    ss.text = [ss.text 'Remaining data have p less than 0.05; using mean and SD. <br>'];
end;
% Infer production rate uncertainty...note linearizes, i.e. assumes
% radioactive decay has negligible effect. Thus, this is appropriate only
% for ages that are young relative to the relevant half-life. 
reldelP = mean(sqrt(dtext(ss.strip.ok).^2 - dtint(ss.strip.ok).^2)./t(ss.strip.ok)); % All should give the same result, 
% but take mean anyway to be slightly error-tolerant
% Propagate into summary value
ss.sumdel_ext = sqrt((ss.sumdel_int).^2 + (reldelP.*ss.sumval).^2);

ss.text = [ss.text ' Summary value is ' sprintf('%0.0f',ss.sumval)];
ss.text = [ss.text ' +/- ' sprintf('%0.0f',ss.sumdel_int) ' (' sprintf('%0.0f',ss.sumdel_ext) ').<br>'];

% Decide whether your moraine belongs to the YD or not. 
% Normal cdf with sumval and external uncertainty
cdfpts = 0.5*(1 + erf(([14700 13000 12900 11700]-ss.sumval)./(sqrt(2).*ss.sumdel_ext)));
ss.prob_ACR = cdfpts(1) - cdfpts(2); ss.prob_YD = cdfpts(3) - cdfpts(4);
if ss.prob_ACR > 0.01 || ss.prob_YD > 0.01;
    ss.text = [ss.text 'If this is a moraine, the probability that it is Younger Dryas age is ' sprintf('%0.2f',ss.prob_YD) '.<br>'];
    ss.text = [ss.text 'The probability it belongs to the Antarctic Cold Reversal is ' sprintf('%0.2f',ss.prob_ACR) '.'];
end;










    
    
    
    
    
    
    
    
        