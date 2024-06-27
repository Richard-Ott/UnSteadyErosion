function out = summarize_cal_data(ages)

% Extract needed data

% Associate sites with individual production rate estimates

sitenames = ages.c.site_name(ages.n.index);
usites = unique(sitenames);
numsites = length(usites);
sitenums = zeros(size(sitenames));

for a = 1:numsites;
    sitenums = sitenums + a.*strcmp(sitenames,usites{a});
end;

% Extract data for scaling scheme and pass to averaging subroutine
% This needs to be fixed to only act on scaling schemes that actually
% exist. An alternative would be just to turn off the short/long flag. 

sfs = {'St','Lm','LSDn'};

for a = 1:length(sfs);
    eval(['in.P = ages.n.calc_P_' sfs{a} ';']);
    eval(['in.delP = ages.n.calc_delP_' sfs{a} ';']);
    eval(['in.minP = ages.n.calc_minP_' sfs{a} ';']);
    eval(['in.delminP = ages.n.calc_delminP_' sfs{a} ';']);
    eval(['in.maxP = ages.n.calc_maxP_' sfs{a} ';']);
    eval(['in.delmaxP = ages.n.calc_delmaxP_' sfs{a} ';']);
    in.sitenums = sitenums;
    in.sfs = sfs{a};
    
    this_summary = do_summary(in);
    
    eval(['out.' sfs{a} ' = this_summary;']);    
end;
    
% ------------------------ end main fctn -------------------------

function s = do_summary(in)

% Subroutine to do averaging for one scaling scheme

% Determine which min/max/exact exist

site_exact = unique(in.sitenums(in.P > 0));
site_max = unique(in.sitenums(in.maxP > 0));
site_min = unique(in.sitenums(in.minP > 0));

% Now compute means/medians/etc by site

sitemean = [];
sitemedian = [];
sitesd = [];
max_sitemean = [];
max_sitemedian = [];
max_sitesd = [];
min_sitemean = [];
min_sitemedian = [];
min_sitesd = [];

for a = site_exact';
    thisP = in.P(in.sitenums == a);
    thisdP = in.delP(in.sitenums == a);
    sitemean(end+1) = mean(thisP);
    sitemedian(end+1) = median(thisP);
    if length(thisP) > 1;
        % If n > 1, take SD
        sitesd(end+1) = std(thisP);
    else
        % If n = 1, just use individual uncert
        sitesd(end+1) = thisdP(1);
    end;        
end;

for a = site_max';
    thisP = in.maxP(in.sitenums == a);
    thisdP = in.delmaxP(in.sitenums == a);
    max_sitemean(end+1) = mean(thisP);
    max_sitemedian(end+1) = median(thisP);
    if length(thisP) > 1;
        % If n > 1, take SD
        max_sitesd(end+1) = std(thisP);
    else
        % If n = 1, just use individual uncert
        max_sitesd(end+1) = thisdP(1);
    end; 
end;

for a = site_min';
    thisP = in.minP(in.sitenums == a);
    thisdP = in.delminP(site_min == a);
    min_sitemean(end+1) = mean(thisP);
    min_sitemedian(end+1) = median(thisP);
    if length(thisP) > 1;
        % If n > 1, take SD
        min_sitesd(end+1) = std(thisP);
    else
        % If n = 1, just use individual uncert
        min_sitesd(end+1) = thisdP(1);
    end; 
end;

% OK. Now we want to average the site averages with the following
% algorithm. One, take the average of the exact measurements. Two, if this
% average violates any of the max/min constraints, add that constraint and
% take the average again. Repeat until no new constraints are added. 

to_avg = sitemean;
to_avg_d = sitesd;
possmins = min_sitemean;
possmins_sites = site_min;
possmins_d = min_sitesd;
possmaxs = max_sitemean;
possmaxs_sites = site_max;
possmaxs_d = max_sitesd;

active_min = [];
active_max = [];

its = 0;
while 1;
    if its > 10;
        % Deal with unintended consequences
        break;
    end;
    s.avgP = mean(to_avg);
    mins = find(possmins >= s.avgP); % If avg is less than any min limits
    maxs = find(possmaxs <= s.avgP); % If avg is greater than any max limits
    if isempty(mins) && isempty(maxs);
        break;
    end;
    % Add active constraints to data to be averaged
    to_avg = [to_avg possmins(mins) possmaxs(maxs)];
    to_avg_d = [to_avg_d possmins_d(mins) possmaxs_d(maxs)];
    % Record which were activated
    active_min = [active_min possmins_sites(mins)];
    active_max = [active_max possmaxs_sites(maxs)];
    % Remove those from further consideration
    possmins = possmins(~mins);
    possmins_d = possmins_d(~mins);
    possmins_sites(~mins);
    possmaxs = possmaxs(~maxs);
    possmaxs_d = possmaxs(~maxs);
    possmaxs_sites(~maxs);
    its = its + 1;
end;

% We have now computed the average (s.avgP). 
% Next we want to compute the uncertainty. 
% We do this from scatter of individual production rate estimates around
% mean computed above. 

s.text = ['<b>' in.sfs ' scaling.</b> Fitting parameter: ' sprintf('%0.3f',s.avgP) '<br>'];

if length(in.P) > 1;
    % Case more than one measurement, calculate population stats
    % Obtain all active data: exact values and active constraints
    add_min_P = []; add_min_dP = [];
    for a = 1:length(active_min);
        add_min_P = [add_min_P in.minP(in.sitenums == active_min(a))'];
        add_min_dP = [add_min_dP in.delminP(in.sitenums == active_min(a))'];
    end;
    add_max_P = []; add_max_dP = [];
    for a = 1:length(active_max);
        add_max_P = [add_max_P in.maxP(in.sitenums == active_max(a))'];
        add_max_dP = [add_max_dP in.delmaxP(in.sitenums == active_max(a))'];
    end;
        
    activeP = [in.P(in.P > 0)' add_min_P add_max_P];
    activedP = [in.delP(in.P > 0)' add_min_dP add_max_dP];
    
    
    % Now calculate statistics on that.    
    s.delPtot = (std(activeP/s.avgP)); % Total scatter (relative units)
    s.delPsite = (std(to_avg./s.avgP)); % site scatter (relative units)
    % various population statistics
    s.chi2tot = sum(((activeP - s.avgP)./activedP).^2); % total. Note this doesn't include active constraints. 
    s.pchi2tot = chisquarecdf_rt(s.chi2tot,length(activeP)-1);
    s.chi2site = sum(((to_avg-s.avgP)./to_avg_d).^2); % by site
    s.pchi2site = chisquarecdf_rt(s.chi2site,length(to_avg)-1);
    s.text = [s.text 'Uncertainty by total scatter: ' sprintf('%0.3f',s.delPtot.*s.avgP) ' (' sprintf('%0.1f',s.delPtot.*100) '%%)<br>'];
    s.text = [s.text 'Chi-squared (all samples) is ' sprintf('%0.2f',s.chi2tot) ' for ' int2str(length(activeP)-1) ' DOF; p = ' sprintf('%0.3f',s.pchi2tot) '<br>'];
    if length(to_avg) > 1;
        s.text = [s.text 'Uncertainty by site-to-site scatter: ' sprintf('%0.3f',s.delPsite.*s.avgP) ' (' sprintf('%0.1f',s.delPsite.*100) '%%)<br>'];
        s.text = [s.text 'Chi-squared (all sites) is ' sprintf('%0.2f',s.chi2site) ' for ' int2str(length(to_avg)-1) ' DOF; p = ' sprintf('%0.3f',s.pchi2site) '<br>'];
    end;
    % Pick an uncertainty to carry through to subsequent calculations
    s.use_uncert = s.avgP.*max([s.delPtot s.delPsite]); 
else
    s.delPtot = in.delP(1);
    s.use_uncert = s.delPtot;
    s.delPsite = 0;
    s.chi2tot = 0;
    s.chi2site = 0;
    s.pchi2tot = 0;
    s.pchi2site = 0;
    s.text = [s.text 'One sample: uncertainty is ' sprintf('%0.2f',s.delPtot.*s.avgP)];
end;

    
    
    


