function out = retention_plot(n_name1,n_name2,ages,plotFlag)

% This code makes He-3 retention diagrams from data structure returned
% from get_ages_v3. One of the nuclides must be He-3; if not, results will
% probably not make any sense. 
%
% Returns a filename with no extension identifying the plot files. 
%
% Greg Balco - Berkeley Geochronology Center - February, 2015.
% Development, not licensed for use or distribution. 
% Changed to GMT5 only Dec 2018

% Load consts

load consts_v3;

% Default to MATLAB plotting not GMT plotting
if nargin < 4; plotFlag = 1; end

% Want He-3 to be the y-coordinate. 
index1 = find(strcmp(n_name1,consts.nuclides));
index2 = find(strcmp(n_name2,consts.nuclides));
if strcmp(n_name1,'N3quartz')
    x.index = index2; y.index = index1;
else
    x.index = index1; y.index = index2;
end
    
% Accumulate relevant apparent ages

% Note this has to accommodate multiple instances of a particular
% nuclide measurement. Here we choose to average. 

x.t = [];x.dt = [];
y.t = [];y.dt = [];

for a = 1:ages.numsamples
    % Extract relevant information. This is redundant. 
    this_sample.nuclide = ages.n.nuclide(ages.n.index == a);
    this_sample.t_St = ages.n.t_St(ages.n.index == a);
    this_sample.delt_ext_St = ages.n.delt_ext_St(ages.n.index == a);
    % Only do this if both present
    if any(strcmp(consts.nuclides{x.index},this_sample.nuclide)) && any(strcmp(consts.nuclides{y.index},this_sample.nuclide))
        % x nuclide, this can be lots of stuff
        x.ni = find(strcmp(consts.nuclides{x.index},this_sample.nuclide));
        if length(x.ni) > 1
            [x.t(end+1),x.dt(end+1)] = ewmean(this_sample.t_St(x.ni),this_sample.delt_ext_St(x.ni));
        else
            x.t(end+1) = this_sample.t_St(x.ni);
            x.dt(end+1) = this_sample.delt_ext_St(x.ni);
        end
        % y nuclide, this should always be He-3-quartz
        y.ni = find(strcmp(consts.nuclides{y.index},this_sample.nuclide));
        if length(y.ni) > 1
            [y.t(end+1),y.dt(end+1)] = ewmean(this_sample.t_St(y.ni),this_sample.delt_ext_St(y.ni));
        else
            y.t(end+1) = this_sample.t_St(y.ni);
            y.dt(end+1) = this_sample.delt_ext_St(y.ni);
        end
    end
    clear this_sample
end

% Now get axis limits. 
maxl = max([max(x.t + x.dt) max(y.t + y.dt)]);
minl = min([min(x.t - x.dt) min(y.t - y.dt)]);
if minl < 10; minl = 10; end

if log10(maxl)-log10(minl) > 1.5
    axx = 'log';
    xmin = 10^(floor(log10(minl)));
    if xmin < 10; xmin = 10; end 
    xmax = 10^(ceil(log10(maxl)));
else
    axx = 'linear';
    xmin = 0;
    xmax1 = maxl.*1.2;
    xsp = 10.^(floor(log10(maxl)));
    xmax = xsp.*ceil(xmax1./xsp);
end


% 1:1 line
xx1 = [xmin xmax]; yy1 = [1 1];
    

if plotFlag == 1
    % Case plotting in MATLAB, not making GMT file
    figure;
    
    hub = plot(xx1,yy1,'k');hold on;

    title([strtok(consts.properName{y.index}) ' / ' strtok(consts.properName{x.index})],'fontsize',14);
    set(gca,'ylim',[0 1.2]);
    xlabel(['Apparent ' strtok(consts.properName{x.index}) ' age'],'fontsize',14);
    ylabel(['Apparent ' strtok(consts.properName{y.index}) ' age'],'fontsize',14);
    set(gca,'fontsize',12); 
    set(gca,'xscale',axx,'xlim',[xmin xmax]);
    
    % Retention curves (from consts file)
    for a = 1:length(consts.RHe.TCs)
        xx = consts.RHe.tv; yy = consts.RHe.R(a,:);
        plot(xx,yy,'g');
    end
    
    % Ellipses
    if ~isempty(x.t)
        plot(x.t,y.t,'r.');
        for a = 1:length(x.t)
            ellipse(x.t(a),x.dt(a),y.t(a),y.dt(a),1,'r');
        end
    end
    
    out = [];
    
elseif plotFlag > 1
    % Case using GMT
    % Generate random number filename
    rand('seed',sum(100*clock));
    uniqueid = int2str(round(rand(1)*1e6));
    
    
    if consts.isLocal == 0
        % Assume running on Linux system
        % everything goes on in the scratch directory.
        fname = [consts.scratchdir 'hrplotv3_' uniqueid];
        gmtname = [fname '.gmt'];
        psname = [fname '.ps'];
        pngname = [fname '.png'];
    else
        % Perhaps running on GB's laptop
        fname = [consts.scratchdir 'hrplotv3_' uniqueid];
        gmtname = [fname '.gmt'];
        psname = [fname '.eps'];
        pngname = [fname '.png'];
    end
    
    % GMT format stuff
    formatstr = '--MAP_ANNOT_OFFSET_PRIMARY=0.1i --MAP_TICK_LENGTH_PRIMARY=-0.05i --MAP_FRAME_PEN=0.5p --FONT_ANNOT_PRIMARY=10p,2,black --FONT_LABEL=14p,3,black';

    
    % open file for write
    fid = fopen(gmtname,'w');

    % Create GMT script
    % Define format strings
    
    lxs = ['"Apparent ' strtok(consts.properName{x.index}) ' age (yr)"'];
    lys = ['"He-3 retention"'];
    
    if strcmp(axx,'log')
        jstring = '-JX5.5il/3.5i';
        rstring = ['-R' sprintf('%0.0f',xmin) '/' sprintf('%0.0f',xmax) '/0/1.2'];
        bstring = ['-Ba1f1g3p:' lxs ':/a0.2g0.2:' lys ':WeSn'];
    else
        jstring = '-JX5.5i/3.5i';
        rstring = ['-R' sprintf('%0.0f',xmin) '/' sprintf('%0.0f',xmax) '/0/1.2'];
        bstring = ['-Ba' sprintf('%0.0f',xsp) 'g' sprintf('%0.0f',xsp) ':' lxs ':/a0.2g0.2:' lys ':WeSn'];
    end
    
    % Tell shell path to GMT if necessary, this is for GB's Mac
    %     if consts.isLocal == 1;
    %         fprintf(fid,'%s\n\n','PATH=$PATH:/sw/bin');
    %         fprintf(fid,'%s\n\n','export PATH');
    %     end;
    
    % Make the GMT basemap
    fprintf(fid,'%s\n','# ----------- Draw the plot axes ------------');
    temp_string = ['gmt psbasemap ' bstring ' ' jstring ' ' rstring ' -K -P ' formatstr ' > ' psname];
    fprintf(fid,'%s\n',temp_string);
    
    % GMT simple-exposure boundaries
    fprintf(fid,'%s\n','# ----------- Plot the total-retention line ------------');
    temp_string = ['gmt psxy ' rstring ' ' jstring ' -P -K -O -W0.5p,0 << EOF >> ' psname];
    fprintf(fid,'%s\n',temp_string);
    fprintf(fid,'%0.5g %0.5g\n',[xx1' yy1']');
    fprintf(fid,'%s\n','EOF');
    
    % GMT retention curves
    fprintf(fid,'%s\n','# ----------- Plot the total-retention line ------------');
    
    for a = 1:length(consts.RHe.TCs)
        temp_string = ['gmt psxy ' rstring ' ' jstring ' -P -K -O -W0.5p,0 << EOF >> ' psname];
        fprintf(fid,'%s\n',temp_string);
        xx2 = consts.RHe.tv; yy2 = consts.RHe.R(a,:);
        fprintf(fid,'%0.5g %0.5g\n',[xx2' yy2']');
        fprintf(fid,'%s\n','EOF');
    end  
    
    % GMT ellipses
    for a = 1:length(x.t)
        % get the x,y for the ellipse
        [ex,ey] = ellipse(x.t(a),x.dt(a),y.t(a),y.dt(a),0);
        % File write 
        fprintf(fid,'%s\n','# ----------- Plot the ellipse ------------');
        temp_string = ['gmt psxy ' rstring ' ' jstring ' -P  -O -K -W0.5p,255/50/50 -G255/200/200 << EOF >> ' psname];
        fprintf(fid,'%s\n',temp_string);
        fprintf(fid,'%0.5g %0.5g\n',[ex',ey']');
        fprintf(fid,'%s\n','EOF');
    end
    
    % GMT center dots
    fprintf(fid,'%s\n','# ----------- Plot the center point of the ellipse ------------');
    temp_string = ['gmt psxy ' rstring ' ' jstring ' -P  -O -Sc0.05i -G255/50/50 << EOF >> ' psname];
    fprintf(fid,'%s\n',temp_string);
    % loop for all the center dots
    for a = 1:length(x.t)
        fprintf(fid,'%0.5g %0.5g\n',[x.t(a) y.t(a)./x.t(a)]');
    end
    fprintf(fid,'%s\n','EOF');

    % close the file
    fclose(fid);
    
    if consts.isLocal == 0
        % Then you are running on hess...
        % Execute the script you just made
        system(['chmod a+x ' gmtname]);
        system(gmtname);

        % convert the ps file to a png using ImageMagick...
        system(['convert ' psname ' -trim -transparent white ' pngname]);   
    else
        % On a Mac
        % Just send the command to the command window for cut and paste.
        % Too much hassle with inconsistent GMT/MATLAB libs. 
        disp(['sh ' gmtname]);
        disp(['open ' psname]);
    end
   
    % Return file name stub
    out = ['hrplotv3_' uniqueid];
    
end
