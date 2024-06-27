function out = generic_splot(n_name1,n_name2,ages,plotFlag)

% This code makes generic two-nuclide diagrams from data structure returned
% from get_ages_v3. This one doesn't do He-3 retention plots, that is a
% separate function, so doesn't expect the upstream script to pass it He-3
% concentrations. 
%
% Returns a filename with no extension identifying the plot files. 
%
% Greg Balco - Berkeley Geochronology Center - February, 2015.
% Revised Dec 2018 to use GMT5 only

% Load consts

load consts_v3;

% Default to MATLAB plotting not GMT plotting
if nargin < 4; plotFlag = 1; end

% Get in order based on order defined in constants file

index1 = find(strcmp(n_name1,consts.nuclides));
index2 = find(strcmp(n_name2,consts.nuclides));

order1 = consts.splot_order(index1);
order2 = consts.splot_order(index2);

% Want numerator to be the one with the lower order. 

if order1 < order2
    x.index = index2; y.index = index1;
else
    x.index = index1; y.index = index2;
end

% Make banana...

L = Lsp(); % ??

x.l = consts.l(x.index); y.l = consts.l(y.index);

tt = logspace(0,log10(25e6),100);
ee = logspace(-7,1,100);

if x.l > 0;
    x.nt = (1./x.l).*(1-exp(-x.l.*tt));
else
    x.nt = tt;
end
x.ne = (1./(x.l + ee./L));

if y.l > 0
    y.nt = (1./y.l).*(1-exp(-y.l.*tt));
else
    y.nt = tt;
end
y.ne = (1./(y.l + ee./L));

xx1 = x.nt; yy1 = y.nt./x.nt;
xx2 = x.ne; yy2 = y.ne./x.ne;
    
% Accumulate relevant normalized nuclide concentrations

% Allocate 

x.Nnorm_St = [];x.delNnorm_St = [];
y.Nnorm_St = [];y.delNnorm_St = [];
x.Nnorm_Lm = [];x.delNnorm_Lm = [];
y.Nnorm_Lm = [];y.delNnorm_Lm = [];
x.Nnorm_LSDn = [];x.delNnorm_LSDn = [];
y.Nnorm_LSDn = [];y.delNnorm_LSDn = [];

for a = 1:ages.numsamples
    % Extract relevant information. This is redundant. 
    this_sample.nuclide = ages.n.nuclide(ages.n.index == a);
    this_sample.Nnorm_St = ages.n.Nnorm_St(ages.n.index == a);
    this_sample.delNnorm_St = ages.n.delNnorm_St(ages.n.index == a);
    % Get other normalizations if present.
    isLm = 0; isLSDn = 0;
    if isfield(ages.n,'Nnorm_Lm')
        isLm = 1;
        this_sample.Nnorm_Lm = ages.n.Nnorm_Lm(ages.n.index == a);
        this_sample.delNnorm_Lm = ages.n.delNnorm_Lm(ages.n.index == a);
    end
    if isfield(ages.n,'Nnorm_LSDn')
        isLSDn = 1;
        this_sample.Nnorm_LSDn = ages.n.Nnorm_LSDn(ages.n.index == a);
        this_sample.delNnorm_LSDn = ages.n.delNnorm_LSDn(ages.n.index == a);
    end
    % Only do this if both present
    if any(strcmp(consts.nuclides{x.index},this_sample.nuclide)) && any(strcmp(consts.nuclides{y.index},this_sample.nuclide))
        % x nuclide
        x.ni = find(strcmp(consts.nuclides{x.index},this_sample.nuclide));
        if length(x.ni) > 1
            % If there are multiple measurements of the same thing, use the
            % error-weighted mean 
            [x.Nnorm_St(end+1),x.delNnorm_St(end+1)] = ewmean(this_sample.Nnorm_St(x.ni),this_sample.delNnorm_St(x.ni));
            if isLm == 1
                [x.Nnorm_Lm(end+1),x.delNnorm_Lm(end+1)] = ewmean(this_sample.Nnorm_Lm(x.ni),this_sample.delNnorm_Lm(x.ni));
            end
            if isLSDn == 1
                [x.Nnorm_LSDn(end+1),x.delNnorm_LSDn(end+1)] = ewmean(this_sample.Nnorm_LSDn(x.ni),this_sample.delNnorm_LSDn(x.ni));
            end
                
        else
            x.Nnorm_St(end+1) = this_sample.Nnorm_St(x.ni);
            x.delNnorm_St(end+1) = this_sample.delNnorm_St(x.ni);
            if isLm == 1
                x.Nnorm_Lm(end+1) = this_sample.Nnorm_Lm(x.ni);
                x.delNnorm_Lm(end+1) = this_sample.delNnorm_Lm(x.ni);
            end
            if isLSDn == 1
                x.Nnorm_LSDn(end+1) = this_sample.Nnorm_LSDn(x.ni);
                x.delNnorm_LSDn(end+1) = this_sample.delNnorm_LSDn(x.ni);
            end
        end
        % y nuclide
        y.ni = find(strcmp(consts.nuclides{y.index},this_sample.nuclide));
        if length(y.ni) > 1
            % Again, if multiple measurements of the same thing, use the
            % error-weighted mean
            [y.Nnorm_St(end+1),y.delNnorm_St(end+1)] = ewmean(this_sample.Nnorm_St(y.ni),this_sample.delNnorm_St(y.ni));
            if isLm == 1
                [y.Nnorm_Lm(end+1),y.delNnorm_Lm(end+1)] = ewmean(this_sample.Nnorm_Lm(y.ni),this_sample.delNnorm_Lm(y.ni));
            end
            if isLSDn == 1
                [y.Nnorm_LSDn(end+1),y.delNnorm_LSDn(end+1)] = ewmean(this_sample.Nnorm_LSDn(y.ni),this_sample.delNnorm_LSDn(y.ni));
            end
        else
            y.Nnorm_St(end+1) = this_sample.Nnorm_St(y.ni);
            y.delNnorm_St(end+1) = this_sample.delNnorm_St(y.ni);
            if isLm == 1
                y.Nnorm_Lm(end+1) = this_sample.Nnorm_Lm(y.ni);
                y.delNnorm_Lm(end+1) = this_sample.delNnorm_Lm(y.ni);
            end
            if isLSDn == 1
                y.Nnorm_LSDn(end+1) = this_sample.Nnorm_LSDn(y.ni);
                y.delNnorm_LSDn(end+1) = this_sample.delNnorm_LSDn(y.ni);
            end
        end
    end
    sample_names{a} = ages.s.sample_name{a};
    clear this_sample
end

% Now get axis limits. 
maxx = max(x.Nnorm_St + x.delNnorm_St);
minx = min(x.Nnorm_St - x.delNnorm_St);
if minx < 10; minx = 10; end

if (log10(maxx)-log10(minx) > 1.5)
    axx = 'log';
    xmin = 10^(floor(log10(minx)));
    if xmin < 10; xmin = 10; end 
    xmax = 10^(ceil(log10(maxx)));
else
    axx = 'linear';
    xmin = 0;
    xmax1 = maxx.*1.2;
    xsp = 10.^(floor(log10(xmax1)));
    xmax = xsp.*ceil(xmax1./xsp);
end

% Determine if we want to use log axes in y

yLogFlag = 0;
tempR = y.Nnorm_St./x.Nnorm_St;
if any(tempR < 0.2); yLogFlag = 1; miny = 10.^floor(log10(min(tempR))); end

if plotFlag == 1
    % Case plotting in MATLAB, not making GMT file
    figure;
    
    hub = plot(xx1,yy1,'k');hold on;
    hlb = plot(xx2,yy2,'k');
    title([strtok(consts.properName{y.index}) ' / ' strtok(consts.properName{x.index})],'fontsize',14);
    set(gca,'ylim',[0 1.2]);
    xlabel(['[' strtok(consts.properName{x.index}) ']^*'],'fontsize',14);
    ylabel(['[' strtok(consts.properName{y.index}) ']^* / [' strtok(consts.properName{x.index}) ']^*'],'fontsize',14);
    set(gca,'fontsize',12); 
    set(gca,'xscale',axx,'xlim',[xmin xmax]);
    if yLogFlag == 1; set(gca,'yscale','log','ylim',[miny 2]); end
    
    % Ellipses
    if ~isempty(x.Nnorm_St)
        plot(x.Nnorm_St,(y.Nnorm_St./x.Nnorm_St),'r.');
        for a = 1:length(x.Nnorm_St)
            ellipse(x.Nnorm_St(a),x.delNnorm_St(a),y.Nnorm_St(a),y.delNnorm_St(a),1,'r');
            % Also text
            text(x.Nnorm_St(a),(y.Nnorm_St(a)./x.Nnorm_St(a)),sample_names{a});
        end
        if isLm == 1
           for a = 1:length(x.Nnorm_Lm)
            ellipse(x.Nnorm_Lm(a),x.delNnorm_Lm(a),y.Nnorm_Lm(a),y.delNnorm_Lm(a),1,'g');
           end
        end
        if isLSDn == 1
           for a = 1:length(x.Nnorm_LSDn)
            ellipse(x.Nnorm_LSDn(a),x.delNnorm_LSDn(a),y.Nnorm_LSDn(a),y.delNnorm_LSDn(a),1,'b');
           end
        end
    end
    
    out = [];
    
elseif plotFlag > 1
    % Case using GMT
    
    % Generate random number filename
    rand('seed',sum(100*clock));
    uniqueid = int2str(round(rand(1)*1e6));
    
    
    if consts.isLocal == 0
        % Running on Linux machine
        % everything goes on in the scratch directory.
        fname = [consts.scratchdir 'splotv3_' uniqueid];
        gmtname = [fname '.gmt'];
        psname = [fname '.ps'];
        pngname = [fname '.png'];
    else
        % Perhaps running on GB's laptop, refer to consts file
        fname = [consts.scratchdir 'splotv3_' uniqueid];
        gmtname = [fname '.gmt'];
        psname = [fname '.eps'];
        pngname = [fname '.png'];
    end
    
    % format commands
    formatstr = '--MAP_ANNOT_OFFSET_PRIMARY=0.1i --MAP_TICK_LENGTH_PRIMARY=-0.05i --MAP_FRAME_PEN=0.5p --FONT_ANNOT_PRIMARY=10p,2,black --FONT_LABEL=14p,3,black';
    
    % open file for write
    fid = fopen(gmtname,'w');

    % Create GMT script
    % Define format strings
    
    lxs = ['"[' strtok(consts.properName{x.index}) ']*"'];
    lys = ['"[' strtok(consts.properName{y.index}) ']* / [' strtok(consts.properName{x.index}) ']*"'];
    
    if strcmp(axx,'log')
        if yLogFlag == 1
            jstring = '-JX5.5il/3.5il';
            rstring = ['-R' sprintf('%0.0f',xmin) '/' sprintf('%0.0f',xmax) '/' sprintf('%0.1e',miny) '/2'];
            bstring = ['-Ba1f1g3p:' lxs ':/a1f1g3p:' lys ':WeSn'];
        else
            jstring = '-JX5.5il/3.5i';
            rstring = ['-R' sprintf('%0.0f',xmin) '/' sprintf('%0.0f',xmax) '/0/1.2'];
            bstring = ['-Ba1f1g3p:' lxs ':/a0.2g0.2:' lys ':WeSn'];
        end
    else
        if yLogFlag == 1
            jstring = '-JX5.5i/3.5il';
            rstring = ['-R' sprintf('%0.0f',xmin) '/' sprintf('%0.0f',xmax) '/' sprintf('%0.1e',miny) '/2'];
            bstring = ['-Ba' sprintf('%0.0f',xsp) 'g' sprintf('%0.0f',xsp) ':' lxs ':/a1f1g3p:' lys ':WeSn'];
        else
            jstring = '-JX5.5i/3.5i';
            rstring = ['-R' sprintf('%0.0f',xmin) '/' sprintf('%0.0f',xmax) '/0/1.2'];
            bstring = ['-Ba' sprintf('%0.0f',xsp) 'g' sprintf('%0.0f',xsp) ':' lxs ':/a0.2g0.2:' lys ':WeSn'];
        end
    end
    
    % Tell shell path to GMT if necessary, this is for GB's Mac
         if consts.isLocal == 1;
             fprintf(fid,'%s\n\n','PATH=$PATH:/sw/bin');
             fprintf(fid,'%s\n\n','export PATH');
         end;
    
    % Make the GMT basemap
    fprintf(fid,'%s\n','# ----------- Draw the plot axes ------------');
    temp_string = ['gmt psbasemap ' bstring ' ' jstring ' ' rstring ' -K -P ' formatstr ' > ' psname];
    fprintf(fid,'%s\n',temp_string);
    
    % GMT simple-exposure boundaries
    fprintf(fid,'%s\n','# ----------- Plot the simple exposure region ------------');
    temp_string = ['gmt psxy ' rstring ' ' jstring ' -P -K -O -W0.5p,0 << EOF >> ' psname];

    fprintf(fid,'%s\n',temp_string);
    fprintf(fid,'%0.5g %0.5g\n',[xx1' yy1']');
    fprintf(fid,'%s\n','>');
    fprintf(fid,'%0.5g %0.5g\n',[xx2' yy2']');
    fprintf(fid,'%s\n','EOF');
    
    % GMT ellipses
    
    if ~isempty(x.Nnorm_Lm)
        for a = 1:length(x.Nnorm_Lm)
            % get the x,y for the ellipse
            [ex,ey] = ellipse(x.Nnorm_Lm(a),x.delNnorm_Lm(a),y.Nnorm_Lm(a),y.delNnorm_Lm(a),0);
            % File write 
            fprintf(fid,'%s\n','# ----------- Plot ellipse for Lm ------------');
            temp_string = ['gmt psxy ' rstring ' ' jstring ' -P  -O -K -W0.5p,100/150/100 -G230/255/230 << EOF >> ' psname];
            fprintf(fid,'%s\n',temp_string);
            fprintf(fid,'%0.5g %0.5g\n',[ex',ey']');
            fprintf(fid,'%s\n','EOF');
        end
    end
    

    
    if ~isempty(x.Nnorm_LSDn)
        for a = 1:length(x.Nnorm_LSDn)
            % get the x,y for the ellipse
            [ex,ey] = ellipse(x.Nnorm_LSDn(a),x.delNnorm_LSDn(a),y.Nnorm_LSDn(a),y.delNnorm_LSDn(a),0);
            % File write 
            fprintf(fid,'%s\n','# ----------- Plot ellipse for LSDn ------------');
            temp_string = ['gmt psxy ' rstring ' ' jstring ' -P  -O -K -W0.5p,100/100/255 -G230/230/255 << EOF >> ' psname];
           
            fprintf(fid,'%s\n',temp_string);
            fprintf(fid,'%0.5g %0.5g\n',[ex',ey']');
            fprintf(fid,'%s\n','EOF');
        end
    end
    
    
        
    for a = 1:length(x.Nnorm_St)
        % get the x,y for the ellipse
        [ex,ey] = ellipse(x.Nnorm_St(a),x.delNnorm_St(a),y.Nnorm_St(a),y.delNnorm_St(a),0);
        % File write 
        fprintf(fid,'%s\n','# ----------- Plot ellipse for St ------------');
        temp_string = ['gmt psxy ' rstring ' ' jstring ' -P  -O -K -W0.5p,255/50/50 -G255/200/200 << EOF >> ' psname]; 
        fprintf(fid,'%s\n',temp_string);
        fprintf(fid,'%0.5g %0.5g\n',[ex',ey']');
        fprintf(fid,'%s\n','EOF');
    end
    
    % GMT center dots, Lm
    fprintf(fid,'%s\n','# ----------- Plot the center points of the ellipses for Lm ------------');
    temp_string = ['gmt psxy ' rstring ' ' jstring ' -P  -O -K -Sc0.05i -W0.5p,100/150/100 << EOF >> ' psname];
    fprintf(fid,'%s\n',temp_string);
    % loop for all the center dots
    for a = 1:length(x.Nnorm_Lm)
        fprintf(fid,'%0.5g %0.5g\n',[x.Nnorm_Lm(a) (y.Nnorm_Lm(a)./x.Nnorm_Lm(a))]');
    end
    fprintf(fid,'%s\n','EOF');
    
    % GMT center dots, LSDn
    fprintf(fid,'%s\n','# ----------- Plot the center points of the ellipses for LSDn ------------');
    temp_string = ['gmt psxy ' rstring ' ' jstring ' -P  -O -K -Sc0.05i -W0.5p,100/100/255 << EOF >> ' psname];
    fprintf(fid,'%s\n',temp_string);
    % loop for all the center dots
    for a = 1:length(x.Nnorm_LSDn)
        fprintf(fid,'%0.5g %0.5g\n',[x.Nnorm_LSDn(a) (y.Nnorm_LSDn(a)./x.Nnorm_LSDn(a))]');
    end
    fprintf(fid,'%s\n','EOF');
    
    % GMT center dots, St
    fprintf(fid,'%s\n','# ----------- Plot the center points of the ellipses for St ------------');
    temp_string = ['gmt psxy ' rstring ' ' jstring ' -P  -O -Sc0.05i -G255/50/50 << EOF >> ' psname];
    fprintf(fid,'%s\n',temp_string);
    % loop for all the center dots
    for a = 1:length(x.Nnorm_St)
        fprintf(fid,'%0.5g %0.5g\n',[x.Nnorm_St(a) (y.Nnorm_St(a)./x.Nnorm_St(a))]');
    end
    fprintf(fid,'%s\n','EOF');
    
    
    

    % close the file
    fclose(fid);
    
    if consts.isLocal == 0
        % Then you are running on a linux machine...
        % Execute the script you just made
        system(['chmod a+x ' gmtname]);
        system(gmtname);

        % convert the ps file to a png using ImageMagick...
        system(['convert ' psname ' -trim -transparent white ' pngname]);       
    else
        % On GB's Mac
        % Just send the command to the command window for cut and paste.
        % Too much hassle with inconsistent GMT/MATLAB libs. 
        disp(['sh ' gmtname]);
        disp(['open ' psname]);
    end
   
    % Return file name stub
    out = ['splotv3_' uniqueid];
end
