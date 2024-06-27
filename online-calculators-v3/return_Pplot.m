function out = return_Pplot(ages,summary)

% This code makes a plot summarizing production rate calibration results. 
% Called from within cal_input_v3.  
%
% args are the output from  get_ages_v3 (that contains all the individual
% sample results) and the output from summarize_cal_data (that contains
% summary info). 
%
% Returns a filename with no extension identifying the plot files. 
%
% Greg Balco - Berkeley Geochronology Center - December, 2016.
% Development, not licensed for use or distribution. 
%
    
% Load consts for environment info
load consts_v3

% ------------ Determine axis limits -------------------

% We are plotting data against elevation, because that usually looks good. 

% Obtain data. Remember we have min, max, and exact data. 

% Sort these 
isexact = ages.n.calc_P_St > 0;
ismin = ages.n.calc_minP_St > 0;
ismax = ages.n.calc_maxP_St > 0;

% Obtain elevation data
elv_exact = ages.s.elv(ages.n.index(isexact));
elv_min = ages.s.elv(ages.n.index(ismin));
elv_max = ages.s.elv(ages.n.index(ismax));

% Determine y axis limits
plot_emax = ceil(max([elv_exact' elv_min' elv_max'])./500).*500;
emax_str = int2str(plot_emax);
if plot_emax <= 500 
    ediv_str = int2str(100);
elseif plot_emax <= 1000
    ediv_str = int2str(200);
else
    ediv_str = int2str(500);
end

% Obtain residuals data

r_St = ages.n.calc_P_St(isexact)./summary.St.avgP;
dr_St = ages.n.calc_delP_St(isexact)./summary.St.avgP;
minr_St = ages.n.calc_minP_St(ismin)./summary.St.avgP;
mindr_St = ages.n.calc_delminP_St(ismin)./summary.St.avgP;
maxr_St = ages.n.calc_maxP_St(ismax)./summary.St.avgP;
maxdr_St = ages.n.calc_delmaxP_St(ismax)./summary.St.avgP;


r_Lm = []; dr_Lm = [];
if isfield(summary,'Lm')
    r_Lm = ages.n.calc_P_Lm(isexact)./summary.Lm.avgP;
    dr_Lm = ages.n.calc_delP_Lm(isexact)./summary.Lm.avgP;
    minr_Lm = ages.n.calc_minP_Lm(ismin)./summary.Lm.avgP;
    mindr_Lm = ages.n.calc_delminP_Lm(ismin)./summary.Lm.avgP;
    maxr_Lm = ages.n.calc_maxP_Lm(ismax)./summary.Lm.avgP;
    maxdr_Lm = ages.n.calc_delmaxP_Lm(ismax)./summary.Lm.avgP;
end

r_LSDn = []; dr_LSDn = [];
if isfield(summary,'LSDn')
    r_LSDn = ages.n.calc_P_LSDn(isexact)./summary.LSDn.avgP;
    dr_LSDn = ages.n.calc_delP_LSDn(isexact)./summary.LSDn.avgP;
    minr_LSDn = ages.n.calc_minP_LSDn(ismin)./summary.LSDn.avgP;
    mindr_LSDn = ages.n.calc_delminP_LSDn(ismin)./summary.LSDn.avgP;
    maxr_LSDn = ages.n.calc_maxP_LSDn(ismax)./summary.LSDn.avgP;
    maxdr_LSDn = ages.n.calc_delmaxP_LSDn(ismax)./summary.LSDn.avgP;
end


% Determine miss limits

temp1 = max([(abs([r_St' r_Lm' r_LSDn']) - 1) (abs([minr_St' minr_Lm' minr_LSDn']) - 1) (abs([minr_St' minr_Lm' minr_LSDn']) - 1)]);

mmax = (1 + ceil(temp1*10)./10)+0.05;
if mmax < 1.25; mmax=1.25; end
mmin = (1 - (mmax-1));
mmax_str = sprintf('%0.2f',mmax);
mmin_str = sprintf('%0.2f',mmin);

str_mmax = num2str(mmax); str_mmin = num2str(mmin);

% ---------- end sort out bounds -----------------------

%% start GMT

% Using GMT on web server

pstr = 'gmt ';
wsep = ',';
mstr = '';
formatstr = '--MAP_ANNOT_OFFSET_PRIMARY=0.1i --MAP_TICK_LENGTH_PRIMARY=-0.05i --FONT_ANNOT_PRIMARY=10p,2,black --FONT=12p,2,black';

% Generate random number filename

uniqueid = int2str(round(rand(1)*1e6));

fnamestub = ['Pplotv3_' uniqueid];

if consts.isLocal == 0
    % Assume running on hess.ess or equivalent
    % everything goes on in the scratch directory.
    fname = ['/var/www/html/scratch/' fnamestub];
    gmtname = [fname '.gmt'];
    psname = [fname '.ps'];
    pngname = [fname '.png'];
else
    % Perhaps running on GB's laptop, refer to consts file
    fname = [consts.scratchdir fnamestub];
    gmtname = [fname '.gmt'];
    psname = [fname '.eps'];
    pngname = [fname '.png'];
end

% open file for write
fid = fopen(gmtname,'w');

% Create GMT script
% Define format strings

lxs = ['"P/Pmean"'];
lys = ['"Elevation (m)"'];

rstring = ['-R' mmin_str '/' mmax_str '/0/' emax_str];
bstring = ['-Ba0.2g0.05:' lxs ':/a' ediv_str 'g' ediv_str ':' lys ':WeSn'];
jstring = '-JX2i/4i';

% Things that change in loops
sfs = {'St','Lm','LSDn'};
colors = {'0','0/200/0','0/0/255'};

% Tell shell path to GMT if necessary, this is for GB's Mac
if consts.isLocal == 1
    fprintf(fid,'%s\n\n','PATH=$PATH:/sw/bin');
    fprintf(fid,'%s\n\n','export PATH');
end

% Do everything three times
for a = 1:3
    cstr = colors{a};
    % Make the GMT basemap
    fprintf(fid,'%s\n','# ----------- Draw the plot axes ------------');
    if a == 1
        temp_string = [pstr 'psbasemap ' bstring ' ' jstring ' ' rstring ' -K -P ' formatstr ' > ' psname];
    else
        temp_string = [pstr 'psbasemap ' strrep(bstring,':WeSn',':weSn') ' ' jstring ' ' rstring ' -X2.1i -K -O -P ' formatstr ' >> ' psname];
    end
    fprintf(fid,'%s\n',temp_string);
    
    % Plot dark line at 1
    fprintf(fid,'%s\n','# ----------- Plot dark line at 1 ------------');
    temp_string = [pstr 'psxy ' rstring ' ' jstring mstr ' -P -O -K -W1p' wsep '0 << EOD >> ' psname];
    fprintf(fid,'%s\n',temp_string);
    fprintf(fid,'%s\n','1 0');
    fprintf(fid,'%s\n',['1 ' emax_str]);
    fprintf(fid,'%s\n','EOD');
   

    % Plot a label on the GMT basemap
    fprintf(fid,'%s\n','# ----------- Label by scaling method ------------');
    temp_string = [pstr 'pstext ' jstring ' -R0/1/0/1 -K -P -O << EOD >> ' psname];
    fprintf(fid,'%s\n',temp_string);
    temp_string = ['0.05 0.95 14 0 2 TL ' sfs{a}];
    fprintf(fid,'%s\n',temp_string);
    fprintf(fid,'%s\n','EOD');       
    
    % Plot error bars for exacts
    fprintf(fid,'%s\n','# ----------- Plot error bars for exact P ------------');
    temp_string = [pstr 'psxy ' rstring ' ' jstring mstr ' -P -O -K -W0.5p' wsep cstr ' << EOD >> ' psname];
    fprintf(fid,'%s\n',temp_string);
    eval(['this_x = r_' sfs{a} ';this_dx = dr_' sfs{a} ';']);
    for b = 1:length(elv_exact);
        if b > 1; fprintf(fid,'%s\n','>');end;
        fprintf(fid,'%0.5g %0.5g\n',[this_x(b)+this_dx(b) elv_exact(b)]); 
        fprintf(fid,'%0.5g %0.5g\n',[this_x(b)-this_dx(b) elv_exact(b)]); 
    end
    fprintf(fid,'%s\n','EOD');
    
    % Plot symbols for exact P
    fprintf(fid,'%s\n','# ----------- Plot symbols for exact P ------------');
    temp_string = [pstr 'psxy ' rstring ' ' jstring mstr ' -P -O -K -Sc0.1i -G' cstr ' << EOD >> ' psname];
    fprintf(fid,'%s\n',temp_string);
    fprintf(fid,'%0.5g %0.5g\n',[this_x elv_exact]');
    fprintf(fid,'%s\n','EOD');
    
    % Plot min Ps
    if ~isempty(elv_min)
        fprintf(fid,'%s\n','# ----------- Plot error bars for min P ------------');
        temp_string = [pstr 'psxy ' rstring ' ' jstring mstr ' -P -O -K -W0.5p' wsep cstr ' << EOD >> ' psname];
        fprintf(fid,'%s\n',temp_string);
        eval(['this_x = minr_' sfs{a} ';this_dx = mindr_' sfs{a} ';']);
        for b = 1:length(elv_min)
            if b > 1; fprintf(fid,'%s\n','>');end
            fprintf(fid,'%0.5g %0.5g\n',[this_x(b)+this_dx(b) elv_min(b)]); 
            fprintf(fid,'%0.5g %0.5g\n',[this_x(b)-this_dx(b) elv_min(b)]); 
        end
        fprintf(fid,'%s\n','EOD');
        
        % Plot symbols for min P
        fprintf(fid,'%s\n','# ----------- Plot symbols for min P ------------');
        temp_string = [pstr 'psxy ' rstring ' ' jstring mstr ' -P -O -K -St0.1i -G' cstr ' << EOD >> ' psname];
        fprintf(fid,'%s\n',temp_string);
        fprintf(fid,'%0.5g %0.5g\n',[this_x elv_min]');
        fprintf(fid,'%s\n','EOD');
    end
    
    % Plot max Ps
    if ~isempty(elv_max)
        fprintf(fid,'%s\n','# ----------- Plot error bars for max P ------------');
        temp_string = [pstr 'psxy ' rstring ' ' jstring mstr ' -P -O -K -W0.5p' wsep cstr ' << EOD >> ' psname];
        fprintf(fid,'%s\n',temp_string);
        eval(['this_x = maxr_' sfs{a} ';this_dx = maxdr_' sfs{a} ';']);
        for b = 1:length(elv_max)
            if b > 1; fprintf(fid,'%s\n','>');end;
            fprintf(fid,'%0.5g %0.5g\n',[this_x(b)+this_dx(b) elv_max(b)]); 
            fprintf(fid,'%0.5g %0.5g\n',[this_x(b)-this_dx(b) elv_max(b)]); 
        end
        fprintf(fid,'%s\n','EOD');
        
        % Plot symbols for max P
        fprintf(fid,'%s\n','# ----------- Plot symbols for max P ------------');
        temp_string = [pstr 'psxy ' rstring ' ' jstring mstr ' -P -O -K -Si0.1i -G' cstr ' << EOD >> ' psname];
        fprintf(fid,'%s\n',temp_string);
        fprintf(fid,'%0.5g %0.5g\n',[this_x elv_max]');
        fprintf(fid,'%s\n','EOD');
    end
end 
    

% Close out the GMT
fprintf(fid,'%s\n',[pstr 'psxy ' rstring ' ' jstring ' -Sc0.01i -G255 -P -O << EOD >> ' psname]);
fprintf(fid,'%s\n','1 1');
fprintf(fid,'%s\n','EOD');


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
    % On GB's Mac
    disp(['sh ' gmtname]);
    disp(['open ' psname]);
    % Running the convert command doesn't seem to work. 
end

% Return file name stub
out = fnamestub;

