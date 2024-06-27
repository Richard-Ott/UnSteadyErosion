function out = return_sat14plot(ages,consts)

% This code makes a saturated C-14 plot. 
%
% first arg is ages structure, second is consts struct
%
% Returns a filename with no extension identifying the plot files. 
%
% Greg Balco - Berkeley Geochronology Center - March, 2018.
% Development, not licensed for use or distribution. 
%

% Load consts if necessary
if nargin < 2
    load consts_v3
end
    
% ----- extract data -----------------------------------

is14 = strcmp(ages.n.nuclide,'N14quartz');
px = ages.n.N(is14);
pdx = ages.n.delN(is14);
py = ages.s.elv(ages.n.index(is14));

% ------------ Determine axis limits -------------------

possy = [500 1000 2000 3000 4000];
possx = [3e5 5e5 1e6 2e6 3e6];

possyi = [100 200 500 500 1000];
possxi = [1e5 1e5 2e5 5e5 5e5];

uselim = min(find(possy > max(py)));
if isempty(uselim); uselim = 4; end

xmax = possx(uselim);
ymax = possy(uselim);
xsp = possxi(uselim);
ysp = possyi(uselim);

%% start GMT

% Using GMT on web server

pstr = 'gmt ';
wsep = ',';
formatstr = '--MAP_ANNOT_OFFSET_PRIMARY=0.1i --MAP_TICK_LENGTH_PRIMARY=-0.05i --MAP_FRAME_PEN=0.5p --FONT_ANNOT_PRIMARY=10p,2,black --FONT_LABEL=14p,3,black';
mstr = '';

% Generate random number filename

%rand('seed',sum(100*clock));
uniqueid = int2str(round(rand(1)*1e6));

fnamestub = ['sat14v3_' uniqueid];

if consts.isLocal == 0
    % Assume running on hess.ess
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

lxs = ['"[C-14] (atoms/g)"'];
lys = ['"Elevation (m)"'];


rstring = ['-R0/' sprintf('%0.0f',xmax) '/0/' sprintf('%0.0f',ymax)];
bstring = ['-Ba' sprintf('%0.0f',xsp) 'g' sprintf('%0.0f',xsp) ':' lxs ':/a' sprintf('%0.0f',ysp) 'g' sprintf('%0.0f',ysp) ':' lys ':WeSn']; % what is this?
jstring = '-JX5.5i/3.5i';

% Tell shell path to GMT if necessary, this is for GB's Mac
if consts.isLocal == 1
    fprintf(fid,'%s\n\n','PATH=$PATH:/sw/bin');
    fprintf(fid,'%s\n\n','export PATH');
end


% Make the GMT basemap
fprintf(fid,'%s\n','# ----------- Draw the plot axes ------------');

temp_string = [pstr 'psbasemap ' bstring ' ' jstring ' ' rstring ' -K -P ' formatstr ' > ' psname];
fprintf(fid,'%s\n',temp_string);

% Plot the saturation curve
thisy = consts.sat14.elv;

fprintf(fid,'%s\n','# ----------- Plot saturation curve ------------');
thisx = consts.sat14.mu.sat14 + consts.sat14.sp.sat14.*consts.refP_LSDn(5);
temp_string = [pstr 'psxy ' rstring ' ' jstring ' -P -O -K -W1p' wsep '0 << EOD >> ' psname];
fprintf(fid,'%s\n',temp_string);
fprintf(fid,'%0.2e %0.0f\n',[thisx' thisy']');
fprintf(fid,'%s\n','EOD');

fprintf(fid,'%s\n','# ----------- Plot saturation curve error envelope ------------');
satx = thisx;
thisx = satx.*(1 - consts.delrefP_LSDn(5)./consts.refP_LSDn(5));
temp_string = [pstr 'psxy ' rstring ' ' jstring ' -P -O -K -W0.5p' wsep '0 << EOD >> ' psname];
fprintf(fid,'%s\n',temp_string);
fprintf(fid,'%0.2e %0.0f\n',[thisx' thisy']');
fprintf(fid,'%s\n','EOD');

thisx = satx.*(1 + consts.delrefP_LSDn(5)./consts.refP_LSDn(5));
temp_string = [pstr 'psxy ' rstring ' ' jstring ' -P -O -K -W0.5p' wsep '0 << EOD >> ' psname];
fprintf(fid,'%s\n',temp_string);
fprintf(fid,'%0.2e %0.0f\n',[thisx' thisy']');
fprintf(fid,'%s\n','EOD');

iccol = '200/200/200';

fprintf(fid,'%s\n','# ----------- Plot 2 ka curve ------------');
x2k = consts.sat14.mu.t2k + consts.sat14.sp.t2k.*consts.refP_LSDn(5);
temp_string = [pstr 'psxy ' rstring ' ' jstring ' -P -O -K -W0.5p' wsep iccol ' << EOD >> ' psname];
fprintf(fid,'%s\n',temp_string);
fprintf(fid,'%0.2e %0.0f\n',[x2k' thisy']');
fprintf(fid,'%s\n','EOD');

fprintf(fid,'%s\n','# ----------- Plot 4 ka curve ------------');
x4k = consts.sat14.mu.t4k + consts.sat14.sp.t4k.*consts.refP_LSDn(5);
temp_string = [pstr 'psxy ' rstring ' ' jstring ' -P -O -K -W0.5p' wsep iccol ' << EOD >> ' psname];
fprintf(fid,'%s\n',temp_string);
fprintf(fid,'%0.2e %0.0f\n',[x4k' thisy']');
fprintf(fid,'%s\n','EOD');

fprintf(fid,'%s\n','# ----------- Plot 6 ka curve ------------');
x6k = consts.sat14.mu.t6k + consts.sat14.sp.t6k.*consts.refP_LSDn(5);
temp_string = [pstr 'psxy ' rstring ' ' jstring ' -P -O -K -W0.5p' wsep iccol ' << EOD >> ' psname];
fprintf(fid,'%s\n',temp_string);
fprintf(fid,'%0.2e %0.0f\n',[x6k' thisy']');
fprintf(fid,'%s\n','EOD');

fprintf(fid,'%s\n','# ----------- Plot 8 ka curve ------------');
x8k = consts.sat14.mu.t8k + consts.sat14.sp.t8k.*consts.refP_LSDn(5);
temp_string = [pstr 'psxy ' rstring ' ' jstring ' -P -O -K -W0.5p' wsep iccol ' << EOD >> ' psname];
fprintf(fid,'%s\n',temp_string);
fprintf(fid,'%0.2e %0.0f\n',[x8k' thisy']');
fprintf(fid,'%s\n','EOD');

fprintf(fid,'%s\n','# ----------- Plot 10 ka curve ------------');
x10k = consts.sat14.mu.t10k + consts.sat14.sp.t10k.*consts.refP_LSDn(5);
temp_string = [pstr 'psxy ' rstring ' ' jstring ' -P -O -K -W0.5p' wsep iccol ' << EOD >> ' psname];
fprintf(fid,'%s\n',temp_string);
fprintf(fid,'%0.2e %0.0f\n',[x10k' thisy']');
fprintf(fid,'%s\n','EOD');

fprintf(fid,'%s\n','# ----------- Plot 15 ka curve ------------');
x15k = consts.sat14.mu.t15k + consts.sat14.sp.t15k.*consts.refP_LSDn(5);
temp_string = [pstr 'psxy ' rstring ' ' jstring ' -P -O -K -W0.5p' wsep iccol ' << EOD >> ' psname];
fprintf(fid,'%s\n',temp_string);
fprintf(fid,'%0.2e %0.0f\n',[x15k' thisy']');
fprintf(fid,'%s\n','EOD');

fprintf(fid,'%s\n','# ----------- Plot 20 ka curve ------------');
x20k = consts.sat14.mu.t20k + consts.sat14.sp.t20k.*consts.refP_LSDn(5);
temp_string = [pstr 'psxy ' rstring ' ' jstring ' -P -O -K -W0.5p' wsep iccol ' << EOD >> ' psname];
fprintf(fid,'%s\n',temp_string);
fprintf(fid,'%0.2e %0.0f\n',[x20k' thisy']');
fprintf(fid,'%s\n','EOD');

fprintf(fid,'%s\n','# ----------- Plot error bars on data ------------');

temp_string = [pstr 'psxy ' rstring ' ' jstring mstr ' -P -O -K -W0.5p' wsep '0/0/255 << EOD >> ' psname];
fprintf(fid,'%s\n',temp_string);
for b = 1:length(px)
	if b > 1; fprintf(fid,'%s\n','>');end
        fprintf(fid,'%0.5g %0.5g\n',[px(b)+pdx(b) py(b)]); 
        fprintf(fid,'%0.5g %0.5g\n',[px(b)-pdx(b) py(b)]); 
end
fprintf(fid,'%s\n','EOD');

fprintf(fid,'%s\n','# ----------- Plot C-14 data ------------');
temp_string = [pstr 'psxy ' rstring ' ' jstring mstr ' -P -O -K -Sc0.1i -W0.5p' wsep '36-.85-.5 -G36-.85-1 << EOD >> ' psname];
fprintf(fid,'%s\n',temp_string);
fprintf(fid,'%0.5g %0.5g\n',[px py]');
fprintf(fid,'%s\n','EOD');


% Plot label text for isochrons
fprintf(fid,'%s\n','# ----------- Plot label text for isochrons ------------');
temp_string = ['gmt pstext -R0/' sprintf('%0.0f',xmax) '/0/1 -JX5.5i/0.3i -P -O -K -F+f10p,Helvetica+jBC -Y3.6i << EOD >> ' psname];
fprintf(fid,'%s\n',temp_string);

% Obtain text locations for labeling isochrons
ts = [2 4 6 8 10 15 20 100];
ils = {'2 ka','4 ka','6 ka','8','10','15','20','(sat)'};

for a = 1:length(ts)
    if a < length(ts)
        eval(['thisx = x' int2str(ts(a)) 'k;']);
    else
        thisx = satx;
    end
    textx = 1.02*interp1(thisy,thisx,ymax);
    temp_string = [sprintf('%0.0f',textx) ' 0.1 ' ils{a}];
    fprintf(fid,'%s\n',temp_string);
end
fprintf(fid,'%s\n','EOD');

% Close out the GMT
fprintf(fid,'%s\n',[pstr 'psxy -R0/1/0/1 -JX0.1i/0.3i -Sc0.01i -G255 -P -O << EOD >> ' psname]);
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

