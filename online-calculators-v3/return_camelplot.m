function out = return_camelplot(t,dt,s,xlimits)

% This code makes a camelplot. Called from within summarize_ages. 
%
% first args t and dt; this would typically be dtint
% third arg is summary data structure produced internally in
% summarize_ages, which is optional. 
%
% fourth arg xlimits allows control of x axis limits. 
%
% Returns a filename with no extension identifying the plot files. 
%
% Greg Balco - Berkeley Geochronology Center - December, 2016.
% Development, not licensed for use or distribution. 
%

% defaults
if nargin < 4
    xlimits = []; % forces compute x axis limits
end

isS = 1;
if nargin < 3
    isS = 0;
end
    
% Load consts for environment info
load consts_v3

% ------------ Determine axis limits -------------------

possx = [1e2 2e2 5e2 1e3 2e3 5e3 1e4 2e4 3e4 5e4 1e5 2e5 5e5 1e6 2e6 5e6 10e6 20e6 50e6];

minAge = min(t - 3.*dt);
maxAge = max(t + 3.*dt);

if isempty(xlimits)
    xlimits(1) = 0; 
    % Find first possible x limit option
    minopts = possx(min(find(possx > maxAge)));% Find first possible x limit option > maxAge
    if isempty(minopts)
        xlimits(2) = 50e6;
    else
        xlimits(2) = minopts; 
    end
end
    
% Also obtain x spacing
possspace = [1e2 1e3 2e3 5e3 1e4 2e4 5e4 1e5 2e5 5e5 1e6 2e6 5e6 10e6];
xsp = possspace(min(find(possspace > maxAge./5)));
    
%% ---------- Do camelplot calculation here

% this is guts of camelplot.m

numsamples = length(t);

x = linspace(minAge,maxAge,500);

for a = 1:numsamples
	xn = (x - t(a)) ./ dt(a);
	yn = exp(-0.5 * xn .^2) ./ (sqrt(2*pi) .* dt(a));  
    % Fancy-pants normalization so that the height of the total curve can't
    % exceed 1. Also corrects decreasing-importance-with-age effect.
    expected = 1 ./ (sqrt(2.*pi) .* t.*0.03); % expected height if 3% uncert
    pdfs(a,:) = yn./(expected(a).*numsamples);
end

% Summate

if numsamples > 1
    allpdf = sum(pdfs);
else
    allpdf = pdfs;
end

if isS == 1   
    if isfield(s,'use')
        if length(pdfs(s.use(s.strip.ok),:)) == 1
            strip_pdf = pdfs(s.use(s.strip.ok),:);
        else
            strip_pdf = sum(pdfs(s.use(s.strip.ok),:));
        end
    else
        if length(s.strip.ok) == 1
            strip_pdf = pdfs(s.strip.ok,:);
        else
            strip_pdf = sum(pdfs(s.strip.ok,:));
        end
    end
end

% ---------- end do camelplot calculation -----------------------

%% start GMT

% Using GMT on web server

pstr = 'gmt ';
wsep = ',';
formatstr = '--MAP_ANNOT_OFFSET_PRIMARY=0.1i --MAP_TICK_LENGTH_PRIMARY=-0.05i --MAP_FRAME_PEN=0.5p --FONT_ANNOT_PRIMARY=10p,2,black';

% Generate random number filename

%rand('seed',sum(100*clock));
uniqueid = int2str(round(rand(1)*1e6));

fnamestub = ['camelv3_' uniqueid];

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

%lxs = ['"Exposure age (yr)"'];
lxs = ['""'];
lys = ['""'];

ymax = max(allpdf).*1.1;
if isS == 1
    ymax = max([allpdf strip_pdf]).*1.1;
end

rstring = ['-R' sprintf('%0.0f',xlimits(1)) '/' sprintf('%0.0f',xlimits(2)) '/0/' sprintf('%0.2e',ymax)];
%bstring = ['-Ba' sprintf('%0.0f',xsp) 'g' sprintf('%0.0f',xsp) ':' lxs ':/:' lys ':weSn']; % what is this?
bstring = ['-Ba' sprintf('%0.0f',xsp) ':' lxs ':/:' lys ':weSn'];


% Tell shell path to GMT if necessary, this is for GB's Mac
if consts.isLocal == 1
    fprintf(fid,'%s\n\n','PATH=$PATH:/sw/bin');
    fprintf(fid,'%s\n\n','export PATH');
end


% Make the GMT basemap
fprintf(fid,'%s\n','# ----------- Draw the plot axes ------------');
jstring = '-JX6i/0.05i';
temp_string = [pstr 'psbasemap ' bstring ' ' jstring ' ' rstring ' -K -P ' formatstr ' > ' psname];
fprintf(fid,'%s\n',temp_string);

jstring = '-JX6i/1.5i';

if isS == 1
    % Plot a label on the GMT basemap
    temp_string = [pstr 'pstext ' jstring ' -R0/1/0/1 -K -P -O << EOD >> ' psname];
    fprintf(fid,'%s\n',temp_string);
    temp_string = ['0 1 12 0 2 TL ' s.label];
    fprintf(fid,'%s\n',temp_string);
    fprintf(fid,'%s\n','EOD');
end
       

% For each nuclide, plot the individual Gaussians and the summary one

% Actually, don't plot the individual Gaussians
% % Individual Gaussians
% for a = 1:numsamples;
%     eval(['thisy = pdf' int2str(a) ';']);
%     % File write 
%     fprintf(fid,'%s\n','# ----------- Plot individual Gaussian ------------');
%     if a == 1; 
%         temp_string = ['psxy ' rstring ' ' jstring ' -P -O -K -W0.5p/255/50/50 -Y0.1i << EOF >> ' psname];
%     else
%         temp_string = ['psxy ' rstring ' ' jstring ' -P -O -K -W0.5p/255/50/50 << EOF >> ' psname];
%     end;
%     fprintf(fid,'%s\n',temp_string);
%     fprintf(fid,'%0.5g %0.5g\n',[x',thisy']');
%     fprintf(fid,'%s\n','EOF');
% end;

cstr = '0';

% Determine what color to plot in...
if isS == 1
    if strfind(s.label,'St')
        cstr = '0';
    elseif strfind(s.label,'Lm')
        cstr = '0/200/0';
    elseif strfind(s.label,'LSDn')
        cstr = '0/0/255';
    end
end

% Plot total density
fprintf(fid,'%s\n','# ----------- Plot KDE for all samples ------------');
temp_string = [pstr 'psxy ' rstring ' ' jstring ' -P -O -K -W0.5p' wsep cstr ' -Y0.1i << EOD >> ' psname];
fprintf(fid,'%s\n',temp_string);
fprintf(fid,'%0.5g %0.5g\n',[x' allpdf']');
fprintf(fid,'%s\n','EOD');

if isS == 1
    % Plot total density
    fprintf(fid,'%s\n','# ----------- Plot KDE for non-outliers ------------');
    temp_string = [pstr 'psxy ' rstring ' ' jstring ' -P -O -K -W1p' wsep cstr ' << EOD >> ' psname];
    fprintf(fid,'%s\n',temp_string);
    fprintf(fid,'%0.5g %0.5g\n',[x' strip_pdf']');
    fprintf(fid,'%s\n','EOD');
    
    % Plot summary value and uncert
    fprintf(fid,'%s\n','# ----------- Plot summary value ------------');
    temp_string = [pstr 'psxy ' rstring ' ' jstring ' -P -O -K -W1p' wsep '100/100/100 << EOD >> ' psname];
    fprintf(fid,'%s\n',temp_string);
    fprintf(fid,'%0.5g %0.5g\n',[[s.sumval s.sumval]' [0 ymax]']');
    fprintf(fid,'%s\n','EOD');
    
    fprintf(fid,'%s\n','# ----------- Plot summary value ------------');
    temp_string = [pstr 'psxy ' rstring ' ' jstring ' -P -O -K -W0.5p' wsep '150/150/150 << EOD >> ' psname];
    fprintf(fid,'%s\n',temp_string);
    fprintf(fid,'%0.5g %0.5g\n',[[s.sumval-s.sumdel_int s.sumval-s.sumdel_int]' [0 ymax]']');
    fprintf(fid,'%s\n','>');
    fprintf(fid,'%0.5g %0.5g\n',[[s.sumval+s.sumdel_int s.sumval+s.sumdel_int]' [0 ymax]']');
    fprintf(fid,'%s\n','>');
    fprintf(fid,'%0.5g %0.5g\n',[[s.sumval+s.sumdel_ext s.sumval+s.sumdel_ext]' [0 ymax]']');
    fprintf(fid,'%s\n','>');
    fprintf(fid,'%0.5g %0.5g\n',[[s.sumval-s.sumdel_ext s.sumval-s.sumdel_ext]' [0 ymax]']');
    fprintf(fid,'%s\n','EOD');
    
end
    


% Close out the GMT
fprintf(fid,'%s\n',[pstr 'psxy ' rstring ' ' jstring ' -Sc0.01i -G255 -P -O << EOD >> ' psname]);
fprintf(fid,'%s\n','1 1');
fprintf(fid,'%s\n','EOD');


% close the file
fclose(fid);

if consts.isLocal == 0;
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
end;

% Return file name stub
out = fnamestub;

