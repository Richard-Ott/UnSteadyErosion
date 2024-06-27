function retstr = calToHTML(versions,summary,Pplot_stubname)

% This generates HTML for a v3 exposure age input page using results of
% production rate calibration. 

load consts_v3;

% This just adds a section to the v3_age_in.php file. 
% Open that file and get text. 

fid = fopen('v3_age_in.html','r');
base_text = char((fread(fid,inf,'char'))');
fclose(fid);

% Start HTML to add. This is a closed set of table rows. 

nl = char(10);

% 1. Deal with hidden vars to pass

s = ['<input type="hidden" name="trace_string" value="' summary.trace_string '">' nl];
s = [s '<input type="hidden" name="calib_name" value="' summary.calibration_name '">' nl];
s = [s '<input type="hidden" name="nuclide_name" value="' char(summary.nuclide) '">' nl];

sfs = consts.sschemes;
for a = 1:length(sfs);
    if isfield(summary,sfs{a});
        eval(['thisP = summary.' sfs{a} '.avgP;']);
        eval(['thisdP = summary.' sfs{a} '.use_uncert;']);
        s = [s '<input type="hidden" name="P_' sfs{a} '" value="' sprintf('%0.3f',thisP) '">' nl];
        s = [s '<input type="hidden" name="delP_' sfs{a} '" value="' sprintf('%0.3f',thisdP) '">' nl];
    end;
end;
    
% Done with hidden vars. Make results table. 

% Header table row
s = [s '<tr><td><hr></td></tr>' nl '<tr><td class=title2><i><b>Production rate calibration results</b></i></td></tr>' nl];

% Table row with version info
s = [s '<tr><td><i>' nl];
s = [s nl '<p><b>Version info:</i></b></p>' nl '<blockquote><i>' nl];
fns = fieldnames(versions);
for a = 1:length(fns);
    thisstr = [nl fns{a} ': ' eval(['versions.' fns{a}]) nl '<br>' nl];
    s = [s thisstr];
end;
s = [s nl '</i></blockquote>' nl '</td></tr>'];

% Table row with calibration data set name

s = [s '<tr><td>' nl];
s = [s '<p><b><i>Calibration data set name:</b> ' summary.calibration_name '</i></p>' nl];
s = [s '</td></tr>'];

% Table row with calibration results

s = [s '<tr><td>' nl];
s = [s nl '<p><b><i>Calibration results:</i></b></p>' nl '<blockquote><i>' nl];

nindex = find(strcmp(summary.nuclide,consts.nuclides));

if nindex == 2 || nindex == 3;
    ns = 'He-3 (px/ol)';
else
    ns = consts.properName{nindex};
end;

s = [s '<p>Nuclide/target actually calibrated: <b> ' ns '</b></p>' nl];
s = [s '<p>Age calculations with all other nuclides are unaffected; they retain default parameters.</p>' nl];

for a = 1:length(sfs);
    if isfield(summary,sfs{a});
        s = [s '<p>'];
        eval(['s = [s summary.' sfs{a} '.text nl];']);
        s = [s '</p>' nl];
    end;
end;

s = [s '</i></blockquote></td></tr>' nl];
s = [s nl '<tr><td><hr></td></tr>'];

% Table row with pPlot figure

if ~isempty(Pplot_stubname);
    s = [s nl '<tr>' nl '<td class=standard>' nl '<center><img src="' consts.plotURL Pplot_stubname '.png">'];
    s = [s nl '</td></tr>' nl '<tr><td class=standard valign=middle>'];
    s = [s '<center><p><i><a href="' consts.plotURL Pplot_stubname '.ps">Postscript file</a> '];
    s = [s ' <a href="' consts.plotURL Pplot_stubname '.gmt">GMT script</a></i></p></center>'];
    s = [s nl '</td>' nl '</tr>'];
    s = [s nl '<tr><td><i><p><center>'];
    % insert explanation of plot here
    s = [s nl '</center></p></i></td></tr>'];
end;


% Now insert that in the file. 

retstr = strrep(base_text,'<! -- Calibration info goes here -->',s);

% Also provide absolute link for form target so this page can be saved locally. 

retstr = strrep(retstr,'/cgi-bin/matweb',['http://' consts.machine_name '/cgi-bin/matweb']);



