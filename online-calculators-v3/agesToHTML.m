function out = agesToHTML(ages,versions,splot_stubnames,tdplot_stubname,compplot_stubname,summary,calib_trace)

% This ingests the age results structure from get_ages_v3 and spits out an
% HTML string that is the results web page.  Various arguments. 
%
% Greg Balco - Berkeley Geochronology Center - February, 2015.
% Development, not licensed for use or distribution. 

load consts_v3; % Needed for lists of nuclides and scaling schemes...

dstring = [];

nl = char(10);

twidth = int2str(1000);
snwidth = int2str(100);

% Write header

out = ['<!DOCTYPE html>' nl '<html>' nl '<head>' nl '<title>Exposure age calculator v3 results</title>' nl ...       
    '<style>' nl '<!--table {}' nl '.title2{font-family:Arial;font-size:12.0pt;}' nl ...
    '.title{font-family:Arial;font-size:10.0pt;}' nl '.standard{font-family:Arial;font-size:8.0pt;}' nl ... 
    '-->' nl  '</style>' nl '</head>'];

out = [out nl nl '<body><center>'];

% Title table
out = [out nl '<table width=' twidth '>' nl '<tr>' nl '<td class="title2" colspan="2">' nl ...
    '<hr><i><b><p>Online exposure age calculator v3 results</p></b></i><hr>' nl '</td>' nl '</tr>' ...	
    '</table>'];

% Version table

out = [out nl nl '<table width=' twidth '>' nl '<tr>' nl '<td class=title width=200 valign=top>' ];
out = [out nl 'Version info:' nl '</td>' nl '<td class=standard>'];
fns = fieldnames(versions);
for a = 1:length(fns);
    thisstr = [nl fns{a} ': ' eval(['versions.' fns{a}]) nl '<br>' nl];
    out = [out thisstr];
end;
out = [out nl '</td>' nl '</tr>'];
out = [out nl '<tr><td colspan=2><hr></td></tr>' nl '</table>'];

% Display diagnostics. This comes from the ages{a}.flags field.

out = [out nl nl '<table width=' twidth '>' nl '<tr>' nl '<td class=title width=200>' nl 'Diagnostics:' nl '</td>'];
dtext = [];
out = [out nl '<td class=standard>'];

dtext = [dtext ages.flags nl];

if isempty(dtext); dtext = 'No diagnostics'; end;
out = [out dtext nl '</td>' nl '</tr>'];
out = [out nl '<tr><td colspan=2><hr></td></tr>' nl '</table>'];

% Calibration data

out = [out nl nl '<table width=' twidth '>' nl '<tr>' nl '<td class=title width=200 valign=top>' nl 'Calibration data:' nl '</td>'];
out = [out nl '<td class=standard>'];
out = [out nl '<p><b>Calibration data set: </b>' calib_trace.calibration_name '</p>'];
out = [out nl '<p><b>Trace string: </b>' calib_trace.trace_string '</p>'];
if ~isempty(calib_trace.nuclide);
    out = [out nl '<p>Nuclide/target actually calibrated: <b>' calib_trace.nuclide '</b></p>'];
    out = [out nl '<p>Age calculations with all other nuclides are unaffected; they retain default parameters.</p>'];
end;    
out = [out nl '</td>' nl '</tr>'];
out = [out nl '<tr><td colspan=2><hr></td></tr>' nl '</table>'];

% Age tables. Loop through each scaling scheme. Within each scaling scheme,
% display all nuclide results for each sample sequentially.  

% Build table header line

ths1 = ['<tr class=standard><td width=' snwidth '><b>Sample name</b><br><hr></td><td><b>Nuclide</b><br><hr></td>'];
ths2 = '<tr class=standard><td></td><td></td>';
for a = 1:length(consts.sschemes); 
    ths1 = [ths1 '<td colspan=3 align=center><b>' consts.sschemes{a} '</b><br><hr></td>'];
    ths2 = [ths2 '<td align=right><b>Age (yr)</b><br><hr></td>'];
    ths2 = [ths2 '<td align=right><b>Interr (yr)</b><br><hr></td>'];
    ths2 = [ths2 '<td align=right><b>Exterr (yr)</b><br><hr></td>'];
end;
ths1 = [ths1 '</tr>'];
ths2 = [ths2 '</tr>'];

ncols = 2 + 3.*length(consts.sschemes);

% Start building table...
out = [out nl nl '<table width=' twidth '>' nl '<tr>' nl '<td class=title colspan=' int2str(ncols) '>' nl '<b>Exposure age results:</b>' nl '</td>' nl '</tr>' nl];
out = [out nl '<tr><td colspan=' int2str(ncols) '></td></tr>'];

out = [out nl ths1 nl ths2 nl]; % end with newline


% Now loop through ages. For each sample in ages, loop through different
% nuclide measurements and write a table row for all scaling schemes for each. 

for a = 1:ages.numsamples;
    nindices = find(ages.n.index == a);
    for b = 1:length(nindices);
        thisi = nindices(b);
        % Begin table row
        out = [out '<tr class=standard>' nl '<td>' ages.s.sample_name{a} '</td>'];
        % Use proper name for nuclide-mineral pair
        ni = find(strcmp(ages.n.nuclide{thisi},consts.nuclides));
        out = [out nl '<td>' consts.properName{ni} '</td>'];
        % Loop through all scaling schemes...
        for c = 1:length(consts.sschemes);
            tfield = ['t_' consts.sschemes{c}];
            ifield = ['delt_int_' consts.sschemes{c}];
            efield = ['delt_ext_' consts.sschemes{c}];
            if isfield(ages.n,tfield);
                % Exists age for this scaling scheme
                out = [out nl '<td align=right>' sprintf('%0.0f',eval(['ages.n.' tfield '(thisi)'])) '</td>'];
                if isfield(ages.n,ifield);
                    out = [out nl '<td align=right>' sprintf('%0.0f',eval(['ages.n.' ifield '(thisi)'])) '</td>'];
                else
                    out = [out nl '<td align=center> -- </td>'];
                end;
                if isfield(ages.n,efield);
                    out = [out nl '<td align=right>' sprintf('%0.0f',eval(['ages.n.' efield '(thisi)'])) '</td>'];
                else
                    out = [out nl '<td align=center> -- </td>'];
                end;
            else
                out = [out nl '<td align=center> -- </td><td align=center> -- </td><td align=center> -- </td>'];
            end;
        end;
    end;
    % After each sample, write blank table row
    out = [out nl '<tr height=3><td colspan=' int2str(ncols) '></td></tr>'];
end;

% Close table
out = [out nl '<tr><td colspan=' int2str(ncols) '><hr></td></tr>'];
out = [out nl '</table>'];

% Now need to put in splots and other plot data
% This is another table

% Start table

out = [out nl nl '<table width=' twidth '>'];





for a = 1:length(splot_stubnames);
    out = [out nl '<tr>' nl '<td class=standard>' nl '<center><img src="' consts.plotURL splot_stubnames{a} '.png">'];
    out = [out nl '</td>' nl '<td class=standard valign=middle>'];
    out = [out '<a href="' consts.plotURL splot_stubnames{a} '.ps">Postscript file</a><br>'];
    out = [out '<a href="' consts.plotURL splot_stubnames{a} '.gmt">GMT script</a></center>'];
    out = [out nl '</td>' nl '</tr>'];
    out = [out nl '<tr><td colspan=2><hr></td></tr>'];
end;

if ~isempty(tdplot_stubname);
    out = [out nl '<tr>' nl '<td class=standard>' nl '<center><img src="' consts.plotURL 'tdplotv3_' tdplot_stubname '.png">'];
    out = [out nl '</td>' nl '<td class=standard valign=middle>'];
    out = [out '<a href="' consts.plotURL tdplot_stubname '.ps">Postscript file</a><br>'];
    out = [out '<a href="' consts.plotURL tdplot_stubname '.gmt">GMT script</a></center>'];
    out = [out nl '</td>' nl '</tr>'];
    out = [out nl '<tr><td colspan=2><hr></td></tr>'];
end;

if ~isempty(compplot_stubname);   
    out = [out nl '<tr>' nl '<td class=standard>' nl '<center><img src="' consts.plotURL compplot_stubname '.png">'];
    out = [out nl '</td>' nl '<td class=standard valign=middle>'];
    out = [out '<a href="' consts.plotURL compplot_stubname '.ps">Postscript file</a><br>'];
    out = [out '<a href="' consts.plotURL compplot_stubname '.gmt">GMT script</a></center>'];
    out = [out nl '</td>' nl '</tr>'];
    out = [out nl '<tr><td colspan=2><hr></td></tr>'];
end;

% Close table

out = [out nl '</table>'];

% If summary data included, make another table with summary data

if ~isempty(summary);
    
    % Start table
    out = [out nl nl '<table width=' twidth '>'];
    out = [out nl '<tr><td class=title>Summary statistics:</td></tr>'];
    
    % Loop; if only one nuclide, 
    numnuclides = length(unique(ages.n.nindex));
    nloop = {'all'};
    if numnuclides > 1;
         nloop = [{'all'} consts.nuclides];
    end;

    for sf = {'St','Lm','LSDn'};
        for n = nloop;
            eval(['thiss = summary.' char(n) '.' char(sf) ';']);           
            if ~isempty(thiss);
                % Make a table row
                if strcmp(char(n),'all');
                    idstr = ['<b>' char(sf) ' scaling, all nuclides</b></br>'];
                else
                    idstr = ['<b> ' char(sf) ' scaling, ' char(n) '</b></br>'];
                end;
                out = [out nl '<tr><td class=standard><p>' idstr '</p><blockquote>' thiss.text '</blockquote></td></tr>'];
            end;
        end;
    end;
    out = [out nl '<tr><td><hr></td><tr></table>'];
    
    % Add camelplots
    if isfield(summary,'camels');
        out = [out nl nl '<table width=' twidth '>']; 
        for a = 1:length(summary.camels);
            out = [out nl '<tr>' nl '<td class=standard>' nl '<center><img src="' consts.plotURL summary.camels{a} '.png">'];
            out = [out nl '<br>'];
            out = [out '<p><a href="' consts.plotURL summary.camels{a} '.ps">Postscript file</a> -- '];
            out = [out '<a href="' consts.plotURL summary.camels{a} '.gmt">GMT script</a></p></center>'];
            out = [out nl '</center></td>' nl '</tr>'];
            out = [out nl '<tr><td><hr></td></tr>'];
        end;
    end;
    
    % Close table
    %out = [out nl '<tr><td><hr></td><tr>']
    out = [out nl '</table>'];
    
end;

% End HTML
out = [out nl '</center></body>' nl '</html>'];



