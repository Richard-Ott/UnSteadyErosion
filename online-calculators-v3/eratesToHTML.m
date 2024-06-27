function out = eratesToHTML(erates,versions,splot_stubnames)

% This ingests the age results structure from get_erates_v3 and spits out an
% HTML string that is the results web page.  Various arguments. 
%
% Greg Balco - Berkeley Geochronology Center - March, 2017.
% Development, not licensed for use or distribution. 

load consts_v3; % Needed for lists of nuclides and scaling schemes...

dstring = [];

nl = char(10);

twidth = int2str(1000);
snwidth = int2str(100);
nnwidth = int2str(90);

% Write header

out = ['<!DOCTYPE html>' nl '<html>' nl '<head>' nl '<title>Erosion rate calculator v3 results</title>' nl ...       
    '<style>' nl '<!--table {}' nl '.title2{font-family:Arial;font-size:12.0pt;}' nl ...
    '.title{font-family:Arial;font-size:10.0pt;}' nl '.standard{font-family:Arial;font-size:8.0pt;}' nl ... 
    '-->' nl  '</style>' nl '</head>'];

out = [out nl nl '<body><center>'];

% Title table
out = [out nl '<table width=' twidth '>' nl '<tr>' nl '<td class="title2" colspan="2">' nl ...
    '<hr><i><b><p>Online erosion rate calculator v3 results</p></b></i><hr>' nl '</td>' nl '</tr>' ...	
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

% Display diagnostics. This comes from the erates.flags field.

out = [out nl nl '<table width=' twidth '>' nl '<tr>' nl '<td class=title width=200>' nl 'Diagnostics:' nl '</td>'];
dtext = [];
out = [out nl '<td class=standard>'];

dtext = [dtext erates.flags nl];

if isempty(dtext); dtext = 'No diagnostics'; end;
out = [out dtext nl '</td>' nl '</tr>'];
out = [out nl '<tr><td colspan=2><hr></td></tr>' nl '</table>'];

% Erosion rate tables. Loop through each scaling method. Within each scaling method,
% display all nuclide results for each sample sequentially.  

num_sschemes = 3; % Use only first three

% Build table header line

ths1 = ['<tr class=standard><td width=' snwidth '><b>Sample name</b><br><hr></td><td width=' nnwidth '><b>Nuclide</b><br><hr></td>'];
ths2 = '<tr class=standard><td></td><td></td>';
for a = 1:num_sschemes; 
    ths1 = [ths1 '<td colspan=4 align=center><b>' consts.sschemes{a} '</b><br><hr></td>'];
    ths2 = [ths2 '<td align=center><b>Erosion<br>rate<br>(g/cm2/yr)</b><br><hr></td>'];
    ths2 = [ths2 '<td align=center><b>Erosion<br>rate<br>(m/Myr)</b><br><hr></td>'];
    ths2 = [ths2 '<td align=center><b>Internal<br>uncert<br>(m/Myr)</b><br><hr></td>'];
    ths2 = [ths2 '<td align=center><b>External<br>uncert<br>(m/Myr)</b><br><hr></td>'];
end;
ths1 = [ths1 '</tr>'];
ths2 = [ths2 '</tr>'];

ncols = 2 + 4.*num_sschemes;

% Start building table...
out = [out nl nl '<table width=' twidth '>' nl '<tr>' nl '<td class=title colspan=' int2str(ncols) '>' nl '<b>Erosion rate results:</b>' nl '</td>' nl '</tr>' nl];
out = [out nl '<tr><td colspan=' int2str(ncols) '></td></tr>'];

out = [out nl ths1 nl ths2 nl]; % end with newline


% Now loop through ages. For each sample in ages, loop through different
% nuclide measurements and write a table row for all scaling schemes for each. 

for a = 1:erates.numsamples;
    nindices = find(erates.n.index == a);
    for b = 1:length(nindices);
        thisi = nindices(b);
        % Begin table row
        out = [out '<tr class=standard>' nl '<td>' erates.s.sample_name{a} '</td>'];
        % Use proper name for nuclide-mineral pair
        ni = find(strcmp(erates.n.nuclide{thisi},consts.nuclides));
        out = [out nl '<td>' consts.properName{ni} '</td>'];
        % Loop through all scaling schemes...
        for c = 1:num_sschemes;
            tfield = ['E_' consts.sschemes{c}];
            ifield = ['delE_int_' consts.sschemes{c}];
            efield = ['delE_ext_' consts.sschemes{c}];
            if isfield(erates.n,tfield);
                % Exists age for this scaling scheme
                out = [out nl '<td align=right>' sprintf('%0.3g',eval(['erates.n.' tfield '(thisi)'])) '</td>'];
                out = [out nl '<td align=right>' sprintf('%0.3g',eval(['erates.n.' tfield '(thisi).*1e4./erates.s.rho(erates.n.index(thisi))'])) '</td>'];
                if isfield(erates.n,ifield);
                    out = [out nl '<td align=right>' sprintf('%0.3g',eval(['erates.n.' ifield '(thisi).*1e4./erates.s.rho(erates.n.index(thisi))'])) '</td>'];
                else
                    out = [out nl '<td align=center> -- </td>'];
                end;
                if isfield(erates.n,efield);
                    out = [out nl '<td align=right>' sprintf('%0.3g',eval(['erates.n.' efield '(thisi).*1e4./erates.s.rho(erates.n.index(thisi))'])) '</td>'];
                else
                    out = [out nl '<td align=center> -- </td>'];
                end;
            else
                out = [out nl '<td align=center> -- </td><td align=center> -- </td><td align=center> -- </td><td align=center> -- </td>'];
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

if ~isempty(splot_stubnames)
    out = [out nl '<tr>' nl '<td class=standard>' nl '<center>In contrast to exposure-age results, here multiple-nuclide diagrams use only St scaling for normalization.</center>' nl '</td>' nl '</tr>'];
end;

for a = 1:length(splot_stubnames);
    out = [out nl '<tr>' nl '<td class=standard>' nl '<center><img src="' consts.plotURL splot_stubnames{a} '.png">'];
    out = [out nl '</td>' nl '<td class=standard valign=middle>'];
    out = [out '<a href="' consts.plotURL splot_stubnames{a} '.ps">Postscript file</a><br>'];
    out = [out '<a href="' consts.plotURL splot_stubnames{a} '.gmt">GMT script</a></center>'];
    out = [out nl '</td>' nl '</tr>'];
    out = [out nl '<tr><td colspan=2><hr></td></tr>'];
end;

% Close table

out = [out nl '</table>'];

% End HTML
out = [out nl '</center></body>' nl '</html>'];



