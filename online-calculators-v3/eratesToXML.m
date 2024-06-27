function out = eratesToXML(erates,versions,splot_stubnames)

% This ingests the age results structure from get_erates_v3 and spits out an
% XML result string. 

if nargin < 3; splot_stubnames = {};end;

load consts_v3; % Needed to define scaling scheme loop

dstring = [];

out =  '<calcs_v3_erosion_data>';

rr = '';

for a = 1:erates.numsamples;
    % Which measurements belong to this sample
    nindices = find(erates.n.index == a);
    % Open results field
    out = strcat(out,'<erosionRateResult>',rr);
    
    % Report sample name 
    out = strcat(out,'<sample_name>',erates.s.sample_name{a},'</sample_name>',rr);
    % Report density
    out = strcat(out,'<sample_density>',sprintf('%0.2f',erates.s.rho(a)),'</sample_density>',rr);
    
    % Loop through nuclide-scaling method possibilities
    for b = 1:length(nindices);
        this_nuclide = char(erates.n.nuclide(nindices(b)));
        % Wrap result for each nuclide
        out = strcat(out,'<nuclide_result><nuclide>',this_nuclide,'</nuclide>');
        % Loop through scaling methods
        for s = 1:length(consts.sschemes);
            % Create XML field names
            Ename = strcat('E_gcm2_',consts.sschemes{s});
            dE_int_name = strcat('delE_gcm2_',consts.sschemes{s});
            dE_ext_name = strcat('delE_gcm2_',consts.sschemes{s});           
            if isfield(erates.n,['E_' consts.sschemes{s}]);
                % Get data to put in them, if it exists
                eval(['this_E = erates.n.E_' consts.sschemes{s} '(nindices(b));']);
                eval(['this_dE_int = erates.n.delE_int_' consts.sschemes{s} '(nindices(b));']);
                eval(['this_dE_ext = erates.n.delE_ext_' consts.sschemes{s} '(nindices(b));']);
                % Write erosion rate fields
                out = strcat(out,'<',Ename,'>',sprintf('%0.4e',this_E),'</',Ename,'>');
                out = strcat(out,'<',dE_int_name,'>',sprintf('%0.4e',this_dE_int),'</',dE_int_name,'>');
                out = strcat(out,'<',dE_ext_name,'>',sprintf('%0.4e',this_dE_ext),'</',dE_ext_name,'>',rr); 
            end;
        end;
        % End wrap
        out = strcat(out,'</nuclide_result>');
    end;
            
    % Close result field
    out = [out '</erosionRateResult>' rr];
    

end;

% Append elements for plot names if present

if ~isempty(splot_stubnames);
    for a = 1:length(splot_stubnames);
        out = [out '<ploturlstub>' consts.plotURL splot_stubnames{a} '</ploturlstub>'];
    end;
end;
    
% Also dump diagnostics

if isempty(erates.flags); 
    dstring = 'No diagnostics'; 
else
    dstring = strrep(erates.flags,'<br>','');
end;

out = [out '<diagnostics>' rr dstring  rr '</diagnostics>' rr];

% Also dump version info

out = [out '<version><wrapper>' versions.wrapper '</wrapper><validate>' ...
    versions.validate '</validate><erates>' versions.erates ...
    '</erates><muons>' versions.muons '</muons><consts>' consts.version ...
    '</consts></version>' rr];

% Close 
out = [out '</calcs_v3_erosion_data>'];



