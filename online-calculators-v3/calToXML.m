function out = calToXML(ages,versions,summary,Pplot_stubname)

% This ingests the age results structure from get_ages_v3 and spits out an
% XML result string. 

if nargin < 4; Pplot_stubname = {}; end
if nargin < 3; summary = [];end

load consts_v3; % Needed to define scaling scheme loop

dstring = [];

out =  '<calcs_v3_cal_data>';

rr = '';

for a = 1:ages.numsamples
    % Which measurements belong to this sample
    nindices = find(ages.n.index == a);
    % Open results field
    out = strcat(out,'<calibrationResult>',rr);
    
    % Report sample name 
    out = strcat(out,'<sample_name>',ages.s.sample_name{a},'</sample_name>',rr);
    
    popts = {'','min','max'};
    
    % Loop through scaling method and exact/min/max possibilities
    for b = 1:length(nindices)
        for p = 1:length(popts)
            this_nuclide = char(ages.n.nuclide(nindices(b)));
            for s = 1:length(consts.sschemes)
                % Create XML field names, which are the same as field names in
                % 'ages' structure
                Pname = strcat('calc_', popts{p},'P_',consts.sschemes{s});
                dPname = strcat('calc_del',popts{p},'P_',consts.sschemes{s});          
                if isfield(ages.n,Pname)
                    % Get data to put in them, if it exists
                    eval(['this_P = ages.n.' Pname '(nindices(b));']);
                    eval(['this_dP = ages.n.' dPname '(nindices(b));']);
                    if this_P > 0 && this_P < Inf
                        % Write XML fields
                        out = strcat(out,'<',Pname,'>',sprintf('%0.3f',this_P),'</',Pname,'>');
                        out = strcat(out,'<',dPname,'>',sprintf('%0.3f',this_dP),'</',dPname,'>');
                    end
                end
            end
        end
    end
            
    % Close result field
    out = [out '</calibrationResult>' rr];
    

end

% Append elements for plot names if present

if ~isempty(Pplot_stubname)
    out = [out '<ploturlstub>' consts.plotURL Pplot_stubname '</ploturlstub>'];
end

% Append summary data if present

if ~isempty(summary)

    if isfield(summary,'trace_string')
        out = [out '<trace_string>' summary.trace_string '</trace_string>'];
    end
    if isfield(summary,'calibration_name')
        out = [out '<calibration_name>' summary.calibration_name '</calibration_name>'];
    end
    if isfield(summary,'nuclide')
        out = [out '<nuclide>' summary.nuclide '</nuclide>'];
    end
    
    for a = 1:length(consts.sschemes)
        if isfield(summary,consts.sschemes{a})
            eval(['this_sumval = summary.' consts.sschemes{a} '.avgP;']);
            eval(['this_sumdel = summary.' consts.sschemes{a} '.use_uncert;']);
            eval(['this_sumtext = summary.' consts.sschemes{a} '.text;']);
            
            out = [out '<summary_value_' consts.sschemes{a} '>' sprintf('%0.3f',this_sumval) '</summary_value_' consts.sschemes{a} '>'];
            out = [out '<summary_uncert_' consts.sschemes{a} '>' sprintf('%0.3f',this_sumdel) '</summary_uncert_' consts.sschemes{a} '>'];
            out = [out '<summary_text_' consts.sschemes{a} '>' this_sumtext '</summary_text_' consts.sschemes{a} '>'];
        end
    end
end
    
% Also dump diagnostics

if isempty(ages.flags)
    dstring = 'No diagnostics'; 
else
    dstring = strrep(ages.flags,'<br>','');
end

out = [out '<diagnostics>' rr dstring  rr '</diagnostics>' rr];

% Also dump version info
% Should probably have a versions.summarize


if ~isempty(versions)
    out = [out '<version>'];
    versions_fields = fieldnames(versions);

    for a = 1:length(versions_fields)
        eval(['this_value = versions.' versions_fields{a} ';']);
        out = [out '<' versions_fields{a} '>' this_value '</' versions_fields{a} '>' rr];
    end
    out = [out '</version>'];
end

% Close 
out = [out '</calcs_v3_cal_data>'];



