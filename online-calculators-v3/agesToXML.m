function out = agesToXML(ages,versions,splot_stubnames,summary)

% This ingests the age results structure from get_ages_v3 and spits out an
% XML result string. 

if nargin < 4; summary = []; end
if nargin < 3; splot_stubnames = {};end
if nargin < 2; versions = []; end

load consts_v3; % Needed to define scaling scheme loop

dstring = [];

out =  '<calcs_v3_age_data>';

rr = '';

for a = 1:ages.numsamples
    % Which measurements belong to this sample
    nindices = find(ages.n.index == a);
    % Open results field
    out = strcat(out,'<exposureAgeResult>',rr);
    
    % Report sample name 
    out = strcat(out,'<sample_name>',ages.s.sample_name{a},'</sample_name>',rr);
    
    % Loop through nuclide-scaling scheme possibilities
    for b = 1:length(nindices)
        this_nuclide = char(ages.n.nuclide(nindices(b)));
        for s = 1:length(consts.sschemes)
            % Create XML field names
            tname = strcat('t',this_nuclide(2:end),'_',consts.sschemes{s});
            dt_int_name = strcat('delt',this_nuclide(2:end),'_int_',consts.sschemes{s});
            dt_ext_name = strcat('delt',this_nuclide(2:end),'_ext_',consts.sschemes{s});
            
            if isfield(ages.n,['t_' consts.sschemes{s}])
                % Get data to put in them, if it exists
                eval(['this_t = ages.n.t_' consts.sschemes{s} '(nindices(b));']);
                eval(['this_dt_int = ages.n.delt_int_' consts.sschemes{s} '(nindices(b));']);
                eval(['this_dt_ext = ages.n.delt_ext_' consts.sschemes{s} '(nindices(b));']);
                % Write age fields
                out = strcat(out,'<',tname,'>',sprintf('%0.0f',this_t),'</',tname,'>');
                out = strcat(out,'<',dt_int_name,'>',sprintf('%0.0f',this_dt_int),'</',dt_int_name,'>');
                out = strcat(out,'<',dt_ext_name,'>',sprintf('%0.0f',this_dt_ext),'</',dt_ext_name,'>',rr);               
            end
            
            % Also write fields with normalized nuclide concentrations
            NNname = strcat('Nnorm',this_nuclide(2:end),'_',consts.sschemes{s});
            dNN_int_name = strcat('delNnorm',this_nuclide(2:end),'_int_',consts.sschemes{s});
            dNN_ext_name = strcat('delNnorm',this_nuclide(2:end),'_ext_',consts.sschemes{s});
            
            if isfield(ages.n,['Nnorm_' consts.sschemes{s}])
                eval(['this_Nnorm = ages.n.Nnorm_' consts.sschemes{s} '(nindices(b));']);
                eval(['this_dNnorm_int = ages.n.delNnorm_' consts.sschemes{s} '(nindices(b));']);
                eval(['this_dNnorm_ext = ages.n.delNnorm_ext_' consts.sschemes{s} '(nindices(b));']);
                % Write age fields
                out = strcat(out,'<',NNname,'>',sprintf('%0.0f',this_Nnorm),'</',NNname,'>');
                out = strcat(out,'<',dNN_int_name,'>',sprintf('%0.0f',this_dNnorm_int),'</',dNN_int_name,'>');
                out = strcat(out,'<',dNN_ext_name,'>',sprintf('%0.0f',this_dNnorm_ext),'</',dNN_ext_name,'>',rr);
            end
        end
    end
            
    % Close result field
    out = [out '</exposureAgeResult>' rr];
    

end;

% Append elements for plot names if present

if ~isempty(splot_stubnames);
    for a = 1:length(splot_stubnames);
        out = [out '<ploturlstub>' consts.plotURL splot_stubnames{a} '</ploturlstub>'];
    end;
end;

% Append summary data if present

if ~isempty(summary);
    
    % Do some XML for summary data here. 
    summary_xml = ageSummaryToXML(summary);
    out = [out rr summary_xml rr];
    
    % Now link to camelplots
    if isfield(summary,'camels');
        out = [out '<camelploturlstubs><St>' consts.plotURL summary.camels{1} '</St>' rr];
        out = [out '<Lm>' consts.plotURL summary.camels{2} '</Lm>' rr];
        out = [out '<LSDn>' consts.plotURL summary.camels{3} '</LSDn>' rr '</camelploturlstubs>' rr];
    end;    
end;
    
% Also dump diagnostics

if isempty(ages.flags); 
    dstring = 'No diagnostics'; 
else
    dstring = strrep(ages.flags,'<br>','');
end;

out = [out '<diagnostics>' rr dstring  rr '</diagnostics>' rr];

% Also dump version info
% Should probably have a versions.summarize

if ~isempty(versions);
    out = [out '<version><wrapper>' versions.wrapper '</wrapper>' rr '<validate>' ...
        versions.validate '</validate>' rr '<get_age>' versions.get_age ...
        '</get_age>' rr '<muons>' versions.muons '</muons><consts>' consts.version ...
        '</consts></version>' rr];
end;

% Close 
out = [out '</calcs_v3_age_data>'];



