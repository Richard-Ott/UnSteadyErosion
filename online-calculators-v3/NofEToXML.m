function out = NofEToXML(Ns,versions)

% This ingests the age results structure from get_NofE_v3 and spits out an
% XML result string. 

load consts_v3; % Needed to define scaling scheme loop

dstring = [];

out =  '<calcs_v3_NofE_data>';

rr = '';

for a = 1:Ns.numsamples
    % Which measurements belong to this sample
    nindices = find(Ns.n.index == a);
    % Open results field
    out = strcat(out,'<NofEResult>',rr);
    
    % Report sample name 
    out = strcat(out,'<sample_name>',Ns.s.sample_name{a},'</sample_name>',rr);
    % Indicate erosion rate for each sample for QC
    out = strcat(out,'<erosion_rate_cm_yr>',sprintf('%0.3e',Ns.s.E(a)),'</erosion_rate_cm_yr>');

    % Loop through nuclide-scaling method possibilities
    for b = 1:length(nindices)
        this_nuclide = char(Ns.n.nuclide(nindices(b)));
        % Wrap result for each nuclide
        out = strcat(out,'<nuclide_result><nuclide>',this_nuclide,'</nuclide>');
        % Loop through scaling methods
        for s = 1:length(consts.sschemes)
            % Create XML field names
            Nname = strcat('Npred',consts.sschemes{s});  
            if isfield(Ns.n,['Npred_' consts.sschemes{s}])
                % Get data to put in them, if it exists
                eval(['this_Npred = Ns.n.Npred_' consts.sschemes{s} '(nindices(b));']);
                % Write erosion rate fields
                out = strcat(out,'<',Nname,'>',sprintf('%0.4e',this_Npred),'</',Nname,'>');
            end          
        end
        % End wrap
        out = strcat(out,'</nuclide_result>');
    end
            
    % Close result field
    out = [out '</NofEResult>' rr];    

end
  
% Also dump diagnostics

if isempty(Ns.flags)
    dstring = 'No diagnostics'; 
else
    dstring = strrep(Ns.flags,'<br>','');
end

out = [out '<diagnostics>' rr dstring  rr '</diagnostics>' rr];

% Also dump version info

out = [out '<version><wrapper>' versions.wrapper '</wrapper><validate>' ...
    versions.validate '</validate><erates>' versions.NofE ...
    '</erates><muons>' versions.muons '</muons><consts>' consts.version ...
    '</consts></version>' rr];

% Close 
out = [out '</calcs_v3_NofE_data>'];



