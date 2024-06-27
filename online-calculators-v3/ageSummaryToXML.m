function s = ageSummaryToXML(summary)

% This code converts the output of summarize_ages.m to an XML string. 
%
% Greg Balco
% Berkeley Geochronology Center
% February, 2018.


%% begin

s = '<summary>';


f1 = fieldnames(summary);
% That returns field names for all nuclides. 

ok_fields = find(~strcmp(f1,'camels'));

for a = 1:length(ok_fields);
    % Each nuclide
    this_nuclide = getfield(summary,f1{ok_fields(a)});
    if ~isempty(this_nuclide.St);
        % there's something there
        s = [s '<' f1{ok_fields(a)} '>']; % open block for nuclide
        f2 = fieldnames(this_nuclide);
        % That should return all scaling methods
        for b = 1:length(f2);
            this_sf = getfield(this_nuclide,f2{b});
            s = [s '<' f2{b} '>']; % open block for scaling factor
            f3 = fieldnames(this_sf);
            for c = 1:length(f3);
                % Loop through all fields in scaling factor
                this_outer_field = getfield(this_sf,f3{c});
                if isstruct(this_outer_field);
                    % Case it's a structure itself
                    s = [s '<strip>'];
                    f4 = fieldnames(this_outer_field);
                    for d = 1:length(f4);
                        this_inner_field = getfield(this_outer_field,f4{d});
                        % These should all be scalars or vectors or char now
                        if length(this_inner_field) == 1;
                            s = [s '<' f4{d} '>' sprintf('%0.4e',this_inner_field) '</' f4{d} '>'];
                        end;
                    end;
                    s = [s '</strip>'];
                elseif ischar(this_outer_field)
                    % Case string
                    temp = strrep(this_outer_field,'<br>','...');
                    s = [s '<' f3{c} '>' temp '</' f3{c} '>'];
                elseif length(this_outer_field) == 1;
                    % Case scalar
                    s = [s '<' f3{c} '>' sprintf('%0.4e',this_outer_field) '</' f3{c} '>'];
                else
                    % Do nothing
                end;
            end;
            s = [s '</' f2{b} '>']; % close block for scaling factor
        end;
        s = [s '</' f1{ok_fields(a)} '>']; % close block for nuclide
    end;
end;

s = [s '</summary>'];


          
    
                
                
