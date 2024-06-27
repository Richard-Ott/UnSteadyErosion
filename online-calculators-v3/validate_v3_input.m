function out = validate_v3_input(text_block)

% This parses and validates a block of text input data in exposure age
% calculator v3 form. Returns a usefully arranged data structure. 
%
% Supports v3 multiple-nuclide-analysis input scheme, which is described
% here:
%
% http://hess.ess.washington.edu/math/docs/v3/v3_input_explained.html
%
% Note that all restandardization of nuclide concentrations happens
% internally to this function. So the output data structure contains
% concentrations that are normalized to the reference standard values, i.e.
% 07KNSTD (Be-10), KNSTD (Al-26), CRONUS-A = 320 (Ne-21), CRONUS-P = 5.02
% (He-3). Later on, all production rates are referred to those stds.
%
% Note that this requires a constants file to do the restandardizations.
% See make_consts_v3.m.
%
% Aug 2015, added support for independent age (true_t) data lines. 
%
% Written by Greg Balco -- Berkeley Geochronology Center
% balcs@bgc.org
% 2014-2015
% 
% Copyright 2014-2015, Greg Balco and Berkeley Geochronology Center 
% All rights reserved
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License, version 2,
% as published by the Free Software Foundation (www.fsf.org).


%% Step 0. Load constants file. Only needed for restandardization info.

load consts_v3;

% Initialize error flag

out.error = 0;

% Version

out.version_validate = 'validate_v3_input.m - 3.0';

%% Step 1. Break up into individual lines using semicolons. 

remains = deblank(text_block);
% Also detab
remains = strrep(remains,char(9),' ');
k = 1;
while 1;
    [lines{k}, remains] = strtok(remains,';');
    if isempty(lines{k}); break; end;
    k = k + 1;
end;
lines = lines(1:(end-1)); % Remove the one that signaled the end

if isempty(lines);
    out.error = 1;
    out.message = 'validate_v3_input.m: Nothing in sample data block';
    return;
elseif length(lines) < 2;
    out.error = 1;
    out.message = 'validate_v3_input.m: Only one line of data - nothing to calculate';
    return;
end;

%% Step 2. Classify lines. There are three kinds: sample, nuclide, and 
% independent-age lines. 

ltype = zeros(size(lines));
for a = 1:length(lines);
    if any(strfind(lines{a},' Be-10 ')) || any(strfind(lines{a},' Al-26 ')) || any(strfind(lines{a},' Ne-21 ')) || any(strfind(lines{a},' He-3 ')) || any(strfind(lines{a},' C-14 ')) || any(strfind(lines{a},' Cl-36 '));
        % Case a line with nuclide concentrations in it, flagged by 0
    elseif any(strfind(lines{a},' true_t '))
        % Case a line with independent age information in it
        ltype(a) = 2;
    else
        % Case must be a sample data line, flag with 1. If something slips
        % through the cracks here it will be caught in processing of each
        % line. 
        ltype(a) = 1;
    end;
end;

s_lines = find(ltype == 1);
n_lines = find(ltype == 0);
t_lines = find(ltype == 2);
out.numsamples = length(s_lines);
out.numnuclides = length(n_lines);
out.numts = length(t_lines);

% Here decide if we are doing calibration data or not

if out.numts > 0 && (out.numts ~= out.numsamples);
    % Case there is not an independent age for each sample
    out.error = 1;
    out.message = 'validate_v3_input.m: Mismatch between number of sample and independent age lines: must be 1:1 correspondence';
    return;
end;

% Preallocate sample, nuclide data structures

out.s.sample_name = cell(out.numsamples,1);
out.s.aa = out.s.sample_name;
out.s.lat = zeros(out.numsamples,1);
out.s.long = out.s.lat;
out.s.elv = out.s.lat; 
out.s.pressure = out.s.lat;
out.s.thick = out.s.lat;
out.s.rho = out.s.lat;
out.s.othercorr = out.s.lat;
out.s.E = out.s.lat;
out.s.yr = out.s.lat; % year of sample collection

out.n.index = zeros(out.numnuclides,1); % Field indexing nuclide measurement to sample number
out.n.nuclide = cell(out.numnuclides,1); % Nuclide/target identifier
out.n.sample_name = out.n.nuclide;
out.n.N = out.n.index; % Properly standardized nuclide concentration
out.n.delN = out.n.index; % Same, uncertainty in

if length(out.numts) > 0;
    % If calibration data set, allocate truet. This should have same
    % indexing as samples.
    out.c.sample_name = cell(out.numts,1);
    out.c.site_name = out.c.sample_name;
    out.c.truet = zeros(out.numts,1);
    out.c.dtruet = out.c.truet;
    out.c.mint = out.c.truet; % Min ages default 0
    out.c.dmint = out.c.truet;
    out.c.maxt = out.c.truet + Inf; % Max ages default Inf
    out.c.dmaxt = out.c.truet;
    out.c.index = out.c.truet;
end;

%% Step 3. Process sample data lines. 

for a = 1:length(s_lines);
    % Get line
    thisline = lines{s_lines(a)};
    % Parse line
    
    remains = deblank(thisline);
    k = 1;
    while 1;
        [parsed_text{k}, remains] = strtok(remains);
        if isempty(parsed_text{k}); break; end;
        k = k + 1;
    end;
    
    % Check for length
    if length(parsed_text) ~= 11;
        out.error = 1;
        out.message = ['validate_v3_input.m: Wrong number of elements in sample line ' int2str(s_lines(a))];
        return;       
    end;
    
    % 1. Sample name, item 1
	
    % test for length
	if length(parsed_text{1}) > 32;
        out.error = 1;
		out.message = ['validate_v3_input.m: Sample name more than 32 characters - line ' int2str(s_lines(a))];
        return;
    end;
	
	% test for illegal characters in sample name
	% this allows letters, numbers, underscores, and dashes only. 
	if isempty(regexp(parsed_text{1},'[^\w-]','once'));
		% pass, do assignment
		out.s.sample_name{a} = parsed_text{1};
    else
		% fail
        out.error = 1;
		out.message = ['validate_v3_input.m: Illegal characters in sample name - line ' int2str(s_lines(a))];
        return;
	end;
    
    % 2. Latitude
	
	% illegal character test -- 
	% all numerical inputs may contain digits, ., e,E +, -. 
	if isempty(regexp(parsed_text{2},'[^\d.eE+-]','once'));
		% pass
        % Convert to number
        temp = str2double(parsed_text{2});
        % Check that worked
		if isnan(temp);
            out.error = 1;
			out.message = ['validate_v3_input.m: Un-numericalizable latitude value - line ' int2str(s_lines(a))];
            return;
		end;
		% test for bounds
		if (temp > 90) || (temp < -90);
            out.error = 1;
    		out.message = ['validate_v3_input.m: Latitude out of bounds - line ' int2str(s_lines(a))];
            return;
		end; 
		% Assign
		out.s.lat(a) = temp; clear temp;
    else
		% fail
        out.error = 1;
		out.message = ['validate_v3_input.m: Illegal characters in latitude - line ' int2str(s_lines(a))];
        return;
	end;
    
    % 3. Longitude

	if isempty(regexp(parsed_text{3},'[^\d.eE+-]','once'));
		% pass
        % Convert to number
        temp = str2double(parsed_text{3});
        % Check that worked
        if isnan(temp);
            out.error = 1;
			out.message = ['validate_v3_input.m: Sample data block - un-numericalizable longitude value - line ' int2str(s_lines(a))];
            return;
		end;
		% test for bounds
		if (temp > 180) || (temp < -180);
            out.error = 1;
    		out.message = ['validate_v3_input.m: Longitude out of bounds - line ' int2str(s_lines(a))];
            return;
		end; 
		% Assign
		out.s.long(a) = temp; clear temp;
    else
		% fail
        out.error = 1;
		out.message = ['validate_v3_input.m: Illegal characters in longitude - line ' int2str(s_lines(a))];
        return;
	end;
    
    % 5. Elv/pressure flag -- get this first as it affects checks for (4)
	
	% must match one of three possible options
	if (strcmp(parsed_text{5},'std') || strcmp(parsed_text{5},'ant') || strcmp(parsed_text{5},'pre') );
		% pass
		out.s.aa{a} = parsed_text{5};
    else
		% fail
        out.error = 1;
		out.message = ['validate_v3_input.m: Unknown elevation/pressure flag - line ' int2str(s_lines(a))];
        return;
	end;

    % 4. Elv/pressure
	
	if isempty(regexp(parsed_text{4},'[^\d.eE+-]','once'));
		% pass
        % Convert to number
        temp = str2double(parsed_text{4});
        % Check that worked
        if isnan(temp);
            out.error = 1;
			out.message = ['validate_v3_input.m: Un-numericalizable elevation/pressure value - line ' int2str(s_lines(a))];
            return;
		end;
		% test for bounds
		if strcmp(out.s.aa{a},'std') || strcmp(out.s.aa{a},'ant')
			if (temp < -500);
                out.error = 1;
                out.message = ['validate_v3_input.m: Elevation too low -- line ' int2str(s_lines(a))];
                return;
			end; 
		elseif strcmp(out.s.aa{a},'pre')
			if (temp > 1080) || (temp < 0);
                out.error = 1;
    			out.message = ['validate_v3_input.m: Pressure out of reasonable bounds -- line ' int2str(s_lines(a))];
                return;
			end;
		end;	
		
		% OK, pack either elv or pressure field
		if strcmp(out.s.aa{a},'std') || strcmp(out.s.aa{a},'ant')
			% store the elevation value
			out.s.elv(a) = temp;
			out.s.pressure(a) = -99;
        elseif strcmp(out.s.aa{a},'pre')
			% store the pressure value
			out.s.elv(a) = -99;
			out.s.pressure(a) = temp;
		end;
        clear temp;
    else
		% fail
        out.error = 1;
		out.message = ['validate_v3_input.m: Illegal characters in elevation/pressure - line ' int2str(s_lines(a))];
        return;
	end;
    
    % 6. Thickness
	
	if isempty(regexp(parsed_text{6},'[^\d.eE+-]','once'));
		% Convert to number
        temp = str2double(parsed_text{6});
        % Check that worked
        if isnan(temp);
            out.error = 1;
			out.message = ['validate_v3_input.m: Un-numericalizable thickness value - line ' int2str(s_lines(a))];
            return;
        end;
		% test for bounds
		if (temp < 0);
            out.error = 1;
    		out.message = ['validate_v3_input.m: Thickness less than zero - line ' int2str(s_lines(a))];
            return;
		end; 
		% Assign
		out.s.thick(a) = temp; clear temp;
    else
		% fail
        out.error = 1;
		out.message = ['validate_v3_input.m: Illegal characters in thickness - line ' int2str(s_lines(a))];
        return;
	end;
    
    % 7. Density
	
	if isempty(regexp(parsed_text{7},'[^\d.eE+-]','once'));
		% pass
        % Convert to number
        temp = str2double(parsed_text{7});
        % Check that worked
        if isnan(temp);
            out.error = 1;
			out.message = ['validate_v3_input.m: Un-numericalizable density value - line ' int2str(s_lines(a))];
            return;
		end;
		% test for bounds
		if (temp < 0);
            out.error = 1;
            out.message = ['validate_v3_input.m: Density less than zero - line ' int2str(s_lines(a))];
            return;
		end; 
		% Assign
		out.s.rho(a) = temp; clear temp;
    else
		% fail
        out.error = 1;
        out.message = ['validate_v3_input.m: Illegal characters in density - line ' int2str(s_lines(a))];
        return;
	end;
    
    % 8. Shielding
	
    if isempty(regexp(parsed_text{8},'[^\d.eE+-]','once'));
		% pass
        % Convert to number
        temp = str2double(parsed_text{8});
        % Check that worked
        if isnan(temp);
            out.error = 1;
			out.message = ['validate_v3_input.m: Un-numericalizable shielding value - line ' int2str(s_lines(a))];
            return;
		end;
		% test for bounds
		if (temp < 0) || (temp > 1);
            out.error = 1;
            out.message = ['validate_v3_input.m: Shielding correction out of range - line ' int2str(s_lines(a))];
            return;
		end; 
		% Assign
		out.s.othercorr(a) = temp; clear temp;
    else
		% fail
        out.error = 1;
        out.message = ['validate_v3_input.m: Illegal characters in shielding correction - line ' int2str(s_lines(a))];
        return;
	end;
    
    % 9. Erosion rate
	
	if isempty(regexp(parsed_text{9},'[^\d.eE+-]','once'));
		% pass
        % Convert to number
        temp = str2double(parsed_text{9});
        % Check that worked
        if isnan(temp);
            out.error = 1;
			out.message = ['validate_v3_input.m: Un-numericalizable erosion rate value - line ' int2str(s_lines(a))];
            return;
		end;
		% test for bounds
		if (temp < 0);
            out.error = 1;
            out.message = ['validate_v3_input.m: Erosion rate less than zero - line ' int2str(s_lines(a))];
            return;
		end; 
		% Assign
		out.s.E(a) = temp; clear temp;
    else
		% fail
        out.error = 1;
        out.message = ['validate_v3_input.m: Illegal characters in erosion rate - line ' int2str(s_lines(a))];
        return;
	end;
    
    % 10. Year of sample collection
	
	if isempty(regexp(parsed_text{10},'[^\d]','once')); % Only allow digits
		% pass
        % Convert to number
        temp = round(str2double(parsed_text{10}));
        
        % Check that worked
        if isnan(temp);
            out.error = 1;
			out.message = ['validate_v3_input.m: Un-numericalizable collection year value - line ' int2str(s_lines(a))];
            return;
		end;
		% various bounds tests and insert default for zeros
        todays_date = datevec(date);
		if (temp > todays_date(1));
            out.error = 1;
            out.message = ['validate_v3_input.m: Collection year is in the future - line ' int2str(s_lines(a))];
            return;
        end;
        
        if (temp == 0);
            % Set default
            temp = consts.default_yr; 
        end;
        
        if (temp < 1700);
            out.error = 1;
            out.message = ['validate_v3_input.m: Collection year out of reasonable bounds - line ' int2str(s_lines(a))];
        end; 
            
		% Assign
		out.s.yr(a) = temp; clear temp;
    else
		% fail
        out.error = 1;
        out.message = ['validate_v3_input.m: Illegal characters in collection year - line ' int2str(s_lines(a))];
        return;
	end;
    
    clear parsed_text;
end;

%% Step 4. Process nuclide concentration data lines. 
% Restandardization happens here. 

for a = 1:length(n_lines);
    % Parse
    % Get line
    thisline = lines{n_lines(a)};
    % Parse line
    remains = deblank(thisline);
    k = 1;
    while 1;
        [parsed_text{k}, remains] = strtok(remains);
        if isempty(parsed_text{k}); break; end;
        k = k + 1;
    end;
    
    % Validate sample name
    % 1. Sample name, item 1
	
    % test for length
	if length(parsed_text{1}) > 32;
        out.error = 1;
		out.message = ['validate_v3_input.m: Sample name more than 32 characters - line ' int2str(n_lines(a))];
        return;
    end;
	
	% test for illegal characters in sample name
	% this allows letters, numbers, underscores, and dashes only. 
	if isempty(regexp(parsed_text{1},'[^\w-]','once'));
		% pass, do assignment
		out.n.sample_name{a} = parsed_text{1};
    else
		% fail
        out.error = 1;
		out.message = ['validate_v3_input.m: Illegal characters in sample name - line ' int2str(n_lines(a))];
        return;
	end;
    
    % Determine which sample to match
    %out.s.sample_name;
    %out.n.sample_name{a};
    temp_index = find(strcmp(out.s.sample_name,out.n.sample_name{a}));
    
    % Error on fail
    if isempty(temp_index);
        out.error = 1;
        out.message = ['validate_v3_input.m: Can''t match line ' int2str(n_lines(a)) ' to sample'];
        return;
    end;
    % Error on match more than one
    if length(temp_index) > 1;
        out.error = 1;
        out.message = ['validate_v3_input.m: Line ' int2str(n_lines(a)) ' appears to match more than one sample'];
        return;
    end;
    % Pass, assign index number
    out.n.index(a) = temp_index;
    
    % Validate nuclide identifier
    if isempty(strcmp(parsed_text{2},{'Be-10','Al-26','Ne-21','Cl-36','He-3','C-14'}));
        out.error = 1;
        out.message = ['validate_v3_input.m: Can''t identify nuclide in line ' int2str(n_lines(a))];
        return;
    end;
    
    % Nuclide-specific actions
    if strcmp(parsed_text{2},'Be-10');
        
        % Case Be-10
        % Check for length
        if length(parsed_text) ~= 7;
            out.error = 1;
            out.message = ['validate_v3_input.m: Wrong number of elements in Be-10 line ' int2str(n_lines(a))];
            return;
        end;
        
        % Validate the mineral ID
        if any(strcmp(consts.ok_mins_10,parsed_text{3}));
            % pass, string is present in array of available stds, record
            out.n.nuclide{a} = ['N10' parsed_text{3}];
        else
            % fail
            out.error = 1;
            out.message = ['validate_v3_input.m: Unknown mineral for Be-10 measurement - line ' int2str(n_lines(a))];
            return;
        end;
        
        % Validate Be standard identifier
        if any(strcmp(consts.be_stds_names,parsed_text{6}));
            % pass, string is present in array of available stds
            current_std10 = parsed_text{6};
        else
            % fail
            out.error = 1;
            out.message = ['validate_v3_input.m: Unknown Be-10 standard identifier - line ' int2str(n_lines(a))];
            return;
        end;
        
        % N10
        
        if isempty(regexp(parsed_text{4},'[^\d.eE+-]','once'));
            % pass
            % Convert to number
            temp = str2double(parsed_text{4});
            % Check that worked
            if isnan(temp);
                out.error = 1;
                out.message = ['validate_v3_input.m: Un-numericalizable Be-10 concentration - line ' int2str(n_lines(a))];
                return;
            end;
            % test for bounds
            if (temp <= 0);
                out.error = 1;
                out.message = ['validate_v3_input.m: Be-10 concentration less than or equal to zero - line ' int2str(n_lines(a))];
                return;
            end; 

            % Restandardize.
            this_Be_std_no = find(strcmp(consts.be_stds_names,current_std10));
            this_Be_std_cf = consts.be_stds_cfs(this_Be_std_no);
            % Store
            out.n.N(a) = temp.*this_Be_std_cf; 
            clear temp;
        else
            % fail
            out.error = 1;
            out.message = ['validate_v3_input.m: Illegal characters in Be-10 concentration - line ' int2str(n_lines(a))];
            return;
        end;
        
        % delN10
        
        if isempty(regexp(parsed_text{5},'[^\d.eE+-]','once'));
            % pass
            % Convert to number
            temp = str2double(parsed_text{5});
            % Check that worked
            if isnan(temp);
                out.error = 1;
                out.message = ['validate_v3_input.m: Un-numericalizable Be-10 uncertainty - line ' int2str(n_lines(a))];
                return;
            end;
            % test for bounds
            if (temp <= 0);
                out.error = 1;
                out.message = ['validate_v3_input.m: Be-10 uncertainty less than or equal to zero - line ' int2str(n_lines(a))];
                return;
            end; 

            % Restandardize
            out.n.delN(a) = temp.*this_Be_std_cf; 
            clear temp;
        else
            % fail
            out.error = 1;
            out.message = ['validate_v3_input.m: Illegal characters in Be-10 uncertainty - line ' int2str(n_lines(a))];
            return;
        end;
     
    elseif strcmp(parsed_text{2},'Al-26');
        % Case Al-26
        % Check for length
        if length(parsed_text) ~= 7;
            out.error = 1;
            out.message = ['validate_v3_input.m: Wrong number of elements in Al-26 line ' int2str(n_lines(a))];
            return;
        end;
        
        % Validate the mineral ID
        if any(strcmp(consts.ok_mins_26,parsed_text{3}));
            % pass, string is present in array of available stds, record
            out.n.nuclide{a} = ['N26' parsed_text{3}];
        else
            % fail
            out.error = 1;
            out.message = ['validate_v3_input.m: Unknown mineral for Al-26 measurement - line ' int2str(n_lines(a))];
            return;
        end;
        
        % Validate Al standard identifier
        if any(strcmp(consts.al_stds_names,parsed_text{6}));
            % pass, string is present in array of available stds
            current_std26 = parsed_text{6};
        else
            % fail
            out.error = 1;
            out.message = ['validate_v3_input.m: Unknown Al-26 standard identifier - line ' int2str(n_lines(a))];
            return;
        end;
        
        % N26
        
        if isempty(regexp(parsed_text{4},'[^\d.eE+-]','once'));
            % pass
            % Convert to number
            temp = str2double(parsed_text{4});
            % Check that worked
            if isnan(temp);
                out.error = 1;
                out.message = ['validate_v3_input.m: Un-numericalizable Al-26 concentration - line ' int2str(n_lines(a))];
                return;
            end;
            % test for bounds
            if (temp <= 0);
                out.error = 1;
                out.message = ['validate_v3_input.m: Al-26 concentration less than or equal to zero - line ' int2str(n_lines(a))];
                return;
            end; 

            % Restandardize.
            this_Al_std_no = find(strcmp(consts.al_stds_names,current_std26));
            this_Al_std_cf = consts.al_stds_cfs(this_Al_std_no);
            % Store
            out.n.N(a) = temp.*this_Al_std_cf; 
            clear temp;
        else
            % fail
            out.error = 1;
            out.message = ['validate_v3_input.m: Illegal characters in Al-26 concentration - line ' int2str(n_lines(a))];
            return;
        end;
        
        % delN26
        
        if isempty(regexp(parsed_text{5},'[^\d.eE+-]','once'));
            % pass
            % Convert to number
            temp = str2double(parsed_text{5});
            % Check that worked
            if isnan(temp);
                out.error = 1;
                out.message = ['validate_v3_input.m: Un-numericalizable Al-26 uncertainty - line ' int2str(n_lines(a))];
                return;
            end;
            % test for bounds
            if (temp <= 0);
                out.error = 1;
                out.message = ['validate_v3_input.m: Al-26 uncertainty less than or equal to zero - line ' int2str(n_lines(a))];
                return;
            end; 

            % Restandardize
            out.n.delN(a) = temp.*this_Al_std_cf; 
            clear temp;
        else
            % fail
            out.error = 1;
            out.message = ['validate_v3_input.m: Illegal characters in Al-26 uncertainty - line ' int2str(n_lines(a))];
            return;
        end;
        
    elseif strcmp(parsed_text{2},'He-3');
        % Case He-3
        % Check for length
        if length(parsed_text) ~= 8;
            out.error = 1;
            out.message = ['validate_v3_input.m: Wrong number of elements in He-3 line ' int2str(n_lines(a))];
            return;
        end;
        
        % Validate He standard identifier
        if any(strcmp(consts.he_stds_names,parsed_text{6}));
            % pass, string is present in array of available stds
            current_std3 = parsed_text{6};
        else
            % fail
            out.error = 1;
            out.message = ['validate_v3_input.m: Unknown He-3 standard identifier - line ' int2str(n_lines(a))];
            return;
        end;
        
        % Validate He standard concentration
        
        if isempty(regexp(parsed_text{7},'[^\d.eE+-]','once'));
            % pass
            % Convert to number
            current_std3_val = str2double(parsed_text{7});
            % Check that worked
            if isnan(current_std3_val);
                out.error = 1;
                out.message = ['validate_v3_input.m: Un-numericalizable He-3 standard concentration value - line ' int2str(s_lines(a))];
                return;
            end;
            % test for bounds
            if (current_std3_val < 0);
                out.error = 1;
                out.message = ['validate_v3_input.m: He-3 standard concentration less than zero - line ' int2str(s_lines(a))];
                return;
            end; 
        else
            % fail
            out.error = 1;
            out.message = ['validate_v3_input.m: Illegal characters in He-3 standard concentration - line ' int2str(s_lines(a))];
            return;
        end;

        
        % N3
        
        if isempty(regexp(parsed_text{4},'[^\d.eE+-]','once'));
            % pass
            % Convert to number
            temp_N3 = str2double(parsed_text{4});
            % Check that worked
            if isnan(temp_N3);
                out.error = 1;
                out.message = ['validate_v3_input.m: Un-numericalizable He-3 concentration - line ' int2str(n_lines(a))];
                return;
            end;
            % test for bounds
            if (temp_N3 <= 0);
                out.error = 1;
                out.message = ['validate_v3_input.m: He-3 concentration less than or equal to zero - line ' int2str(n_lines(a))];
                return;
            end;      
        else
            % fail
            out.error = 1;
            out.message = ['validate_v3_input.m: Illegal characters in He-3 concentration - line ' int2str(n_lines(a))];
            return;
        end;
        
        % delN3
        
        if isempty(regexp(parsed_text{5},'[^\d.eE+-]','once'));
            % pass
            % Convert to number
            temp_delN3 = str2double(parsed_text{5});
            % Check that worked
            if isnan(temp_delN3);
                out.error = 1;
                out.message = ['validate_v3_input.m: Un-numericalizable He-3 uncertainty - line ' int2str(n_lines(a))];
                return;
            end;
            % test for bounds
            if (temp_delN3 <= 0);
                out.error = 1;
                out.message = ['validate_v3_input.m: He-3 uncertainty less than or equal to zero - line ' int2str(n_lines(a))];
                return;
            end; 
        else
            % fail
            out.error = 1;
            out.message = ['validate_v3_input.m: Illegal characters in He-3 uncertainty - line ' int2str(n_lines(a))];
            return;
        end;
        
        
        % Restandardize.
        if strcmp(current_std3,'NONE');
            % No restandardization
            out.n.N(a) = temp_N3;
            out.n.delN(a) = temp_delN3;
        else;
            % Restandardize to reference values
            this_He_std_no = find(strcmp(consts.he_stds_names,current_std3));
            this_He_std_ref = consts.he_stds_ref(this_He_std_no);
            out.n.N(a) = temp_N3.*(this_He_std_ref./current_std3_val); 
            out.n.delN(a) = temp_delN3.*(this_He_std_ref./current_std3_val); 
            clear temp_N3 temp_delN3 current_std3_val;
        end;
        
        % Finally, validate the mineral ID
        if any(strcmp(consts.ok_mins_3,parsed_text{3}));
            % pass, string is present in array of available stds, record
            out.n.nuclide{a} = ['N3' parsed_text{3}];
        else
            % fail
            out.error = 1;
            out.message = ['validate_v3_input.m: Unknown mineral for He-3 measurement - line ' int2str(n_lines(a))];
            return;
        end;
        
    elseif strcmp(parsed_text{2},'Ne-21');
        % Case Ne-21
        % Check for length
        if length(parsed_text) ~= 8;
            out.error = 1;
            out.message = ['validate_v3_input.m: Wrong number of elements in Ne-21 line ' int2str(n_lines(a))];
            return;
        end;
        
        % Validate Ne standard identifier
        if any(strcmp(consts.ne_stds_names,parsed_text{6}));
            % pass, string is present in array of available stds
            current_std21 = parsed_text{6};
        else
            % fail
            out.error = 1;
            out.message = ['validate_v3_input.m: Unknown Ne-21 standard identifier - line ' int2str(n_lines(a))];
            return;
        end;
        
        % Validate Ne standard concentration
        
        if isempty(regexp(parsed_text{7},'[^\d.eE+-]','once'));
            % pass
            % Convert to number
            current_std21_val = str2double(parsed_text{7});
            % Check that worked
            if isnan(current_std21_val);
                out.error = 1;
                out.message = ['validate_v3_input.m: Un-numericalizable Ne-21 standard concentration value - line ' int2str(s_lines(a))];
                return;
            end;
            % test for bounds
            if (current_std21_val < 0);
                out.error = 1;
                out.message = ['validate_v3_input.m: Ne-21 standard concentration less than zero - line ' int2str(s_lines(a))];
                return;
            end; 
        else
            % fail
            out.error = 1;
            out.message = ['validate_v3_input.m: Illegal characters in Ne-21 standard concentration - line ' int2str(s_lines(a))];
            return;
        end;

        
        % N21
        
        if isempty(regexp(parsed_text{4},'[^\d.eE+-]','once'));
            % pass
            % Convert to number
            temp_N21 = str2double(parsed_text{4});
            % Check that worked
            if isnan(temp_N21);
                out.error = 1;
                out.message = ['validate_v3_input.m: Un-numericalizable Ne-21 concentration - line ' int2str(n_lines(a))];
                return;
            end;
            % test for bounds
            if (temp_N21 <= 0);
                out.error = 1;
                out.message = ['validate_v3_input.m: Ne-21 concentration less than or equal to zero - line ' int2str(n_lines(a))];
                return;
            end;      
        else
            % fail
            out.error = 1;
            out.message = ['validate_v3_input.m: Illegal characters in Ne-21 concentration - line ' int2str(n_lines(a))];
            return;
        end;
        
        % delN21
        
        if isempty(regexp(parsed_text{5},'[^\d.eE+-]','once'));
            % pass
            % Convert to number
            temp_delN21 = str2double(parsed_text{5});
            % Check that worked
            if isnan(temp_delN21);
                out.error = 1;
                out.message = ['validate_v3_input.m: Un-numericalizable Ne-21 uncertainty - line ' int2str(n_lines(a))];
                return;
            end;
            % test for bounds
            if (temp_delN21 <= 0);
                out.error = 1;
                out.message = ['validate_v3_input.m: Ne-21 uncertainty less than or equal to zero - line ' int2str(n_lines(a))];
                return;
            end; 
        else
            % fail
            out.error = 1;
            out.message = ['validate_v3_input.m: Illegal characters in Ne-21 uncertainty - line ' int2str(n_lines(a))];
            return;
        end;
        
        
        % Restandardize.
        if strcmp(current_std21,'NONE');
            % No restandardization
            out.n.N(a) = temp_N21;
            out.n.delN(a) = temp_delN21;
        else;
            % Restandardize to reference values
            this_Ne_std_no = find(strcmp(consts.ne_stds_names,current_std21));
            this_Ne_std_ref = consts.ne_stds_ref(this_Ne_std_no);
            out.n.N(a) = temp_N21.*(this_Ne_std_ref./current_std21_val); 
            out.n.delN(a) = temp_delN21.*(this_Ne_std_ref./current_std21_val); 
            clear temp_N21 temp_delN21 current_std21_val;
        end;
        
        % Finally, validate the mineral ID
        if any(strcmp(consts.ok_mins_21,parsed_text{3}));
            % pass, string is present in array of available stds, record
            out.n.nuclide{a} = ['N21' parsed_text{3}];
        else
            % fail
            out.error = 1;
            out.message = ['validate_v3_input.m: Unknown mineral for Ne-21 measurement - line ' int2str(n_lines(a))];
            return;
        end;
        
    elseif strcmp(parsed_text{2},'C-14');
        % Case C-14
        
        % Check for length
        if length(parsed_text) ~= 6;
            out.error = 1;
            out.message = ['validate_v3_input.m: Wrong number of elements in C-14 line ' int2str(n_lines(a))];
            return;
        end;
        
        % N14
        
        if isempty(regexp(parsed_text{4},'[^\d.eE+-]','once'));
            % pass
            % Convert to number
            temp = str2double(parsed_text{4});
            % Check that worked
            if isnan(temp);
                out.error = 1;
                out.message = ['validate_v3_input.m: Un-numericalizable C-14 concentration - line ' int2str(n_lines(a))];
                return;
            end;
            % test for bounds
            if (temp <= 0);
                out.error = 1;
                out.message = ['validate_v3_input.m: C-14 concentration less than or equal to zero - line ' int2str(n_lines(a))];
                return;
            end;  
            % Assign
            out.n.N(a) = temp; clear temp;
        else
            % fail
            out.error = 1;
            out.message = ['validate_v3_input.m: Illegal characters in C-14 concentration - line ' int2str(n_lines(a))];
            return;
        end;
        
        % delN14
        
        if isempty(regexp(parsed_text{5},'[^\d.eE+-]','once'));
            % pass
            % Convert to number
            temp = str2double(parsed_text{5});
            % Check that worked
            if isnan(temp);
                out.error = 1;
                out.message = ['validate_v3_input.m: Un-numericalizable C-14 uncertainty - line ' int2str(n_lines(a))];
                return;
            end;
            % test for bounds
            if (temp <= 0);
                out.error = 1;
                out.message = ['validate_v3_input.m: C-14 uncertainty less than or equal to zero - line ' int2str(n_lines(a))];
                return;
            end; 
            % Assign
            out.n.delN(a) = temp; clear temp;
        else
            % fail
            out.error = 1;
            out.message = ['validate_v3_input.m: Illegal characters in C-14 uncertainty - line ' int2str(n_lines(a))];
            return;
        end;
        
        % Finally, validate the mineral ID
        if any(strcmp(consts.ok_mins_14,parsed_text{3}));
            % pass, string is present in array of available stds, record
            out.n.nuclide{a} = ['N14' parsed_text{3}];
        else
            % fail
            out.error = 1;
            out.message = ['validate_v3_input.m: Unknown mineral for C-14 measurement - line ' int2str(n_lines(a))];
            return;
        end;

    elseif strcmp(parsed_text{2},'Cl-36');
        % Case Cl-36, fail
        out.error = 1;
        out.message = ['validate_v3_input.m: Can''t do anything with Cl-36 data in line ' int2str(n_lines(a))];
        return;
        
    else
        % Unidentified nuclide, fail
        out.error = 1;
        out.message = ['validate_v3_input.m: Unrecognized nuclide in line ' int2str(n_lines(a))];
        return;
        
    end;
        
    clear parsed_text;
end;

%% Step 5. Validate independent age lines. 
for a = 1:length(t_lines);
    % Parse
    % Get line
    thisline = lines{t_lines(a)};
    % Parse line
    remains = deblank(thisline);
    k = 1;
    while 1;
        [parsed_text{k}, remains] = strtok(remains);
        if isempty(parsed_text{k}); break; end;
        k = k + 1;
    end;
    
    
    % Validate sample name
    % 1. Sample name, item 1
	
    % test for length
	if length(parsed_text{1}) > 32;
        out.error = 1;
		out.message = ['validate_v3_input.m: Sample name more than 32 characters - line ' int2str(n_lines(a))];
        return;
    end;
	
	% test for illegal characters in sample name
	% this allows letters, numbers, underscores, and dashes only. 
	if isempty(regexp(parsed_text{1},'[^\w-]','once'));
		% pass, do assignment
		out.c.sample_name{a} = parsed_text{1};
    else
		% fail
        out.error = 1;
		out.message = ['validate_v3_input.m: Illegal characters in sample name - line ' int2str(n_lines(a))];
        return;
	end;
    
    % Determine which sample to match
    temp_index = find(strcmp(out.s.sample_name,out.c.sample_name{a}));
    
    % Error on fail
    if isempty(temp_index);
        out.error = 1;
        out.message = ['validate_v3_input.m: Can''t match line ' int2str(t_lines(a)) ' to sample'];
        return;
    end;
    % Error on match more than one
    if length(temp_index) > 1;
        out.error = 1;
        out.message = ['validate_v3_input.m: Line ' int2str(t_lines(a)) ' appears to match more than one sample'];
        return;
    end;
    % Pass. temp_index should tell you where in the out.t structure to put
    % the data. Also record the index. 
    
    out.c.index(a) = temp_index;
    
    % 3. Site name
    
    % test for illegal characters in site name
	% this allows letters, numbers, underscores, and dashes only. 
	if isempty(regexp(parsed_text{3},'[^\w-]','once'));
		% pass, do assignment
		out.c.site_name{a} = parsed_text{3};
    else
		% fail
        out.error = 1;
		out.message = ['validate_v3_input.m: Illegal characters in site name - line ' int2str(n_lines(a))];
        return;
	end;
    
    % 4-7. Ages and del-ages. 
    b = 4;
    while 1;
        if ~isempty(parsed_text{b});
            % Case something there
            thistext = parsed_text{b};
            if (isempty(regexp(thistext,'[^\d]','once'))) || strcmp(thistext,'Inf'); % Only allow digits or 'Inf'. 
                % pass, convert to number
                temp = round(str2double(thistext));
                
                % Check that worked
                if isnan(temp);
                    out.error = 1;
                    out.message = ['validate_v3_input.m: Un-numericalizable independent age/uncert - line ' int2str(s_lines(a)) ' field ' int2str(b)];
                return;
                end;
                
                if (temp < 0);
                    out.error = 1;
                    out.message = ['validate_v3_input.m: Independent age/uncert less than zero - line ' int2str(s_lines(a)) ' field ' int2str(b)];
                end; 
                % Assign
                temp_t(b-3) = temp; clear temp;
            else
                % fail
                out.error = 1;
                out.message = ['validate_v3_input.m: Illegal characters in independent age/uncert - line ' int2str(s_lines(a)) ' field ' int2str(b)];
                return;
            end;  
            
        else
            % Case is empty
            break;
        end;
        
        if b > 7;
            break;
        end;
        b = b + 1;
    end;

	if length(temp_t) == 2; 
        % Case exact age
        out.c.truet(temp_index) = temp_t(1);
        out.c.dtruet(temp_index) = temp_t(2);
    elseif length(temp_t) == 4;
        % Case min/max age
        out.c.mint(temp_index) = temp_t(1);
        out.c.dmint(temp_index) = temp_t(2);
        out.c.maxt(temp_index) = temp_t(3);
        out.c.dmaxt(temp_index) = temp_t(4);
    else
        out.error = 1;
        out.message = ['validate_v3_input.m: Wrong number of fields in independent age/uncert - line ' int2str(s_lines(a))];
        return;
    end;
    clear temp_t;
    clear temp_index;
end;


%% Step 6. Various match checks

% Sample with no nuclide measurements
for a = 1:out.numsamples;
    if isempty(find(out.n.index == a));
        % Case can't find any nuclide measurements indexed to this sample
        out.error = 1;
        out.message = ['validate_v3_input.m: Can''t match any measurements to sample in line ' int2str(s_lines(a))];
        return;
    end;
end;

% If calibration data set, should be one and only one true age line for
% each sample.

if isfield(out,'c');
    scindex = sort(out.c.index);
    if any(~(scindex == (1:length(out.c.index))'));
        out.error = 1;
        out.message = ['validate_v3_input.m: Mismatch between sample and calibration lines?'];
        return;
    end;
end;





%% Done. 






