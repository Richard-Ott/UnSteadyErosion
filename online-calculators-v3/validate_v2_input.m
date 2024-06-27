function out = validate_v2_input(text_block,option)

% This reads, parses, and validates a block of text input in
% exposure-age-calculator v2 format. Returns a data structure with the
% relevant information in useful form. 
%
% One important thing about this is that it restandardizes nuclide
% concentrations internally. So returns nuclide concentrations on reference
% standard for all nuclides. This is 07KNSTD for Be-10 and KNSTD for Al-26.
% 
% This doesn't work with v2 calibration data, i.e. with the extra truet and
% deltruet at the end. 
%
% Modified 201703 to also process v2 form erosion rate input. The only
% difference there is that there is one less item in each line. 

% Deal with erosion flag, if present

if nargin < 2;
    eFlag = 0;
else
    if strcmp(option,'erosion');
        eFlag = 1;
    else
        error('validate_v2_input.m: unrecognized optional arg');
    end;
end;

% Parse into non-white-space elements

parsed_text = textscan(text_block,'%s');
parsed_text = parsed_text{1}; % Peculiarity of textscan

% Now a text array called parsed_text contains all the 
% separate items from the text block as array elements.
numitems = length(parsed_text);

% Here define the correct number of items per row -- 
% standard v2.2 input has 15 items for exposure age, 14 for erosion
if eFlag == 1;
    numcols = 14;
else
    numcols = 15; 
end;
    
% Load constants file. 
% Only necessary here for matching standard names. 
% Make sure to use correct file name. 
load consts_v3;

% Version

out.version_validate = 'validate_v2_input.m - 3.0';

%% ------------ DATA CHECKING AND LOADING FOR UNKNOWNS -----------------

% Initialize error flag

out.error = 0;

% check for correct number of items --

if mod(numitems,numcols) ~= 0;
	out.error = 1;
    out.message = 'validate_v2_input.m: Wrong total number of data elements in sample data block -- check for missing data, extra data, or white space within a single element';
    return;
end;

% if passed that, get number of samples -- 

out.numsamples = numitems./numcols;

% And preallocate sample data structure

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
out.s.yr = out.s.lat; % Year of sample collection; will have default


% Initialize nuclide data structure
out.n.index = []; % Field indexing nuclide measurement to sample number
out.n.nuclide = {}; % Nuclide/target identifier
out.n.N = []; % Properly standardized nuclide concentration
out.n.delN = []; % Same, uncertainty in

% Data checking loop. 
% Select a row, check the 15 strings to see if they are permissible, then 
% turn them into numbers. Check if the numbers are permissible. Finally,
% store them in an array for each input variable. 

for a = 1:out.numsamples;
	si = (a-1)*numcols; % starting index
	
	% 1. Sample name. 
	
	ino = 1; % ino = item number - sample name is item no. 1
	
	% test for length
	if length(parsed_text{si+ino}) > 24;
        out.error = 1;
		out.message = 'validate_v2_input.m: Sample data block - sample name more than 24 characters';
        return;
    end;
	
	% test for illegal characters in sample name
	% this allows letters, numbers, underscores, and dashes only. 
	if isempty(regexp(parsed_text{si+ino},'[^\w-]','once'));
		% pass, do assignment
		out.s.sample_name{a} = parsed_text{si+ino};
    else
		% fail
        out.error = 1;
		out.message = ['validate_v2_input.m: Sample data block - illegal characters in sample name - line ' int2str(a)];
        return;
	end;
		
	% 2. Latitude
	
	ino = 2;
	
	% illegal character test -- 
	% all numerical inputs may contain digits, ., e,E +, -. 
	if isempty(regexp(parsed_text{si+ino},'[^\d.eE+-]','once'));
		% pass
        % Convert to number
        temp = str2double(parsed_text{si+ino});
        % Check that worked
		if isnan(temp);
            out.error = 1;
			out.message = ['validate_v2_input.m: Sample data block - un-numericalizable latitude value - line ' int2str(a)];
            return;
		end;
		% test for bounds
		if (temp > 90) || (temp < -90);
            out.error = 1;
    		out.message = 'validate_v2_input.m: Sample data block - latitude out of bounds';
            return;
		end; 
		% Assign
		out.s.lat(a) = temp; clear temp;
    else
		% fail
        out.error = 1;
		out.message = ['validate_v2_input.m: Sample data block - illegal characters in latitude - line ' int2str(a)];
        return;
	end;

	
	% 3. Longitude
	
	ino = 3;
	
	if isempty(regexp(parsed_text{si+ino},'[^\d.eE+-]','once'));
		% pass
        % Convert to number
        temp = str2double(parsed_text{si+ino});
        % Check that worked
        if isnan(temp);
            out.error = 1;
			out.message = ['validate_v2_input.m: Sample data block - un-numericalizable longitude value - line ' int2str(a)];
            return;
		end;
		% test for bounds
		if (temp > 180) || (temp < -180);
            out.error = 1;
    		out.message = ['Sample data block - longitude out of bounds - line ' int2str(a)];
            return;
		end; 
		% Assign
		out.s.long(a) = temp; clear temp;
    else
		% fail
        out.error = 1;
		out.message = ['Sample data block - illegal characters in longitude - line ' int2str(a)];
        return;
	end;

	
	% 5. Elv/pressure flag -- get this first as it affects checks for (4)
	
	ino = 5;
	
	% must match one of three possible options
	if (strcmp(parsed_text{si+ino},'std') || strcmp(parsed_text{si+ino},'ant') || strcmp(parsed_text{si+ino},'pre') );
		% pass
		out.s.aa{a} = parsed_text{si+ino};
    else
		% fail
        out.error = 1;
		out.message = ['Sample data block - unknown elevation/pressure flag - line ' int2str(a)];
        return;
	end;

	% 4. Elv/pressure
	
	ino = 4;
	
	if isempty(regexp(parsed_text{si+ino},'[^\d.eE+-]','once'));
		% pass
        % Convert to number
        temp = str2double(parsed_text{si+ino});
        % Check that worked
        if isnan(temp);
            out.error = 1;
			out.message = ['validate_v2_input.m: Sample data block - un-numericalizable elevation/pressure value - line ' int2str(a)];
            return;
		end;
		% test for bounds
		if strcmp(out.s.aa{a},'std') || strcmp(out.s.aa{a},'ant')
			if (temp < -500);
                out.error = 1;
                out.message = ['validate_v2_input.m: Sample data block - elevation too low -- line ' int2str(a)];
                return;
			end; 
		elseif strcmp(out.s.aa{a},'pre')
			if (temp > 1080) || (temp < 0);
                out.error = 1;
    			out.message = ['validate_v2_input.m: Sample data block - pressure out of reasonable bounds -- line ' int2str(a)];
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
		out.message = ['validate_v2_input.m: Sample data block - illegal characters in elevation/pressure - line ' int2str(a)];
        return;
	end;
	
	% 6. Thickness
	
	ino = 6;
	
	if isempty(regexp(parsed_text{si+ino},'[^\d.eE+-]','once'));
		% Convert to number
        temp = str2double(parsed_text{si+ino});
        % Check that worked
        if isnan(temp);
            out.error = 1;
			out.message = ['validate_v2_input.m: Sample data block - un-numericalizable thickness value - line ' int2str(a)];
            return;
        end;
		% test for bounds
		if (temp < 0);
            out.error = 1;
    		out.message = ['validate_v2_input.m: Sample data block - thickness less than zero - line ' int2str(a)];
            return;
		end; 
		% Assign
		out.s.thick(a) = temp; clear temp;
    else
		% fail
        out.error = 1;
		out.message = ['validate_v2_input.m: Sample data block - illegal characters in thickness - line ' int2str(a)];
        return;
	end;
	
	% 7. Density
	
	ino = 7;
	
	if isempty(regexp(parsed_text{si+ino},'[^\d.eE+-]','once'));
		% pass
        % Convert to number
        temp = str2double(parsed_text{si+ino});
        % Check that worked
        if isnan(temp);
            out.error = 1;
			out.message = ['validate_v2_input.m: Sample data block - un-numericalizable density value - line ' int2str(a)];
            return;
		end;
		% test for bounds
		if (temp < 0);
            out.error = 1;
            out.message = ['validate_v2_input.m: Sample data block - density less than zero - line ' int2str(a)];
            return;
		end; 
		% Assign
		out.s.rho(a) = temp; clear temp;
    else
		% fail
        out.error = 1;
        out.message = ['validate_v2_input.m: Sample data block - illegal characters in density - line ' int2str(a)];
        return;
	end;
	
	% 8. Shielding
	
	ino = 8;
	
	if isempty(regexp(parsed_text{si+ino},'[^\d.eE+-]','once'));
		% pass
        % Convert to number
        temp = str2double(parsed_text{si+ino});
        % Check that worked
        if isnan(temp);
            out.error = 1;
			out.message = ['validate_v2_input.m: Sample data block - un-numericalizable shielding value - line ' int2str(a)];
            return;
		end;
		% test for bounds
		if (temp < 0) || (temp > 1);
            out.error = 1;
            out.message = ['validate_v2_input.m: Sample data block - shielding correction out of range - line ' int2str(a)];
            return;
		end; 
		% Assign
		out.s.othercorr(a) = temp; clear temp;
    else
		% fail
        out.error = 1;
        out.message = ['validate_v2_input.m: Sample data block - illegal characters in shielding correction - line ' int2str(a)];
        return;
	end;
	
	
	% 9. Erosion rate if processing exposure age data
	
    if eFlag == 0;
        ino = 9;

        if isempty(regexp(parsed_text{si+ino},'[^\d.eE+-]','once'));
            % pass
            % Convert to number
            temp = str2double(parsed_text{si+ino});
            % Check that worked
            if isnan(temp);
                out.error = 1;
                out.message = ['validate_v2_input.m: Sample data block - un-numericalizable erosion rate value - line ' int2str(a)];
                return;
            end;
            % test for bounds
            if (temp < 0);
                out.error = 1;
                out.message = ['validate_v2_input.m: Sample data block - erosion rate less than zero - line ' int2str(a)];
                return;
            end; 
            % Assign
            out.s.E(a) = temp; clear temp;
        else
            % fail
            out.error = 1;
            out.message = ['validate_v2_input.m: Sample data block - illegal characters in erosion rate - line ' int2str(a)];
            return;
        end;
    end; % close if-erosion loop
	
    % Also put default sample collection date in structure
    
    out.s.yr(a) = consts.default_yr;
    
    % Now define offset to distinguish between erosion/exposure age input
    eOffset = 0;
    if eFlag == 1; 
        eOffset = -1;
    end;
    
    % 12. Be-10 standardization -- get this first 
	
	ino = 12 + eOffset;
	
	% must match something in stds structure
    
    if any(strcmp(consts.be_stds_names,parsed_text{si+ino}));
        % pass, string is present in array of available stds
        current_std10 = parsed_text{si+ino};
    else
        % fail
        out.error = 1;
        out.message = ['validate_v2_input.m: Sample data block - unknown Be-10 reference standard identifier - line ' int2str(a)];
        return;
    end;
    
    
	% 10. N10
	
	ino = 10 + eOffset;
	
	if isempty(regexp(parsed_text{si+ino},'[^\d.eE+-]','once'));
		% pass
        % Convert to number
        temp = str2double(parsed_text{si+ino});
        % Check that worked
        if isnan(temp);
            out.error = 1;
			out.message = ['validate_v2_input.m: Sample data block - un-numericalizable Be-10 concentration - line ' int2str(a)];
            return;
		end;
		% test for bounds
		if (temp < 0);
            out.error = 1;
    		out.message = ['validate_v2_input.m: Sample data block - Be-10 concentration less than zero - line ' int2str(a)];
            return;
		end; 
        
        % Restandardize. Note that if "0" is entered as a placeholder this
        % will suppress calculation of Be-10 exposure age. 
        this_Be_std_no = find(strcmp(consts.be_stds_names,current_std10));
        this_Be_std_cf = consts.be_stds_cfs(this_Be_std_no);
        temp = temp.*this_Be_std_cf; 
        
        % If greater than zero, assign a measurement record.  
        % Assign
        if temp > 0;
            out.n.index(end+1) = a;
            out.n.nuclide{end+1} = 'N10quartz';
            out.n.N(end+1) = temp; clear temp;
            % Keep track of what is where
            out.isN10quartz(a) = true;
        else;
            out.isN10quartz(a) = false;
        end;
    else
		% fail
        out.error = 1;
		out.message = ['validate_v2_input.m: Sample data block - illegal characters in Be-10 concentration - line ' int2str(a)];
        return;
	end;
	
	% 11. delN10
	
	ino = 11 + eOffset;
	
	if isempty(regexp(parsed_text{si+ino},'[^\d.eE+-]','once'));
		% pass
        % Convert to number
        temp = str2double(parsed_text{si+ino});
        % Check that worked
        if isnan(temp);
            out.error = 1;
			out.message = ['validate_v2_input.m: Sample data block - un-numericalizable Be-10 uncertainty - line ' int2str(a)];
            return;
		end;
		% test for bounds
		if (temp < 0);
            out.error = 1;
            out.message = ['validate_v2_input.m: Sample data block - Be-10 uncertainty less than zero - line ' int2str(a)];
            return;
		end; 
    else
		% fail
        out.error = 1;
        out.message = ['validate_v2_input.m: Sample data block - illegal characters in Be-10 uncertainty - line ' int2str(a)];
        return;
	end;
    
    % Now we need to decide whether to record the uncertainty or not. 
    if out.isN10quartz(a);
        % Case a concentration was recorded, need an uncertainty
        if temp <= 0;
            % If a concentration was recorded, can't have zero uncertainty
            out.error = 1;
            out.message = ['validate_v2_input.m: If Be-10 concentration entered, uncertainty should be > 0 - line ' int2str(a)];
            return;
        end;
        % Now OK to record uncertainty
        % Restandardize
        temp = temp.*this_Be_std_cf; 
        % Record
        out.n.delN(end+1) = temp;
    else;
        % Case no concentration recorded - can't have uncertainty without
        % concentration
        if temp > 0;
            out.error = 1;
            out.message = ['validate_v2_input.m: Sample data block - Be-10 uncertainty has no corresponding concentration - line ' int2str(a)];
            return;   
        end;
    end;
        
	clear current_std10
    
    % 15. Al-26 standardization -- get this first 
	
	ino = 15 + eOffset;
	
    % must match something in stds structure
    
    if any(strcmp(consts.al_stds_names,parsed_text{si+ino}));
        % pass, string is present in array of available stds
        current_std26 = parsed_text{si+ino};
    else
        % fail
        out.error = 1;
        out.message = ['validate_v2_input.m: Sample data block - unknown Al-26 reference standard identifier - line ' int2str(a)];
        return;
    end;
    
	% 13. N26
	
	ino = 13 + eOffset;
	
	if isempty(regexp(parsed_text{si+ino},'[^\d.eE+-]','once'));
		% pass
        % Convert to number
        temp = str2double(parsed_text{si+ino});
        % Check that worked
        if isnan(temp);
            out.error = 1;
			out.message = ['validate_v2_input.m: Sample data block - un-numericalizable Al-26 concentration - line ' int2str(a)];
            return;
		end;
		% test for bounds
		if (temp < 0);
            out.error = 1;
    		out.message = ['validate_v2_input.m: Sample data block - Al-26 concentration less than zero - line ' int2str(a)];
            return;
		end; 
        
        % Restandardize
        this_Al_std_no = find(strcmp(consts.al_stds_names,current_std26));
        this_Al_std_cf = consts.al_stds_cfs(this_Al_std_no);
        temp = temp.*this_Al_std_cf; 
        
        % If greater than zero, assign a measurement record.  
        % Assign
        if temp > 0;
            out.n.index(end+1) = a;
            out.n.nuclide{end+1} = 'N26quartz';
            out.n.N(end+1) = temp; clear temp;
            % Keep track of what is where
            out.isN26quartz(a) = true;
        else;
            out.isN26quartz(a) = false;
        end;
        
        
    else
		% fail
        out.error = 1;
		out.message = ['validate_v2_input.m: Sample data block - illegal characters in Al-26 concentration - line ' int2str(a)];
        return;
	end;
	
	% 14. delN26
	
	ino = 14 + eOffset;
	
	if isempty(regexp(parsed_text{si+ino},'[^\d.eE+-]','once'));
		% pass
        % Convert to number
        temp = str2double(parsed_text{si+ino});
        % Check that worked
        if isnan(temp);
            out.error = 1;
			out.message = ['validate_v2_input.m: Sample data block - un-numericalizable Al-26 uncertainty - line ' int2str(a)];
            return;
		end;
		% test for bounds
		if (temp < 0);
            out.error = 1;
            out.message = ['validate_v2_input.m: Sample data block - Al-26 uncertainty less than zero - line ' int2str(a)];
            return;
		end; 
    else
		% fail
        out.error = 1;
        out.message = ['validate_v2_input.m: Sample data block - illegal characters in Al-26 uncertainty - line ' int2str(a)];
        return;
	end;
    
    % Now we need to decide whether to record the uncertainty or not. 
    if out.isN26quartz(a);
        % Case a concentration was recorded, need an uncertainty
        if temp <= 0;
            % If a concentration was recorded, can't have zero uncertainty
            out.error = 1;
            out.message = ['validate_v2_input.m: If Al-26 concentration entered, uncertainty should be > 0 - line ' int2str(a)];
            return;
        end;
        % Now OK to record uncertainty
        % Restandardize
        temp = temp.*this_Al_std_cf; 
        % Record
        out.n.delN(end+1) = temp;
    else
        % Case no concentration recorded - can't have uncertainty without
        % concentration
        if temp > 0;
            out.error = 1;
            out.message = ['validate_v2_input.m: Sample data block - Al-26 uncertainty has no corresponding concentration - line ' int2str(a)];
            return;   
        end;
    end;
    
    clear current_std26
    % Done with input data format/bounds checking
	
end;

% Check for nuclide concentration matches

% Find cases where neither were submitted
nodatacheck = find((~out.isN10quartz) & (~out.isN26quartz));
if ~isempty(nodatacheck);
    out.error = 1;
    out.message = [];
    for a = 1:length(nodatacheck);
        out.message = [out.message '<BR> validate_v2_input.m: Sample data block - need either Be-10 or Al-26 concentration - line ' int2str(nodatacheck(a))];
    end;
    return;
end;
    
% I forget what the below is for. 
out.numnuclides = length(find(out.isN10quartz)) + length(find(out.isN26quartz));

% Transpose data to match validate_v3_input. Not exactly sure why this is
% necessary. 

out.n.N = out.n.N';
out.n.delN = out.n.delN';
out.n.index = out.n.index';
out.n.nuclide = out.n.nuclide';


% ---------- DONE DATA CHECKING AND LOADING FOR UNKNOWNS -----------------


   
   
