function retstr = age_input_v3(in)

% This is the wrapper script for the v3 exposure age calculator.
% Copyright: Greg Balco, Berkeley Geochronology Center
% February, 2015
% Not licensed for use or distribution. 

% Load constants

load consts_v3;

% Get in correct directory

if consts.isLocal == 0;
    % Case running on Linux web server
    cd /var/www/html/math/v3
else
    % Case on GB's Mac, no cd but declare globals for diagnostics
    global ages;
end;

% Define version

versions.wrapper = '3.0.2';


% ------------------------- VALIDATE INPUT DATA --------------------

% Things to validate are as follows. 
% 1. in.text_block - input text
% 2. in.reportType - 'HTML' or 'XML'
% 3. in.requesting_ip - ip address
% 4. in.plotFlag - 'yes' or 'no'
% 5. in.resultType - 'short' or 'long'
% 6. in.summary - 'yes' or 'no'

% Validate and parse main text block

% This has to correctly parse either v2 or v3 input data. 
% We recognize v3 input data by looking for semicolons, which is a bit
% hokey, but whatever. 

% Detect semicolons to determine where to pass data
if isempty(strfind(in.text_block,';'));
    % No semicolons; processing v2 input
    d = validate_v2_input(in.text_block);
else
    % Semicolons present; processing v3 input
    d = validate_v3_input(in.text_block);         
end;

% Return errors from validation as appropriate

if d.error;
    if strcmp(in.reportType,'XML');
        retstr = dump_error_XML(d.message);
    else
        retstr = dump_error_HTML(d.message);
    end;
    if consts.isLocal == 1;
        disp(d.message);
    end;
    return;
end;

% Validate report type field

if ~(strcmp(in.reportType,'XML') || strcmp(in.reportType,'HTML'));
    message = ('Unrecognized value for report type field');
    retstr = dump_error_XML(message);
    if consts.isLocal == 1;
        % Dump to command window
        disp(message);
    end;
    return;
end;

% Field to determine whether plotting is needed for XML output

internalPlotFlag = 'no'; % default no plotting

if isfield(in,'plotFlag');
    if strcmp(in.plotFlag,'yes') || strcmp(in.plotFlag,'no')
        % Pass; set internal flag accordingly
        internalPlotFlag = in.plotFlag;
    else
        message = ('Unrecognized value for XML plot field');
        retstr = dump_error_XML(message);
        if consts.isLocal == 1;
            % Dump to command window
            disp(message);
        end;
        return;
    end;
end;

% Field to determine whether summary stats and camelplot are wanted

internalSummaryFlag = 'no'; % default no summary

if isfield(in,'summary');
    if strcmp(in.summary,'yes')
        % Pass; set internal flag accordingly
        internalSummaryFlag = 'yes'; 
    elseif strcmp(in.summary,'no')
        % do nothing, already set
    else
        message = ('Unrecognized value for summary switch field');
        retstr = dump_error_XML(message);
        if consts.isLocal == 1;
            % Dump to command window
            disp(message);
        end;
        return;
    end;
end;
      
% Field to determine whether to calculate short or long output
% short is just St scaling; long is other schemes

control.resultType = 'short'; % default

if isfield(in,'resultType');
    if strcmp(in.resultType,'short') || strcmp(in.resultType,'long')
        % Pass
        control.resultType = in.resultType;
    else
        message = ('Unrecognized value for result type field');
        retstr = dump_error_XML(message);
        if consts.isLocal == 1;
            % Dump to command window
            disp(message);
        end;
        return;
    end;
end;

% Data checking for calibration data

calib_trace.calibration_name = 'Default calibration data set'; % Default if nothing passed
calib_trace.trace_string = [];
calib_trace.nuclide = [];

if isfield(in,'trace_string');
    if isempty(regexp(in.trace_string,'[^\w\s-(),.:=]','once'));
		% pass, do assignment
		calib_trace.trace_string = in.trace_string;
    else
		% fail
        message = ['cal_input_v3.m: bad characters in trace string'];
        if strcmp(in.reportType,'XML');
            retstr = dump_error_XML(message);
        elseif strcmp(in.reportType,'HTML');
            retstr = dump_error_HTML(message);
        end;
        return;
	end;
end;

if isfield(in,'calib_name');
    if isempty(regexp(in.calib_name,'[^\w\s-(),.:]','once'));
		% pass
        calib_trace.calibration_name = in.calib_name;
    else
		% fail
        message = ['cal_input_v3.m: bad characters in calibration data set name'];
        if strcmp(in.reportType,'XML');
            retstr = dump_error_XML(message);
        elseif strcmp(in.reportType,'HTML');
            retstr = dump_error_HTML(message);
        end;
        return;
	end;
end;

if isfield(in,'nuclide_name');
    if any(strcmp(in.nuclide_name,consts.nuclides));
		% pass
        nuclide_index = find(strcmp(in.nuclide_name,consts.nuclides));
        if nuclide_index == 2 || nuclide_index == 3;
            calib_trace.nuclide = 'He-3 (px/ol)';
        else
            calib_trace.nuclide = consts.properName{nuclide_index};
        end;
    else
		% fail
        message = ['cal_input_v3.m: unknown name for calibrated nuclide'];
        if strcmp(in.reportType,'XML');
            retstr = dump_error_XML(message);
        elseif strcmp(in.reportType,'HTML');
            retstr = dump_error_HTML(message);
        end;
        return;
	end;
end;

% Loop through the scaling factors and extract stuff
P = [];
sfs = {'St','Lm','LSDn'};
items = {'P_','delP_'};
for sf = 1:length(sfs);
    for item = 1:length(items);
        % Get the item
        eval(['ok = isfield(in,''' [items{item} sfs{sf}] ''');']);
        if ok;
            eval(['this_item_str = in.' [items{item} sfs{sf}] ';']);
            % Check  - allow digits and decimal point
            if isempty(regexp(this_item_str,'[^\d.]','once'));
                % pass
                % test for bounds
                if (str2double(this_item_str) < 0);
                    error(['age_input_v3.m: ' [items{item} sfs{sf}] ' less than zero']);
                end; 
                % cover other eventualities
                if isnan(str2double(this_item_str));
                    error(['age_input_v3.m - ' [items{item} sfs{sf}] ' can''t be numericalized']);
                end;
                % Pass all - assign
                eval(['P.' [items{item} sfs{sf}] ' = str2double(this_item_str);']);
            else
                % fail
                message = ['age_input_v3.m: inappropriate characters in ' [items{item} sfs{sf}]];
                if strcmp(in.reportType,'XML');
                    retstr = dump_error_XML(message);
                elseif strcmp(in.reportType,'HTML');
                    retstr = dump_error_HTML(message);
                end;
                return;    
            end;       
        end;
    end;
end;

% ---------------- DONE VALIDATING INPUT DATA ----------------------

% ------------ CALCULATE EXPOSURE AGES -----------------------------

% Now we should have a properly validated input data set. 
% Pass this to the age calculation script. 

% remember d is the validated input data block. 

control.cal = 0; % This is the code for exposure ages. 
% There is another simpler one for calibration. 

% If a calibrated production rate was passed, must insert these values into
% the consts structure. 



if ~isempty(P);
    % Check we're not using long form results if we only passed calibration
    % data
    % for St. 
    if ~isfield(P,'P_Lm') || ~isfield(P,'P_LSDn');
        control.resultType = 'short';
    end;
    % Assignments
    % Note that if the calibration is for He-3/px or He-3/ol, it gets put
    % in both of those slots. 
    if nuclide_index == 2 || nuclide_index == 3;
        nuclide_index = [2 3];
    end;
    if isfield(P,'P_St');
        consts.refP_St(nuclide_index) = P.P_St.*ones(size(nuclide_index));
        consts.delrefP_St(nuclide_index) = P.delP_St.*ones(size(nuclide_index));
    end;
    if isfield(P,'P_Lm');
        consts.refP_Lm(nuclide_index) = P.P_Lm.*ones(size(nuclide_index));
        consts.delrefP_Lm(nuclide_index) = P.delP_Lm.*ones(size(nuclide_index));
    end;
    if isfield(P,'P_LSDn');
        consts.refP_LSDn(nuclide_index) = P.P_LSDn.*ones(size(nuclide_index));
        consts.delrefP_LSDn(nuclide_index) = P.delP_LSDn.*ones(size(nuclide_index));
    end;
end;

% Now calculate ages. 

ages = get_ages_v3(d,control,consts); % Note pass consts. 

versions.get_age = ages.version_ages;
versions.muons = ages.version_muons;
versions.validate = d.version_validate;
versions.consts = consts.version;

% --------------- DONE CALCULATING EXPOSURE AGES ---------------------

% ------------------- WRITE TO LOG -----------------------------------

% Get IP address
this_ip = getenv('REMOTE_ADDR');

if consts.isLocal == 0;  
    % Write to the log file. 
    for a = 1:ages.numnuclides;
        temp_log_string = [char(ages.s.sample_name{ages.n.index(a)}) ' ' sprintf('%0.4f',ages.s.lat(ages.n.index(a))) ' ' sprintf('%0.4f',ages.s.long(ages.n.index(a)))];
        what_to_log = [this_ip ' ' in.reportType ' ' temp_log_string];
        write_to_log_v3('age_input_v3.m',what_to_log);
    end;
end;
    
% ------------------ DONE WRITING TO LOG -----------------------------

% --------- MAKE NUCLIDE-NUCLIDE PLOTS (for qtz pairs) ---------------

% Note this only executes if internalPlotFlag is 'yes'. 

if strcmp(internalPlotFlag,'yes');

    % Available pairs are 26q/10q, 26q/21q, 10q/21q, 10px/3px, 10px/3ol
    % (3q/10q or 3q/26q or 3q/14q or 3q/21q) in that priority order 
    % (14q/10q or 14q/26q or 14q/21q) in that priority order
    % You can represent that as a matrix...
    %
    %       N3q    N3px    N3ol    N10q    N14q    N21q    N26q   N10px
    % N3q   No     No      No      A1      A4      A3      A2     No
    % N3px         No      No      No      No      No      No     Yes
    % N3ol                 No      No      No      No      No     Yes
    % N10q                         No      B1      Yes     Yes    No
    % N14q                                 No      B3      B2     No
    % N21q                                         No      Yes    No
    % N26q                                                 No     No
    % N10px                                                       No
    
    % Note: Sometimes 3q/14q doesn't make any sense. Should probably make that
    % always false. 
    
    % First we need to determine which pairs exist. 
    
    pairs = false(8,8);
    % Scan through all samples and determine which pairs they have.
    for a = 1:ages.numsamples;
        these_nuclides = ages.n.nuclide(ages.n.index == a);
        t = false(1,8);
        for b = 1:length(these_nuclides);
            t = t | strcmp(these_nuclides{b},consts.nuclides);
        end;
        tt =repmat(t,8,1);
        pairs = pairs | (tt & tt');
    end;
    % Remove diagonal and lower triangular
    pairs = triu(pairs) & (~diag(ones(1,8)));
    % Remove N3px and N3ol...not any more because they may go with N10px
    % pairs([2 3],:) = false(2,7);
    % pairs(:,[2 3]) = false(7,2);

    % Prioritize only one involving He-3-quartz. Note that this is not
    % completely inclusive - if one sample has 3/10 and another has 3/26, you
    % won't get to see the 3/26 plot. 
    if pairs(1,4) % if N3q and N10q
        pairs(1,5:7) = [false false false]; % set (3,14), (3,21) and (3,26) off
    elseif pairs (1,7) % if N3q and N26q
        pairs(1,[5 6]) = [false false]; % set (3,14) and (3,21) off
    elseif pairs (1,6) % if N3q and N21
        pairs(1,5) = false; % turn off (3,14)
    end

    % Prioritize only one involving C-14-quartz
    if pairs(4,5)
        pairs(5,[6 7]) = [false false];
    elseif pairs (5,7)
        pairs(5,6) = false;
    end

    % Find indices
    [npairs1,npairs2] = find(pairs);

    splot_stubnames = {};

    % Now make plots
    for a = 1:length(npairs1)
        if strcmp(consts.nuclides{npairs1(a)},'N3quartz') || strcmp(consts.nuclides{npairs2(a)},'N3quartz')
            % Note: must fix this to check that there are some nonzero ages
            % for whatever nuclide is not He-3. 
            % if any(ages.n.t_St(find(ages.n.nindex == npairs1(a))) > 0) & any(ages.n.t_St(find(ages.n.nindex == npairs2(a))) > 0)
            % Does that do it? 
            % Plotflag > 1 forces GMT plotting, use 1 for MATLAB plotting
            try
                splot_stubnames{end+1} = retention_plot(consts.nuclides{npairs1(a)},consts.nuclides{npairs2(a)},ages,2);
            catch
                splot_stubnames{end+1} = '';
            end
        else
            % Plotflag > 1 forces GMT plotting, use 1 for MATLAB plotting
            try
                splot_stubnames{end+1} = generic_splot(consts.nuclides{npairs1(a)},consts.nuclides{npairs2(a)},ages,2);
            catch
                splot_stubnames{end+1} = '';
            end
        end
    end

    % ----------------- DONE MAKING NUCLIDE-NUCLIDE PLOTS ----------------
    
    % -------------- MAKE C-14 SATURATION PLOT ---------------------------
    % This only happens if (i) we have C-14 data, and (ii) we are in
    % Antarctica. 
    % It just appends it to the other splots. 
    
    if any(strcmp(ages.n.nuclide,'N14quartz'))
        if any(ages.s.lat > -60)
            % Case some samples not in Antarctica, do nothing
        else
            % Case all samples in Antarctica, do something
            splot_stubnames{end+1} = return_sat14plot(ages,consts);
        end
    end

    % -------------------- MAKE TIME-DEPENDENCE PLOTS --------------------

    if d.numsamples == 1;
        % Make time-dependence plot
        tdplot_stubname = []; % Placeholder
    else
        tdplot_stubname = [];
    end;

    % --------------- DONE MAKING TIME-DEPENDENCE PLOTS ------------------

    % ---------------------- MAKE SF-COMPARISON PLOT ---------------------

    if d.numsamples == 1;
        % Make scaling scheme comparison plot
        compplot_stubname = []; % Placeholder
    else
        compplot_stubname = [];
    end;

    % ------------------- DONE MAKING SF-COMPARISON PLOT -----------------
    
else
    % This is case internalPlotFlag is 'no'; set dummy variables for plots
    % that would have been made
    splot_stubnames = {};
    tdplot_stubname = [];
    compplot_stubname = [];
end;

% ----------------- SUMMARIZE DATA AND MAKE CAMELPLOT IF NEEDED ----------

if strcmp(internalSummaryFlag,'yes');
    % Calculate stats for each nuclide separately as well as all nuclides,
    % all samples. This also makes camelplots. 
    if strcmp(internalPlotFlag,'yes');
        summary = summarize_ages(ages,1);
    else
        summary = summarize_ages(ages,0);
    end;
else
    summary = [];    
end;

% ------------------- DONE WITH SUMMARY AND CAMELPLOT --------------------

% -------------------- RETURN HTML OR XML ----------------------------


if strcmp(in.reportType,'XML'); 
    retstr = agesToXML(ages,versions,splot_stubnames,summary);
else
    retstr = agesToHTML(ages,versions,splot_stubnames,tdplot_stubname,compplot_stubname,summary,calib_trace);
end;


% -------------------------------- DONE ------------------------------








