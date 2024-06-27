function retstr = cal_input_v3(in)

% This is the wrapper script for the v3 production rate calibration code.
% Copyright: Greg Balco, Berkeley Geochronology Center
% December, 2016
% Not licensed for use or distribution. 

% Load constants

load consts_v3;

% Get in correct directory

if consts.isLocal == 0;
    % Case running on Linux web server
    cd /var/www/html/math/v3
else
    % Case on GB's Mac, no cd but declare globals for diagnostics
    global out;
end;

% Define version

versions.wrapper = '3.0.2';

% ------------------------- VALIDATE INPUT DATA --------------------

% Things to validate are as follows. 
% 1. in.text_block - input text
% 2. in.reportType - 'HTML' or 'XML' 
% 3. in.requesting_ip - ip address
% 4. in.plotFlag - 'yes' or 'no'

% Validate and parse main text block
% Must be v3 input format, error if v2

% Detect semicolons to determine where to pass data
if isempty(strfind(in.text_block,';'));
    % No semicolons; processing v2 input: error
    message = 'cal_input_v3.m: can''t process version 2 input format';
    if strcmp(in.reportType,'XML');
        retstr = dump_error_XML(message);
    elseif strcmp(in.reportType,'HTML');
        retstr = dump_error_HTML(message);
    end;
    return;
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

% Note we also need to make sure that we are only dealing with one nuclide.

nucs_present = unique(d.n.nuclide);
if length(nucs_present) > 1;
    % Case more than one nuclide in input data
    % Only allowable case is 'N3olivine' and 'N3pyroxene' together
    if length(nucs_present) == 2 && any(strcmp(nucs_present,'N3olivine')) && any(strcmp(nucs_present,'N3pyroxene'));
        % That is OK
    else
        % Not OK, error
        if strcmp(in.reportType,'XML');
            retstr = dump_error_XML('cal_input_v3.m: more than one nuclide/target in input data');
        else
            retstr = dump_error_HTML('cal_input_v3.m: more than one nuclide/target in input data');
        end;
        if consts.isLocal == 1;
            disp(d.message);
        end;
        return;
    end
end;

% Validate report type field

if ~(strcmp(in.reportType,'XML') || strcmp(in.reportType,'HTML'));
    message = ('cal_input_v3.m: Unrecognized value for report type field');
    retstr = dump_error_XML(message);
    if consts.isLocal == 1;
        % Dump to command window
        disp(message);
    end;
    return;
end;

% Field to determine whether plotting is needed for XML output

internalPlotFlag = 'no'; % default no plotting if no field present

if isfield(in,'plotFlag');
    if strcmp(in.plotFlag,'yes') || strcmp(in.plotFlag,'no')
        % Pass; set internal flag accordingly
        internalPlotFlag = in.plotFlag;
    else
        message = ('Unrecognized value for plot flag');
        retstr = dump_error_XML(message);
        if consts.isLocal == 1;
            % Dump to command window
            disp(message);
        end;
        return;
    end;
end;

% Calibration data set name

if isfield(in,'calibration_name');
    if length(in.calibration_name) > 70;
        message = 'cal_input_v3.m: too many characters in calibration data set name';
        if strcmp(in.reportType,'XML');
            retstr = dump_error_XML(message);
        elseif strcmp(in.reportType,'HTML');
            retstr = dump_error_HTML(message);
        end;
        return;
    end;
        
    if isempty(regexp(in.calibration_name,'[^\w\s-(),.]','once'));
		% pass, do assignment
		calibration_name = in.calibration_name;
    else
		% fail
        message = 'cal_input_v3.m: Not-allowed characters in calibration data set name';
        if strcmp(in.reportType,'XML');
            retstr = dump_error_XML(message);
        elseif strcmp(in.reportType,'HTML');
            retstr = dump_error_HTML(message);
        end;
        return;
    end;
else
    calibration_name = 'Unknown calibration data set name';
end;

% ------------------- DONE VALIDATING INPUT DATA ----------------------

% ------------ CALCULATE PRODUCTION RATES -----------------------------

% Now we should have a properly validated input data set. 
% Pass this to the age calculation script. 

% remember d is the validated input data block. 

control.cal = 1; 
control.resultType = 'long'; % Always all SF for this purpose

% Do the calculation

out = get_ages_v3(d,control,consts); % Note pass consts. 

% Record version info

versions.get_age = out.version_ages;
versions.muons = out.version_muons;
versions.validate = d.version_validate;
versions.consts = consts.version;

% Don't forget to make the trace string

trace_string = ['get_age - ' versions.get_age ' muons - ' versions.muons ' validate - ' versions.validate ' consts - ' versions.consts];

% --------------- DONE CALCULATING PRODUCTION RATES ---------------------

% ---------------------- WRITE TO LOG -----------------------------------

% Get IP address
this_ip = getenv('REMOTE_ADDR');
this_client = getenv('REMOTE_HOST');
if ~isempty(this_client);
    this_ip = this_client;
end;

if consts.isLocal == 0;  
    % Write to the log file. 
    for a = 1:out.numnuclides;
        temp_log_string = [char(out.s.sample_name{out.n.index(a)}) ' ' sprintf('%0.4f',out.s.lat(out.n.index(a))) ' ' sprintf('%0.4f',out.s.long(out.n.index(a)))];
        what_to_log = [this_ip ' ' in.reportType ' ' temp_log_string];
        write_to_log_v3('cal_input_v3.m',what_to_log);
    end;
end;
    
% ------------------ DONE WRITING TO LOG -----------------------------

% ----------------- DO SUMMARY STATS ---------------------------------

summary = summarize_cal_data(out);

% Label summary
summary.trace_string = trace_string;
summary.calibration_name = calibration_name;
summary.nuclide = out.n.nuclide{1}; 
% Note if both He-3/px and He-3/ol, this has to be dealt with later.


% ----------------- DONE WITH SUMMARY STATS ----------------------------

% ----------------------- DO SUMMARY PLOT ------------------------------

if strcmp(internalPlotFlag,'yes');
	% Make plot of all production rates
    Pplot_stubname = return_Pplot(out,summary);    
else
    % Case internalPlotFlag is 'no'; set empty variable for plot
    % that would have been made
    Pplot_stubname = [];
end;

% ----------------- DONE WITH SUMMARY PLOT ----------------------------

% -------------------- RETURN XML OR HTML -----------------------------

if strcmp(in.reportType,'XML');
    retstr = calToXML(out,versions,summary,Pplot_stubname);
elseif strcmp(in.reportType,'HTML');
    retstr = calToHTML(versions,summary,Pplot_stubname);
end;

% -------------------------------- DONE ------------------------------








