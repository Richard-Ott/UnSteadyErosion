function retstr = erosion_input_v3(in)

% This is the wrapper script for the v3 erosion rate calculator.
% Copyright: Greg Balco, Berkeley Geochronology Center
% March, 2017
% Not licensed for use or distribution. 

% Load constants

load consts_v3;

% Get in correct directory

if consts.isLocal == 0;
    % Case running on Linux web server
    cd /var/www/html/math/v3
else
    % Case on GB's Mac, no cd but declare globals for diagnostics
    global erates;
end;

% Define version

versions.wrapper = '3.0';

% ------------------------- VALIDATE INPUT DATA --------------------

% Things to validate are as follows. 
% 1. in.text_block - input text
% 2. in.reportType - 'HTML' or 'XML'
% 4. in.plotFlag - 'yes' or 'no'
% 5. in.resultType - 'short' or 'long'

% Validate and parse main text block

% This has to correctly parse either v2 or v3 input data. 
% We recognize v3 input data by looking for semicolons, which is a bit
% hokey, but whatever. 

% Detect semicolons to determine where to pass data
if isempty(strfind(in.text_block,';'));
    % No semicolons; processing v2 input
    d = validate_v2_input(in.text_block,'erosion');
else
    % Semicolons present; processing v3 input
    % v3 erosion rate input is just the same as exposure age, but whatever
    % is in the erosion rate place is ignored.
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

% Return error if any nuclides that we can't use for erosion rate estimates
% were submitted. 

if any(strcmp(d.n.nuclide,'N3quartz'));
    if strcmp(in.reportType,'XML');
        retstr = dump_error_XML('Can''t compute erosion rates from He-3 in quartz.');
    else
        retstr = dump_error_HTML('Can''t compute erosion rates from He-3 in quartz.');
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

% Field to determine whether splotting is needed for XML output

internalPlotFlag = 'no'; % default no splotting

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

% ---------------- DONE VALIDATING INPUT DATA ----------------------

% ------------ CALCULATE EROSION RATES -----------------------------

% Now we should have a properly validated input data set (d). 
% Pass this to the erosion rate calculation script. 

erates = get_erates_v3(d,control,consts); % Note pass consts. 

% Assemble version info for reporting
versions.erates = erates.version_erates;
versions.muons = erates.version_muons;
versions.validate = d.version_validate;
versions.consts = consts.version;

% --------------- DONE CALCULATING EROSION RATES ---------------------

% ------------------- WRITE TO LOG -----------------------------------

% Get IP address
this_ip = getenv('REMOTE_ADDR');

if consts.isLocal == 0;  
    % Write to the log file. 
    for a = 1:erates.numnuclides;
        temp_log_string = [char(erates.s.sample_name{erates.n.index(a)}) ' ' sprintf('%0.4f',erates.s.lat(erates.n.index(a))) ' ' sprintf('%0.4f',erates.s.long(erates.n.index(a)))];
        what_to_log = [this_ip ' ' in.reportType ' ' temp_log_string];
        write_to_log_v3('erosion_input_v3.m',what_to_log);
    end;
end;
    
% ------------------ DONE WRITING TO LOG -----------------------------

% --------- MAKE NUCLIDE-NUCLIDE PLOTS (for qtz pairs) ---------------

% Note this only executes if internalPlotFlag is 'yes'. 

if strcmp(internalPlotFlag,'yes');

    % Available pairs are 26q/10q, 26q/21q, 10q/21q, 10px/3px, 10px/3ol,
    % (14/10 or 14/26 or 14/21) in that priority order
    % You can represent that as a matrix...
    %
    %       N3q    N3px    N3ol    N10q    N14q    N21q    N26q  N10px
    % N3q   No     No      No      No      No      No      No    No
    % N3px         No      No      No      No      No      No    Yes
    % N3ol                 No      No      No      No      No    Yes
    % N10q                         No      B1      Yes     Yes   No
    % N14q                                 No      B3      B2    No
    % N21q                                         No      Yes   No
    % N26q                                                 No    No
    % N10px                                                      No
    
    % He-3-in-quartz should never be involved. 
    
    % First we need to determine which pairs exist. 
    
    pairs = false(8,8);
    % Scan through all samples and determine which pairs they have.
    for a = 1:erates.numsamples;
        these_nuclides = erates.n.nuclide(erates.n.index == a);
        t = false(1,8);
        for b = 1:length(these_nuclides);
            t = t | strcmp(these_nuclides{b},consts.nuclides);
        end;
        tt =repmat(t,8,1);
        pairs = pairs | (tt & tt');
    end;
    % Remove diagonal and lower triangular
    pairs = triu(pairs) & (~diag(ones(1,8)));
    
    % False all N3px and N3ol...no, don't do this any more because N10px is
    % allowed
    %pairs([2 3],:) = false(2,7);
    %pairs(:,[2 3]) = false(7,2);
    
    % False all N3quartz
    pairs(1,:) = false(1,8);
    pairs(:,1) = false(8,1);

    % Prioritize only one involving C-14-quartz
    if pairs(4,5);
        pairs(5,[6 7]) = [false false];
    elseif pairs (5,7);
        pairs(5,6) = false;
    end;

    % Find indices
    [npairs1,npairs2] = find(pairs);

    splot_stubnames = {};

    % Now make plots
    for a = 1:length(npairs1);
        % A He retention plot should never be called in this script. 
        % Plotflag > 1 forces GMT plotting, use 1 for MATLAB plotting
        try
            splot_stubnames{end+1} = generic_splot(consts.nuclides{npairs1(a)},consts.nuclides{npairs2(a)},erates,2);
        catch
            splot_stubnames{end+1} = '';
        end;
    end;

    % ----------------- DONE MAKING NUCLIDE-NUCLIDE PLOTS ----------------    
else
    % This is case internalPlotFlag is 'no'; set dummy variable for plots
    % that would have been made
    splot_stubnames = {};
end;

% -------------------- RETURN HTML OR XML ----------------------------


if strcmp(in.reportType,'XML');
    retstr = eratesToXML(erates,versions,splot_stubnames);
else
    retstr = eratesToHTML(erates,versions,splot_stubnames);
end;


% -------------------------------- DONE ------------------------------








