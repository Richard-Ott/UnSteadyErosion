function retstr = NofE_input_v3(in)

% This is the wrapper script for use of the erosion rate code to generate
% predicted N as a function of erosion rate. 
%
% The purpose of this script is to service some external piece of code that
% wants to invert something for an erosion rate, but doesn't want to do the
% production rate calculations. 
% 
% This is a lot simpler than the erosion rates wrapper script, because it
% doesn't have to invert for the erosion rate. 
%
% Copyright: Greg Balco, Berkeley Geochronology Center
% May 23, 2022
% Not licensed for use or distribution. 

% Load constants

load consts_v3;

% Get in correct directory

if consts.isLocal == 0
    % Case running on Linux web server
    cd /var/www/html/math/v3
else
    % Case most likely on GB's Mac
    % Do nothing
end

% Define version

versions.wrapper = '3.0-dev';

% ------------------------- VALIDATE INPUT DATA --------------------

% This only needs to validate in.text_block, because report type is always
% XML, there is no plotting, and resultType is always 'long.'   
% Also, this only has to validate v3 input data, v2 is not expected. 

d = validate_v3_input(in.text_block);         

% Return errors from validation as appropriate. Note always uses XML. 

if d.error
    retstr = dump_error_XML(d.message);
    if consts.isLocal == 1
        disp(d.message);
    end
    return;
end

% Return error if any inapplicable nuclides were submitted. 

if any(strcmp(d.n.nuclide,'N3quartz'))
    retstr = dump_error_XML('NofE_input_v3.m: Can''t compute erosion rates from He-3 in quartz.');
    if consts.isLocal == 1
        disp('NofE_input_v3.m: Can''t compute erosion rates from He-3 in quartz.');
    end
    return;
end


% ---------------- DONE VALIDATING INPUT DATA ----------------------

% ----------------------- CALCULATE NS  -----------------------------

% Now we should have a properly validated input data set (d). 
% Pass this to the erosion rate calculation script. 

Ns = get_NofE_v3(d,consts); % Note pass consts

% Assemble version info for reporting
versions.NofE = Ns.version_NofE;
versions.muons = Ns.version_muons;
versions.validate = d.version_validate;
versions.consts = consts.version;

% ------------------ DONE CALCULATING NS -----------------------------

% ------------------- WRITE TO LOG -----------------------------------

% Get IP address
this_ip = getenv('REMOTE_ADDR');

if consts.isLocal == 0
    % Write to the log file. 
    for a = 1:Ns.numnuclides
        temp_log_string = [char(Ns.s.sample_name{Ns.n.index(a)}) ' ' sprintf('%0.4f',Ns.s.lat(Ns.n.index(a))) ' ' sprintf('%0.4f',Ns.s.long(Ns.n.index(a)))];
        what_to_log = [this_ip ' XML ' temp_log_string];
        write_to_log_v3('NofE_input_v3.m',what_to_log);
    end
end
    
% ------------------ DONE WRITING TO LOG -----------------------------

% ------------------------- RETURN XML -------------------------------

retstr = NofEToXML(Ns,versions);

% -------------------------------- DONE ------------------------------








