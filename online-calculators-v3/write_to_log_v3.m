function write_to_log_v3(whoIsLogging,whatToLog)

% This writes to the v3 log file.
%
% Syntax: write_to_log_v3(whoIsLogging,whatToLog)
% both args are strings

str = sprintf('%02.0f ',clock);
str = [str ' ' whoIsLogging ' ' whatToLog];

fid = fopen('/var/www/html/scratch/log_v3.log','a');
fprintf(fid,'%s\n',str);
fclose(fid);