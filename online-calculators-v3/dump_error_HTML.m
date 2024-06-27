function out = dump_error_HTML(instr)

% This just dumps an error string in an intelligible HTML form.

% Assemble HTML return string

out = ['<html><head><meta http-equiv="content-type" content="text/html;charset=utf-8" />' ...
    '<title>v3: error</title></head>' ...
    '<body>'];

out = [out ' <table class=standard width=1000><tr><td><hr></td><tr>' ...
    '<tr><td>' instr '</td></tr>' ...
    '<tr><td><hr></td></tr></table>'];

out = [out '</body></html>'];