function outstr = dump_error_XML(errorString)

% the following is scheme to dump errors in appropriate XML format

outstr = '<calcs_v3_age_data>';

% Put in some zero results to prevent potential failure of results display code
outstr = [outstr '<exposureAgeResult>'];
outstr = [outstr '<t10St>0</t10St><delt10_int_St>0</delt10_int_St>0<delt10_ext_St>0</delt10_ext_St>'];
outstr = [outstr '<t26St>0</t26St><delt26_int_St>0</delt26_int_St>0<delt26_ext_St>0</delt26_ext_St>'];
outstr = [outstr '</exposureAgeResult>'];

% Put error string in diagnostics element
outstr = [outstr '<diagnostics>' errorString '</diagnostics>'];

outstr = [outstr '</calcs_v3_age_data>'];