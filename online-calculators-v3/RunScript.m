% script to Cronus v3

load consts_v3.mat

in = validate_v3_input('PH-1	30	0	500	std	0	2.65	1	0.0000	1999;PH-1	Be-10	quartz	11130	3717	KNSTD;PH-1	C-14	quartz	43268	1238');

control.resultType =  'short';

output = get_erates_v3(in,control,consts);