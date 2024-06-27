function SAMS = cosmosampleread(varargin)

% COSMOSAMPLEREAD load input data
%
% Syntax
%
%     SAMS = cosmosampleread
%     SAMS = cosmosampleread(fname)
%
% Description
%
%     COSMOSAMPLEREAD is used to load cosmogenic nuclide sample data from 
%     an Excel spreadsheet. The data needs to be stored in the first spread
%     sheet and have the following columns:
%       ID          - Sample ID  
%       Lat         - Latitude (decimal degrees) of sampling position
%       Lon         - Longitude (decimal degrees) of sampling position
%       N10         - 10Be concentration (atoms/g)
%       N10sigma    - 10Be concentration uncertainty (atoms/g)
%       N10standard - 10Be AMS standard (see CRONUS documentation)
%       Age         - Age of sample (years; use 0 for modern river samples)
%       Density     - Density of bedrock in drainage area (g/cm3)
%
% Input arguments
%
%     fname     name of the Excel-file (xls/xlsx)
%
% Output arguments
%
%     SAMS      Sample structure array that contains
%     .ID       Sample ID (character array)
%     .Lat      Latitude of sampling point (decimal degrees)
%     .Lon      Longitude of sampling point (decimal degrees)
%     .N10      10Be concentration (atoms/g)
%     .N10sig   10Be concentration uncertainty (atoms/g)
%     .N10std   10Be standard (character array)
%     .Age      Sample age (years)
%     .Density  Bedrock density (g/cm3)
%
% Example
% 
%   SAMS = cosmosampleread('dibiase_data.xlsx');
%
%
% See also: xlsread
%
%
% Author: Dirk Scherler (scherler[at]gfz-potsdam.de)
% Date: 16. August, 2017


% Check whether GUI needed for selecting file
if nargin == 0
    [fn,pn] = uigetfile({'*.xlsx';'*.xls'});
    source = [pn,fn];
else
    source = varargin{1};
end


% Read sample data from xls/xlsx-file
[num,txt,raw] = xlsread(source,1);
% T = readtable(source); EDIT: replace function
fields = txt(1,:); % the first row needs to contain the column names
xx = strcmp(fields,'');
fields = fields(~xx);
fieldtag = true(size(fields));
raw = raw(2:end,:);
SAMS = struct;
nsamples = size(num,1);
ix = strcmp(raw,'nan');
raw(ix) = {NaN};


% Set required fields
requiredfields = {'ID','Lat','Lon','N10','N10sigma','N10standard','Age'};
for i = 1 : length(requiredfields)
    istr = strcmpi(fields,requiredfields{i});
    if sum(istr)<1
        error('TT.cosmo: field ''%s'' not found in file ''%s''.\n',requiredfields{i},source);
    elseif sum(istr)>1
        error('TT.cosmo: field ''%s'' appears more than once in file ''%s''.\n',requiredfields{i},source);
    end
    cv = raw(:,istr);
    [SAMS(1:nsamples).(requiredfields{i})] = deal(cv{:});
    fieldtag(istr) = false;
end

% Add optional fields
optionalfields = fields(fieldtag);
for i = 1 : length(optionalfields)
    istr = strcmpi(fields,optionalfields{i});
    cv = raw(:,istr);
    [SAMS(1:nsamples).(optionalfields{i})] = deal(cv{:});
    fieldtag(istr) = 0;
end

% Check standards
allstds = {'07KNSTD','KNSTD','NIST_Certified','NIST_30000','NIST_30200',...
    'NIST_30300','NIST_30600','NIST_27900','BEST433','S555','S2007',...
    'BEST433N','S555N','S2007N','LLNL31000','LLNL10000','LLNL3000',...
    'LLNL1000','LLNL300'};
for i = 1 : length(SAMS)
    if sum(strcmpi(SAMS(i).N10standard,allstds))<1
        error('TT.cosmo: sample ''%s'' uses unknown standard ''%s''.\n',SAMS(i).ID,SAMS(i).N10standard);
    end
end
 
    
end



