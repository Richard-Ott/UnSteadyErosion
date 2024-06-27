function varargout = cosmowatersheds(SAMS,DEM,varargin)

% COSMOWATERSHEDS find upstream area for sample locations
%
% Syntax
%
%     OUT = cosmowatersheds(SAMS,DEM)
%     [OUT,MS] = cosmowatersheds(SAMS,DEM)
%     [OUT,MS] = cosmowatersheds(SAMS,DEM,'pn','pv',...)
%
% Description
%
%     COSMOWATERSHEDS finds contributing areas in the digital elevation 
%     model DEM for sample locations provided in the sample structure SAMS.
%     These are stored in a sample structure that is identical to SAMS, but
%     which contains additional fields. In standard work flow, the output
%     variable can be the same as the input variable SAMS and will have the
%     extra fields simply added to the structure.
%     Optional parameter name and value pairs are used to control the
%     automated search for the upstream drainage areas. If the drainage
%     area is known (e.g., when recalculating erosion rates for published
%     samples), the Area shall be stored under the field name "Area" in the
%     sample structure SAMS. It will then be used to search for stream
%     locations with similar upstream areas. The parameter 'areatol'
%     indicates the fraction deviation from this area that is tolerated
%     during the search. If the drainage area is unknown, the parameter
%     'minarea' controls the extent of the stream network.
%
% Input arguments
%
%     DEM       digital elevation model DEM (class GRIDobj)
%     SAMS      structure array with sample information
%
%     Parameter name/value pairs {default}
%
%     'minarea': scalar {1e6}
%     Minimum upstream area of channel network
%
%     'areatol': scalar {0.3}
%     Fractional area tolerance when expected upstream areas are provided
%
%     'maxsnap': scalar {3000}
%     Maximum snapping distance when relocating sample points to nearby
%     streams
%
%
% Output arguments
%
%     OUT      structure array that contains
%     .WSflag       flag indicating watershed delineation (1 = completed)
%     .snapdist     snapping distance
%     .OutletLat    Latitude of outlet
%     .OutletLon    Longitude of outlet
%     .OutletX      x coordinate of outlet
%     .OutletY      y coordinate of outlet
%     .WSArea       drainage area of watershed
%     .WSedgepixels number of watershed pixels that touch the edge of the
%                   input DEM (gives an assessment whether truncated or not)
%     .WSDEM        GRIDobj of the watershed area       
%     .WSDEMsource  name of the DEM source file
%     .WSLon        Longitude of watershed boundary
%     .WSLat        Latitude of watershed boundary
%
%     MS       geographic mapping structure of watershed boundaries
%
% Example
%
%     SAMS = cosmosampleread('dibiase_data.xlsx');
%     DEM0 = cosmogetdem(SAMS,'demtype','SRTMGL3','buffer',0.1);
%     DEM = reproject2utm(DEM0,90);
%     SAMS = cosmowatersheds(DEM,SAMS);
%     
%
% See also:
%     cosmosampleread, cosmogetdem
%
%
% Author: Dirk Scherler (scherler[at]gfz-potsdam.de)
% Date: 16. August, 2017



% Parse Inputs
p = inputParser;
p.FunctionName = 'cosmowatersheds';
p.addRequired('S', @(x) isstruct(x));
p.addRequired('DEM', @(x) isa(x,'GRIDobj'));

p.addParamValue('areafield','none',@(x) ischar(x));
p.addParamValue('areaunit','km',@(x) ismember(x,{'km','m'}));
p.addParamValue('areatol',0.3,@(x) isscalar(x));
p.addParamValue('maxsnap',3000,@(x) isscalar(x));
p.addParamValue('minarea',1e6,@(x) isscalar(x));
p.addParamValue('plot',false,@(x) islogical(x));
p.addParamValue('sinks',false,@(x) islogical(x));
p.addParamValue('sinkdepth',50,@(x) isnumeric(x));


%addParamValue(p,'area',true,@(x) isscalar(x));

p.parse(SAMS,DEM,varargin{:});
DEM   = p.Results.DEM;
SAMS = p.Results.S;
%areafield = p.Results.areafield;
areaunit = p.Results.areaunit;
areatol = p.Results.areatol;
maxsnap = p.Results.maxsnap;
minarea = p.Results.minarea;
doplot = p.Results.plot;
sinks = p.Results.sinks;
sinkd = p.Results.sinkdepth;

switch areaunit
    case 'km'
        areafac = 1e6;
    case 'm'
        areafac = 1;
end

nsamples = length(SAMS);
if ~isfield(SAMS,'Area')
    [SAMS(1:nsamples).Area] = deal(0);
end

% Process DEM
fprintf(1,'TT.COSMO: Finding contributing areas.\n');
tic
[x,y] = getcoordinates(DEM);

if (sinks)
    DEMfs = fillsinks(DEM);
    SINKS = DEMfs-DEM;
    S = SINKS>sinkd;
    FD = FLOWobj(DEM,'preprocess','carve','sinks',S);
else
    FD = FLOWobj(DEM,'preprocess','carve');
end
FA = flowacc(FD).*FD.cellsize.*FD.cellsize;
STR0 = FA>minarea;
STR0 = STREAMobj(FD,STR0);


% Initialize output variable
OUT = SAMS;
if ~isfield(OUT,'WSflag')
    [OUT.WSflag] = deal(true);
end
flag = [OUT.WSflag];

% Find outlets that are contained within the current DEM
if isfield(SAMS,'Lat')
    LAT = [SAMS.Lat];
    LON = [SAMS.Lon];
elseif isfield(SAMS,'X')
    LAT = [SAMS.Y];
    LON = [SAMS.X];
else
    error('TT.COSMO: Cannot find coordinates in sample structure!');
end
[out_x,out_y] = mfwdtran(DEM.georef.mstruct,LAT,LON);
ix = out_x < max(x) & out_x > min(x) & out_y < max(y) & out_y > min(y) & flag;
inside_x = out_x(ix);
inside_y = out_y(ix);
inside_ix = find(ix);

if (doplot)
    figure(1)
    imagesc(DEM)
end

% Loop over outlets within the DEM
for k = 1 : length(inside_ix)
    
    this_x = inside_x(k);
    this_y = inside_y(k);
    this_ix = inside_ix(k);
    this_area = SAMS(this_ix).Area.*areafac;
    
    % Snap to stream
    if this_area>0
        STR = FA>(1-areatol)*this_area & FA<(1+areatol)*this_area;
        STR = STREAMobj(FD,STR);
    else
        STR = STR0;
    end
    [xn,yn,ix,snapdist,~] = snap2stream(STR,this_x,this_y,'maxdist',maxsnap);
    
    if (doplot)
        hold on
        plot(STR,'k-')
        plot(this_x,this_y,'ro')
        hold off
    end
    
    if ~isnan(ix)
        
        % If we reached this far, it seems we find a stream with an
        % upstream area that is either smaller than the minimum area or 
        % similar to what has been published and not too far from where 
        % the outlet is supposed to be located.
        % Now, let's get the catchment
        
        OUT(this_ix).snapdist = snapdist;
        
        % Relocated outlet
        [latn,lonn] = minvtran(DEM.georef.mstruct,xn,yn);
        OUT(this_ix).OutletLat = latn;
        OUT(this_ix).OutletLon = lonn;
        OUT(this_ix).OutletX = xn;
        OUT(this_ix).OutletY = yn;
        
        % Drainage area
        DB = dependencemap(FD,ix);
        OUT(this_ix).WSArea = sum(sum(DB.Z)).*DB.cellsize.*DB.cellsize;
        
        % Edge-touching pixels
        BW = DB;
        BW.Z(2:end-1,2:end-1) = false;
        OUT(this_ix).WSedgepixels = sum(BW.Z(:));
        
        
        % DEM
        try
            DEMcrop = crop(DEM,DB,nan);
        catch ME
            warning('TT.COSMO: DEM extraction failed during cropping. Continuing with empty DEM.')
        	DEMcrop = [];
        end
        OUT(this_ix).WSDEM = DEMcrop;
        OUT(this_ix).WSDEMsource = DEM.name;
        
        if not(isempty(DEMcrop))
            OUT(this_ix).WSflag = false;
        end
        
        % Catchment boundary
        DB.Z = double(DB.Z);
        DB.Z(DB.Z==0) = nan;
        M = GRIDobj2polygon(DB);
        [WSLat,WSLon] = minvtran(DB.georef.mstruct,M.X,M.Y);
        OUT(this_ix).WSLon = WSLon;
        OUT(this_ix).WSLat = WSLat;
        
        if this_area>0
            fprintf(1,'TT.COSMO: Catchment ID %s done. Provided area: %1.2f km^2, extracted area: %1.2f km^2.\n',...
                OUT(this_ix).ID,this_area./1e6,OUT(this_ix).WSArea./1e6);
        else
            fprintf(1,'TT.COSMO: Catchment ID %s done. Extracted area: %1.2f km^2.\n',...
                OUT(this_ix).ID,OUT(this_ix).WSArea./1e6);
        end
        
        if (doplot)
            hold on
            plot(xn,yn,'go')
            plot(M.X,M.Y,'w-')
            hold off
        end
        
    else
        fprintf(1,'TT.COSMO: Catchment ID %s watershed extraction failed. Snapping distance likely too long.\n',...
            OUT(this_ix).ID);
        
        OUT(this_ix).WSflag = false;
        OUT(this_ix).snapdist = NaN;
        OUT(this_ix).OutletLat = NaN;
        OUT(this_ix).OutletLon = NaN;
        OUT(this_ix).OutletX = NaN;
        OUT(this_ix).OutletY = NaN;
        OUT(this_ix).WSArea = NaN;
        OUT(this_ix).WSedgepixels = NaN;
        OUT(this_ix).WSDEM = NaN;
        OUT(this_ix).WSDEMsource = NaN;
        OUT(this_ix).WSLon = NaN;
        OUT(this_ix).WSLat = NaN;
        
    end % snapping condition
    
end % loop over outlets within the present DEM
toc

% Assemble output
varargout{1} = OUT;
if nargout>1
    % Create a mapping structure that can be exported as a shapefile
    MS = struct;
    for k = 1:length(OUT)
        MS(k).Geometry = 'Polygon';
        MS(k).ID = OUT(k).ID;
        MS(k).X = OUT(k).WSLon;
        MS(k).Y = OUT(k).WSLat;
        MS(k).N10 = OUT(k).N10;
        MS(k).N10sigma = OUT(k).N10sigma;
        MS(k).Area = OUT(k).WSArea;
    end 
    varargout{2} = MS;
end

