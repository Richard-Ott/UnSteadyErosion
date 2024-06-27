function out = get_LSDnSF(in,datadir)

% Interpolates into precalculated scaling factor grids to obtain LSDn
% scaling factors for specified Rc, pressure, and solar parameter vectors. 
%
% This code ingests an input structure with the following information:
%   cell array sfdata.t contains time vectors
%   cell array sfdata.dipRc contains Rc at corresponding times
%   cell array sfdata.S contains solar parameter
%   cell array nuclide is nuclide code
%   vector pressure is atmospheric pressure
%
% Adds a cell array containing corresponding vectors of scaling factor
% and returns the structure. 

% Load grid file

load ([datadir 'LSDn2015_grid']);

% First decide which data files to load

if any(strcmp('N3quartz',in.nuclide));
    % Load He-3-in-quartz scaling block
    load([datadir 'LSDn2015_N3quartz']);
end;

if any(strcmp('N3pyroxene',in.nuclide)) || any(strcmp('N3olivine',in.nuclide));
    % Load He-3-in-ol/px scaling block
    load([datadir 'LSDn2015_N3pxol']);
end;

if any(strcmp('N10quartz',in.nuclide));
    % Load Be-10-in-quartz scaling block
    load([datadir 'LSDn2015_N10quartz']);
end;

if any(strcmp('N14quartz',in.nuclide));
    % Load C-14-in-quartz scaling block
    load([datadir 'LSDn2015_N14quartz']);
end;

if any(strcmp('N21quartz',in.nuclide));
    % Load Ne-21-in-quartz scaling block
    load([datadir 'LSDn2015_N21quartz']);
end;

if any(strcmp('N26quartz',in.nuclide));
    % Load Al-26-in-quartz scaling block
    load([datadir 'LSDn2015_N26quartz']);
end;

if any(strcmp('N10pyroxene',in.nuclide));
    % Load Be-10-in-quartz scaling block; this is same for quartz and px
    if ~exist('sfN10quartz','var');
        load([datadir 'LSDn2015_N10quartz']);
    end
    % Create an identical structure for px so that the below will run
    sfN10pyroxene = sfN10quartz;
end;

% Now do interpolation.    

for a = 1:length(in.tmax);
    % Choose appropriate nuclide
    eval(['d = sf' in.nuclide{a} ';']);
    
    % Make pressure vector. This is planning ahead for time-dependent
    % pressure changes. 
    if iscell(in.pressure); 
        % Case in.pressure is cell array with time-dependent p
        temp_p = in.pressure{a};
    else
        % Case in.pressure is vector with non-time-dependent p
        temp_p = zeros(size(in.tmax{a})) + in.pressure(a);
    end;
    
     % Catch excessive Rc.
    bigRc = find(in.Rc{a} > 20);
    if ~isempty(bigRc);
        in.Rc{a}(bigRc) = 20 + zeros(size(bigRc));
    end;
    
    % Do interpolation...with linear approximation for solar variability
    % assign vector to cell in field in.LSDn
    in.LSDn{a} = interp2(sfgrid.Rc,sfgrid.p,d.b,in.Rc{a},temp_p) + in.S{a}.*interp2(sfgrid.Rc,sfgrid.p,d.m,in.Rc{a},temp_p);
end;

% Return

out = in;
    
    
    
    
    
    
    
    