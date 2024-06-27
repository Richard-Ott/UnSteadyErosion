function out = get_LmSF(in,datadir)

% Interpolates into precalculated scaling factor grids to obtain Lm
% scaling factors for specified Rc and pressure vectors. 
%
% This code ingests an input structure with the following information:
%   cell array sfdata.t contains time vectors
%   cell array sfdata.dipRc contains Rc at corresponding times
%   vector pressure is atmospheric pressure
%
% Adds a cell array containing corresponding vectors of scaling factor
% and returns the structure. 

% Load grid file

load([datadir 'LSDn2015_grid.mat']);

% Load sf file

load([datadir 'Lm2015_all.mat']);
 

% Now do interpolation.    

for a = 1:length(in.tmax);
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
      
    % Do interpolation.
    
    
    % assign vector to cell in field in.Lm
    in.Lm{a} = interp2(sfgrid.Rc,sfgrid.p,sf,in.Rc{a},temp_p);  
    
end;

% Return

out = in;
    
    
    
    
    
    
    
    