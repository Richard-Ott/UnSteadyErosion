function out = intNmu(thise,thisp,pp,logee,Nmu)

% Objective function for muon inventory interpolation

if isa(thise, 'single')
    thise=double(thise); 
end % this is necessary because I store inversion results as single to save memory

% Try once, should work
out = interp2(pp,logee,Nmu,thisp,log10(thise));

if isnan(out);
    % Presumably exceeds ee bounds
    if log10(thise) > max(max(logee)); 
        thise = 10.^(max(max(logee)));
    elseif log10(thise) < min(min(logee)); 
        thise = 10.^(min(min(logee)));
    end;
    out = interp2(pp,logee,Nmu,thisp,log10(thise));
end;

if isnan(out);
    % Must exceed bounds in p
    if thisp > max(max(pp));
        thisp = max(max(pp));
    elseif thisp < min(min(pp));
        thisp = min(min(pp));
    end;
    out = interp2(pp,logee,Nmu,thisp,log10(thise));
end;

if isnan(out);
    error('intNmu.m: interp2 gives NaN after bounds checks');
end;


