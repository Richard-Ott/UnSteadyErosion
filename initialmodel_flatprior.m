function mini  = initialmodel_flatprior(prior_ranges,nWalks)
% This function determines an initial model based on the provided prior
% ranges. Initial models are drawn randomly from within the prior range.
% Richard Ott, 2024

    [nparas,~] =  size(prior_ranges);
    
    mini = nan(nparas, nWalks);
    for i = 1:nWalks
        mini(:,i) = prior_ranges(:,1) + rand(nparas,1) .* (prior_ranges(:,2) - prior_ranges(:,1));
    end

end