function H = conc_modelledVSobserved(models,N10,dN10,N14,dN14)
% plots the compairson of modelled vs observed N10 and N14 concentrations
% as output from the gwmcmc ensemlbe sampler against the measured
% concentrations

hold on; box on; grid on;
xlabel(['N10 (atoms/g)']);
ylabel(['N14 (atoms/g)']);

    
for ns=1:length(N10)
    
    
    Be = [];
    C = [];
    for nw = 1:model.Nwalk  % loop through all random walks and concatenate all modelled concentrations
        I = find((model.walker{nw}.status == 1)&(model.walker{nw}.restot <= maxres));
%          I = find((model.walker{nw}.status > -1)&(model.walker{nw}.restot <= maxres));  % get indices of samples
         % get modelled conentrations
        Be = [Be(:);model.walker{nw}.gm(n0+1,I)'];
        C = [C(:);model.walker{nw}.gm(n0+2,I)'];
    end
    
    Beint = linspace(min(Be)-1e4,max(Be+1e4),40)';
    Cint = linspace(min(C)-1e4,max(C)+1e4,40)';
    
    % make contour plot of modelled data
    N = hist3([C,Be],{Cint Beint});
    N = N/sum(N(:));
    [X,Y] = meshgrid(Beint,Cint);
    contour(X,Y,N,40);
    
    % plot original data value
    errorbar(model.data{ns}.N10,model.data{ns}.N14,model.data{ns}.dN10,'horizontal','.k');
    errorbar(model.data{ns}.N10,model.data{ns}.N14,model.data{ns}.dN14,'vertical','.k');
    
    %subplot(np(1),np(2),n0+2);
    %hold on; box on; grid on;
  
    
end
legend({'modelled concentrations','measured concentrations'})
title('modelled vs measured nuclide concentrations')


end