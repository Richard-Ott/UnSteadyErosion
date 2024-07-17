function H = conc_modelledVSobserved(m,N10,dN10,N14,dN14)
% plots the compairson of modelled vs observed N10 and N14 concentrations
% as output from the gwmcmc ensemlbe sampler against the measured
% concentrations

% currently only for N10 and N14 combination. Should be updated to include
% 26Al

hold on; box on; grid on;

Nmodels= [];
for i = 1:length(m)
    Nmodels = [ Nmodels, m{i}.gm(:,m{i}.status == 1)];
end

nSamp = size(Nmodels,1)/2;

H= figure;
set(gcf,'units','normalized','position',[.2,.25,.4,.6]);

hold on; box on; grid on;
xlabel('N10 (atoms/g)');
ylabel('N14 (atoms/g)');

%[np,nn] = numSubplots(model.Nsnr*model.Nds);
    
for i = 1:nSamp
      
    Be = Nmodels(i,:);
    C =  Nmodels(i+nSamp,:);
    
    Beint = linspace(min(Be),max(Be),40)';
    Cint = linspace(min(C),max(C),40)';
    
    N = hist3([C',Be'],{Cint Beint});
    N = N/sum(N(:));
    [X,Y] = meshgrid(Beint,Cint);
    contour(X,Y,N,40);
    
    errorbar(N10(i),N14(i),dN14(i),dN14(i),dN10(i),dN10(i),'Color','k');
end

end