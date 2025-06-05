function walker = Egholm_MCMC(nWalks,Nobs,dNobs,mini,prior_range,forward_model,Nlogical)
Nobs=Nobs(Nlogical);
dNobs=dNobs(Nlogical);

%save some general MCMC paramers
burnin = 1e4;%5e4;              % burn in iterations (accepted)
Nmod = 1e5;%2e5;                % number of accpected samples
Nmax = 2e5;%2e6;                % maximum number of models
Nmp  = size(mini,1);       % number of model parameters
Temp = 1;
Cobs = diag(dNobs.^2);
Cobsinv = inv(Cobs);

% % real samples need a higher variance tolerance to get going, so I start
% % with 5 times the original error before switching after the MCMC got a
% % foothold. This is only needed for Crete samples for some reason
% CobsSTART = diag((dNobs*3).^2);
% CobsinvSTART = inv(CobsSTART);


du0 = prior_range(:,2) - prior_range(:,1);
umin = prior_range(:,1);
umax = prior_range(:,2);


%loop walkers
wb = waitbar(0,'Opening the box of Pandora...');
for nw = 1:nWalks

    %walker starting point - for initial parameter vector
    u = mini(:,nw);
                
 
    % initialize
    minres = 1e20;       % ???????????
    res_current = 1e20;  % residual inital guess
    restot = 0;          % residual total (observational error)
    acount = 0;          % accepted samples after burn-in
    bcount = 0;          % accepted samples during burn-in
    rcount = 0;          % rejected samples
    accrat = 0;          % acceptance ratio
    status = zeros(Nmax,1);
    erosion_rec = zeros(Nmax,1);
    up_rec = zeros(Nmp,Nmax);
    u_rec = zeros(Nmp,Nmax);
    N10_rec = zeros(Nmax,1);
    N14_rec = zeros(Nmax,1);
    restot_rec = zeros(Nmax,1);
    accrat_rec = zeros(Nmax,1);
    k_rec = zeros(Nmax,1);
    duR = zeros(Nmp,1);

    accfac = 1e-2;    % ???????????? accceptance factor?
    ktarget = 0.025;  % ???????????? target step size? seems not used
    
    %run models
    mi = 0;           % counter variable for model proposals
    k = 0.01;         % stepsize, 0.01, 0.05

    while ((mi < Nmax)&(acount < Nmod))
        
        mi = mi + 1;
    
       %***** step length ******
       
       %acceptance ratio
        if (mi > 100)
            accrat = (sum(abs(status((mi-100):(mi-1))))+1)/100;
        elseif (bcount < burnin)
            accrat = 0.1;
        else 
            accrat = 0.3;
        end
        
        %burnin
        if (bcount < 0.5*burnin) 
            
            k = k*((1-accfac) + accfac*accrat/0.1);     
            
        elseif (bcount < burnin) 

            k = k*((1-accfac) + accfac*accrat/0.2);
            
        elseif (acount < Nmod)

            k = k*((1-accfac) + accfac*accrat/0.3);    
                        
        end        
        
        if (bcount < burnin) Temp = 1.0 + 20.0*(burnin-bcount)/burnin;
        else Temp = 1.0;
        end
            
        % weird RO edit ------------------------------------------
        % if k < 0.005  % make sure step sizes dont get too small, this seems to happen
        %     k = 0.005;
        % end
        %********* propose new parameters ************
        
        %random step
        du = 0.5*randn(Nmp,1).*du0(:);
        
        %proposed model
        up = u(:) + k*du(:);
     
        
        %retake
        while (any(up(:) < umin(:))|any(up(:) > umax(:)))
        
            %random step
            du = 0.5*randn(Nmp,1).*du0(:);

            %proposed model
            up = u(:) + k*du(:);
            
            % THIS IS AN RO EDIT, NA DNOT A GOOD ONE
            if any(up<0)  % ideally I shouldnt interfere with random sampling, but correcting values if the algorithm wants to propose negative parameters speeds up the search by orders of magnitude
                up(up<0)= umin(up<0);
            end
            
        end
        
        
        %********** Forward model *****************
        gmp = forward_model(up);
               
        %Acceptance critieria
        % if (acount + bcount) > 1000
        restot = (Nobs(:)-gmp(:))'*Cobsinv*(Nobs(:)-gmp(:))/Temp;
        % else
            % restot = (Nobs(:)-gmp(:))'*CobsinvSTART*(Nobs(:)-gmp(:))/Temp;
        % end
        rfrac = exp(-0.5*restot)/exp(-0.5*res_current);
        alpha = rand(1);
    
    
        %if model is accepted
        if ((alpha < rfrac)|(mi == 1))
        
            u = up;
            gm = gmp;
            res_current = restot;
        
            %accepted after burnin
            if (bcount > burnin)
        
                status(mi) = 1;
                acount = acount + 1;
                
            %if accepted during burnin    
            else
                status(mi) = -1;
                bcount = bcount + 1;          
            end
        
        %rejected
        else
            status(mi) = 0;
            rcount = rcount + 1;                     
        end
       
        %save things
        up_rec(:,mi) = up(:);
        u_rec(:,mi) = u(:);
        gm_rec(:,mi) = gm(:);
        restot_rec(mi) = res_current;
        accrat_rec(mi) = accrat;
        k_rec(mi) = k;              
    end

    
    walker{nw}.status = status(1:mi);
    walker{nw}.up = up_rec(:,1:mi);
    walker{nw}.u = u_rec(:,1:mi);
    walker{nw}.gm = gm_rec(:,1:mi);
    walker{nw}.restot = restot_rec(1:mi);
    walker{nw}.acount = acount;
    walker{nw}.bcount = bcount; 
    walker{nw}.rcount = rcount;
    walker{nw}.accrate = accrat_rec(1:mi);
    walker{nw}.kstep = k_rec(1:mi);
    
    waitbar(nw/nWalks,wb)
    disp([num2str(acount), ' accpected samples'])  % tell user how many samples were accepted
    disp(['acceptance ratio = ' ,num2str(acount/mi)])  % tell user how many samples were accepted
end
close(wb)

end