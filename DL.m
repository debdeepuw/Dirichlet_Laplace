function[pmeantht,pmedtht,thtout]=DL(y,a,burn,nrun,thin,plt,save_samples)
%% Objective: estimation of a high-dim sparse normal mean theta_0  %%
%% Model: y ~ N(theta, I_n) %%
%% prior for theta: theta_j ~ N(0, psi_j phi_j^2 tau^2), psi_j i.i.d. exp(1/2), 
%% phi ~ Dir(a,..., a), tau ~ gamma(n*a, 1/2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% This function employs the algorithm proposed in Section 2.4 "Dirichlet Laplace priors for optimal shrinkage" by
%% Bhattacharya et. al. (2015) published in the Journal of the American Statistical Association, 110 (512): 1479-1489. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Input: y=response, a n*1 vector %%
%%        a= parameter for the Dirichlet
%%        burn= number of burnin MCMC samples %%
%%        nrun= total number of posterior draws, only nrun-burn of them will be saved %%
%%        thin= thinning parameter of the chain %%
%%        plt= binary indicator for plotting posterior median with credible intervals
%%        (blue: datapoints, red: posterior median of theta, black lines: 95% credible intervals
%%        save_samples= binary indicator whether posterior samples should be saved in a file or not %%


%% Output: 
%%         pmeantht= posterior mean of theta, a n by 1 vector%%
%%         pmedtht=posterior median of theta, a n by 1 vector %%
%%         thtout=posterior samples of theta - a matrix of size nrun-burn by n%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


addpath('aux_files');

% -- define MCMC stuff -- %
effsamp = (nrun-burn)/thin;

% -- define parameters -- %
n=length(y); psi = ones(n,1); T = ones(n,1); phi = T/n; S = ones(n,1); tau = 1; 


% -- define output files -- %
thtout = zeros(effsamp,n);
phiout = zeros(effsamp,n);


% -- start gibbs sampling -- %

for i = 1:nrun
    
    % -- update theta -- %
    temp1 = tau^2*(psi.*(phi.^2));
    sigtht = 1./(1 + 1./temp1); mutht = sigtht.*y;
    theta = normrnd(mutht,sqrt(sigtht));
    
    % -- update psi -- %
    mu = tau*(phi./abs(theta));     
    V = normrnd(0,1,[n,1]).^2; U = unifrnd(0,1,[n,1]);
    temp2 = mu.*V; temp3 = sqrt(4*temp2 + temp2.^2);
    W = mu + 0.5*(V.*(mu.^2)) - 0.5*(mu.*temp3);
    locs = (U <= mu./(mu + W));
    psitil = locs.*W + (1-locs).*((mu.^2)./W); psi = 1./psitil;
    
    % -- update tau -- %
    tau = randraw('gig',[n*a-n, 2*sum(abs(theta)./phi), 1],[1, 1]); % Note: order of rho & chi in randraw is opposite
    
    % -- update phi -- %
    uT = unifrnd(0,exp(-1./(2*S))); lb = 1./(2*log(1./uT));
    Flb = gamcdf(lb,1-a,1./abs(theta)); uu = unifrnd(Flb,1); uu = min(uu,1-(1e-20));
    S = gaminv(uu,1-a,1./abs(theta));   
    T = 1./S; phi = T/sum(T); phi(phi <= (1e-20)) = 1e-20;
    
    
        
    if mod(i,100) == 0
        disp(i)
    end
    
    if i > burn && mod(i, thin)== 0
        thtout((i-burn)/thin,:) = theta';
        phiout((i-burn)/thin,:) = phi';
        
    end
end

%% Compute posterior mean, median  and credible intervals of theta %%

pmeantht = mean(thtout)'; pmedtht = median(thtout)'; 
puqtht = quantile(thtout,0.975)'; plqtht = quantile(thtout,0.025)';

%% Plot posterior median with data points and credible intervals %%

if plt
    
    figure(1);  ff = figure(1); set(ff,'Position',[2 300 900 400]); set(gcf,'PaperPositionMode','auto');
    grid = (1:n)'; plot([1,n],[0,0],'k','LineWidth',2.5); hold on;
    plot(grid,pmedtht,'r.','MarkerSize',15); hold on;
    plot(grid,y,'b.','MarkerSize',14); hold on;
    for j=1:n
        fill([j,j,j],[plqtht(j),puqtht(j),plqtht(j)],[0.75 0.75 0.75]); hold on;
    end
    set(gca,'FontSize',14,'FontWeight','bold','tickdir','out');
end

%% save samples %%
if save_samples
       dlmwrite('post_tht_samples.txt',thtout);
end
