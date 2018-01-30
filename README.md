# Dirichlet_Laplace prior for optimal shrinkage

## Objective: estimation of a high-dim sparse normal mean theta0  

Model: y ~ N(theta, I_n)

prior for theta: theta_j ~ N(0, psi_j phi_j^2 tau^2), psi_j i.i.d. exp(1/2), 

phi ~ Dir(a,..., a), tau ~ gamma(n*a, 1/2)

## Use the runDL.m script to simulate data and run DL as follows

```matlab

%% simulate data  %%
n=500;                                   % dimension
A=5;                                     % signal strength
qn=10;                                   % # of non-zero signals
theta0=[A*ones(qn,1);  zeros(n- qn, 1)]; % create theta0
y=normrnd(theta0, 1);                    % Draw y ~ N(theta0, I_n)


%% set hyperparameter a of Dirichlet Laplace prior %%
a=1/40;

%% MCMC variables
nrun = 10000; burn = 5000; thin = 1; 

%% Call the function DL %%
plt=1;                                   % Plot data points, posterior mean and credible intervals
save_samples=0;                          % Do not save samples in a separate text file

[pmeantht,pmedtht,thtout]=DL(y,a,burn,nrun,thin,plt,save_samples); 

```
## Illustrative plot 

