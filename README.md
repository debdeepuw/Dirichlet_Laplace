# Dirichlet_Laplace
Dirichlet Laplace prior for a Normal Means Model
Objective: estimation of a high-dim sparse normal mean theta_0  %%
Model: $y ~ N(\theta, I_n)$
prior for theta: theta_j ~ N(0, psi_j phi_j^2 tau^2), psi_j i.i.d. exp(1/2), 
phi ~ Dir(a,..., a), tau ~ gamma(n*a, 1/2)
