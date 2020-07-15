function [IRGN] = setup_IRGN()

IRGN.lambda=1e-5;
IRGN.lambda_min=6e-1;
IRGN.gamma=1e-1;
IRGN.qgamma=2;
IRGN.qlambda=0.7;
IRGN.gamma_max=1;

IRGN.GNsteps=20;
IRGN.PDsteps=100;
IRGN.tau_init=1/sqrt(0.5*(18.0 + sqrt(33)));
IRGN.beta=400;
IRGN.theta_init=1;
IRGN.delta=1;
IRGN.mu=0.5;
IRGN.tol=5e-3;
IRGN.stag=1;

IRGN.alfa0=2;
IRGN.alfa1=1;

end

