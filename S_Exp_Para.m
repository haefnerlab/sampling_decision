function P = S_Exp_Para()

P.number_orientations=2;
P.prior_task=[0 1];
P.number_locations=1;
P.dimension_X=16;
P.dimension_G=4;
P.kappa_O=[1 0]; % attended and unattended
P.kappa_G=3;
P.tau=5; % ??
P.stimulus_regime='static';  
P.stimulus_contrast=zeros(1,P.number_orientations);
P.number_repetitions=100; % number of trials
P.number_burn_in=20; % number of burn-in samPles
P.number_samples_to_use=80; %Number  of non-burn samples to be used for evidence
P.number_samples_per_evidence=5;
P.phi_O=[[0:2:2*(P.number_orientations-1)];[1:2:2*P.number_orientations-1]]*pi/2/P.number_orientations;
P.pO=ones(1,P.number_orientations)/P.number_orientations;
P.pL=ones(1,P.number_locations)/P.number_locations;
P.tauStyle=1;
P.sigmaStyle=1;
P.odds_inc=1/P.number_samples_per_evidence;
P.nT = 2; % Number of possible tasks
end