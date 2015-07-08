function P = S_Exp_Para(mode)

switch mode
  case 'paper-2AFC-corr'
    P.number_orientations=2;
    P.prior_task=[0 1]; % [cardinal, oblique]
    P.number_locations=1;
    P.dimension_X=1024;
    P.dimension_G=256;
    P.kappa_O=[1 0]; % attended and unattended
    P.kappa_G=3;
    P.delta=5; % ??
    P.stimulus_regime='static';
    P.stimulus_contrast=zeros(1,P.number_orientations);
    P.number_repetitions=5; % number of trials
    P.number_burn_in=0; % number of burn-in samPles
    P.number_samples_to_use=100; %Number  of non-burn samples to be used for evidence
    P.number_samples_per_evidence=5;
    P.phi_O=[[0:2:2*(P.number_orientations-1)];[1:2:2*P.number_orientations-1]]*pi/2/P.number_orientations;
    P.pO=ones(1,P.number_orientations)/P.number_orientations;
    P.pL=ones(1,P.number_locations)/P.number_locations;
    P.tauStyle=1;
    P.sigmaStyle=1e-10;
    P.odds_inc=1/P.number_samples_per_evidence;
    P.nT = 2; % Number of possible tasks

  case 'paper-2AFC-PK'
    % to be added
    
    
  otherwise
    warning('invalid option');
end