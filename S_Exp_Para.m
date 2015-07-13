function P = S_Exp_Para(mode)

switch mode
  case 'paper-2AFC-corr'
    P.number_orientations=2;
    P.prior_task=[0 1]; % [cardinal, oblique]
    P.nL=1;
    P.nX=1024;
    P.nG=256;
    P.kappa_O=[1 0]; % attended and unattended
    P.kappa_G=3;
    P.delta=15; % strength of X-G coupling for corr & CPs
    P.stimulus_regime='static';
    P.stimulus_contrast=zeros(1,P.number_orientations);
    P.number_repetitions=16; % number of trials
    P.n_zero_signal=20; % number of frames before onset of stimulus
    P.number_burn_in=0; % number of burn-in samPles
    P.number_samples_to_use=100; %Number  of non-burn samples to be used for evidence
    P.number_samples_per_evidence=5; % for dynamic-switching-signal-blocked
    P.phi_O=[[0:2:2*(P.number_orientations-1)];[1:2:2*P.number_orientations-1]]*pi/2/P.number_orientations;
    P.pO=ones(1,P.number_orientations)/P.number_orientations;
    P.pL=ones(1,P.number_locations)/P.number_locations;
    P.tauStyle=0.1; % initial value for style
    P.sigmaStyle=0.1;
    P.odds_inc=1/20; % with 80 samples/sec, every 250ms independent signal
    P.nT = 2; % Number of possible tasks

  case 'paper-2AFC-PK'
    P.number_orientations=2;
    P.prior_task=[0 1]; % [cardinal, oblique]
    P.nL=1;
    P.nX=1024;
    P.nG=256;
    P.kappa_O=[1 0]; % attended and unattended
    P.kappa_G=3;
    P.delta=80; % strength of X-G coupling
    P.stimulus_regime='dynamic-switching-signal-blocked';
    P.stimulus_contrast=zeros(1,P.number_orientations);
    P.number_repetitions=16; % number of trials; 
    P.n_zero_signal=50; % number of frames before onset of stimulus
    P.number_burn_in=0; % number of burn-in samPles
    P.number_samples_to_use=100; %Number  of non-burn samples to be used for evidence
    P.number_samples_per_evidence=2; % for dynamic-switching-signal-blocked
    P.phi_O=[[0:2:2*(P.number_orientations-1)];[1:2:2*P.number_orientations-1]]*pi/2/P.number_orientations;
    P.pO=ones(1,P.number_orientations)/P.number_orientations;
    P.pL=ones(1,P.number_locations)/P.number_locations;
    P.tauStyle=1;
    P.sigmaStyle=1e-10;
    P.odds_inc=1/20; % with 80 samples/sec, every 250ms independent signal
    P.nT = 2; % Number of possible tasks
    
  case 'paper-corr-performance'
    P.number_orientations=2;
    P.prior_task=[1 0]; % [cardinal, oblique]
    P.number_locations=1;
    P.dimension_X=512;
    P.dimension_G=64;
    P.kappa_O=[1 0]; % attended and unattended
    P.kappa_G=3;
    P.delta=0.08; % strength of X-G coupling for corr & CPs, 80/1024
    P.stimulus_regime='static';
    P.stimulus_contrast=zeros(1,P.number_orientations);
    P.number_repetitions=512; % number of trials
    P.n_zero_signal=20; % number of frames before onset of stimulus
    P.number_burn_in=0; % number of burn-in samPles
    P.number_samples_to_use=100; %Number  of non-burn samples to be used for evidence
    P.number_samples_per_evidence=5; % for dynamic-switching-signal-blocked
    P.phi_O=[[0:2:2*(P.number_orientations-1)];[1:2:2*P.number_orientations-1]]*pi/2/P.number_orientations;
    P.pO=ones(1,P.number_orientations)/P.number_orientations;
    P.pL=ones(1,P.number_locations)/P.number_locations;
    P.tauStyle=0.1; % initial value for style
    P.sigmaStyle=0.1;
    P.odds_inc=1/20; % with 80 samples/sec, every 250ms independent signal
    P.nT = 2; % Number of possible tasks

  case 'test-contrast'
    P.number_orientations=2;
    P.prior_task=[1 0]; % [cardinal, oblique]
    P.number_locations=1;
    P.dimension_X=256;
    P.dimension_G=64;
    P.kappa_O=[1 0]; % attended and unattended
    P.kappa_G=3;
    P.delta=80; % strength of X-G coupling in paper
    P.stimulus_regime='static';
    P.stimulus_contrast=zeros(1,P.number_orientations);
    P.number_repetitions=128; % number of trials
    P.n_zero_signal=20; % number of frames before onset of stimulus
    P.number_burn_in=0; % number of burn-in samPles
    P.number_samples_to_use=100; %Number  of non-burn samples to be used for evidence
    P.number_samples_per_evidence=5; % for dynamic-switching-signal-blocked
    P.phi_O=[[0:2:2*(P.number_orientations-1)];[1:2:2*P.number_orientations-1]]*pi/2/P.number_orientations;
    P.pO=ones(1,P.number_orientations)/P.number_orientations;
    P.pL=ones(1,P.number_locations)/P.number_locations;
    P.tauStyle=1;
    P.sigmaStyle=1e-10;
    P.odds_inc=1/20; % with 80 samples/sec, every 250ms independent signal
    P.nT = 2; % Number of possible tasks
    
  otherwise
    warning('invalid option');
end