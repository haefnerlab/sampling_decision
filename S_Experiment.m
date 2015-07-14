function out = S_Experiment()
P = S_Exp_Para('paper-2AFC-corr');
% P.fct:
% 2x2: two locations, two orientations (NIPS)
% 1xN: 1 location, many orientations
% nxN: n locations, many orientations (Current!)
if P.number_orientations<0, 
    task='detection'; 
    P.number_orientations=2; 
else task='discrimination'; end

    % Generative model
    P.nx=Get_Varargin;
    if isempty(P.nx),       P.nx=32; end
    p=C_Projection('nxN',P.nx,P.dimension_X,P.dimension_G,P.number_locations);
    P.fct='t-l-op-g-s';
    P.task=task;
    P.phi_x = p.phi_x;
    P.phi_g = p.phi_g;
    P.G = p.G;
    P.x = p.x;
    P.y=p.y;
    P.ny=p.ny;
    
    % Sampling setup
    P.n_samples=P.number_burn_in+P.number_samples_to_use;
    
    
    %Split P into Generative Model G, Sampling parameters S, and Input I.
    
    G.number_orientations = P.number_orientations;
    G.prior_task = P.prior_task;
    G.number_locations = P.number_locations
    G.dimension_X = P.dimension_X;
    G.dimension_G = P.dimension_G;
    G.ny = P.ny;
    G.nx = P.nx;
    G.kappa_O = P.kappa_O;
    G.kappa_G = P.kappa_G;
    I.stimulus_regime = P.stimulus_regime;
    I.stimulus_contrast = P.stimulus_contrast;
    S.number_repetitions = P.number_repetitions;
    S.zero_signal = P.n_zero_signal;
    S.number_burn_in = P.number_burn_in;
    S.number_samples_to_use = P.number_samples_to_use;
    S.number_samples_per_evidence = P.number_samples_per_evidence;
    S.phi_O = P.phi_O;
    S.pO = P.pO;
    S.pL= P.pL;
    S.nT = P.nT;
    S.tauStyle = P.tauStyle;
    S.sigmaStyle = P.sigmaStyle;
    G.delta = P.delta;
    G.phi_x = P.phi_x;
    G.task = P.task;
    S.odds_inc = P.odds_inc;
    G.phi_g = P.phi_g;
    G.G = P.G;
    I.x = P.x;
    I.y = P.y;
    S.n_samples = P.n_samples;
    
    
    
    
    
    
    
    % Input
    I.fct='nx2';
    I.n_zero_signal=S.zero_signal;
    switch I.stimulus_regime
      case {'static','blank'}
        I.n_frames=1;
        S.access=ones(1,S.n_samples);
      case 'static-delayed'
        I.n_frames=2;
        S.access=ones(1,S.n_samples);
        S.access(I.n_zero_signal+1:end)=2;
      case 'dynamic-delayed'
        I.n_frames=S.n_samples;
        if I.n_zero_signal>=I.n_frames
          error('n_zero_signal>=n_frames!');
        end
        S.access=1:S.n_samples;
      case 'dynamic-switching-signal'
        I.n_frames=S.n_samples;
        if I.n_zero_signal>=I.n_frames
          error('n_zero_signal>=n_frames!');
        end
        S.access=1:S.n_samples;
      case 'dynamic-switching-signal-blocked'
        zs=I.n_zero_signal+1; ns=S.n_samples; spe=S.number_samples_per_evidence;
        S.access=[1:zs zs+ceil((1:ns-zs)/spe)];
        I.n_frames=numel(unique(S.access));
        if I.n_zero_signal>=I.n_frames
          error('n_zero_signal>=n_frames!');
        end
    end
    %zs=I.n_zero_signal+1;
    %ns=P.n_samples;
    %spe=P.number_samples_per_evidence;
    I.n_frames=numel(unique(S.access));
    if max(abs(I.stimulus_contrast))>0 && G.prior_task(1)<G.prior_task(2)
      warning('non-zero stim contrasts imply Task=1, not 2!!');
    end
    % Sanity checks
    if G.kappa_O(2)>=G.kappa_O(1), warning('are you sure about P.kappa?'); end

    % Repeating Gibbs sampling P.number_repetitions times
    out.X=zeros(S.number_repetitions,G.number_locations*G.dimension_X, S.n_samples);

    % the below two lines are my way of dealing with Matlab's lack of
    % macros. use comments to switch between serial and parallel processing

    %PARALLEL SUPPORT
    %parfor i=1:P.number_repetitions,
    warning('serial!!'); for i=1:S.number_repetitions, ProgressReport(10,i-1,S.number_repetitions);

    %parfor i=1:P.number_repetitions
    %warning('serial!!'); for i=1:P.number_repetitions, ProgressReport(10,i-1,P.number_repetitions);

      if mod(i,20)==0
        disp(['Computing Repetition ' num2str(i) ' / ' num2str(S.number_repetitions)]);
      end
      switch I.stimulus_regime
        case 'static'
          signal=I.stimulus_contrast;
        case 'static-delayed'
          signal(1,:,:)=zeros(size(I.stimulus_contrast));
          signal(2,:,:)=I.stimulus_contrast;
        case 'dynamic-delayed'
          if G.number_locations>1, error('not implemented'); end
          signal=zeros(I.n_frames,1,2);
          for k=I.n_zero_signal+1:I.n_frames
            signal(k,1,:)=I.stimulus_contrast;
          end
        case {'dynamic-switching-signal','dynamic-switching-signal-blocked'}
          if G.number_locations>1, error('not implemented'); end
          signal=zeros(I.n_frames,G.number_locations,2);
          for k=I.n_zero_signal+1:I.n_frames
            on=1+binornd(1,0.5);
            signal(k,1,on)=I.stimulus_contrast(on);
          end
      end
      switch I.stimulus_regime
        case 'blank'
          Y=zeros(G.ny,G.number_locations*G.ny);
        otherwise
          Y=InputImage(I.fct,G.number_locations,G.ny,signal);
      end
      
      %Perform a trial 
      [aux_X{i} aux_G{i} aux_O{i} aux_L{i} aux_S{i} aux_T{i}]=...
        Sampling_Gibbs_InPlace_Fast(G, S, I ,Y,S.access,I.n_zero_signal);

     for j=1:G.dimension_X*G.number_locations
        switch I.stimulus_regime
          case {'static','blank'}
            Signal{i}{j}=mean(mean(Y.*reshape(G.G(:,j),G.ny,G.nx)));
          case {'static-delayed','dynamic-delayed','dynamic-switching-signal','dynamic-switching-signal-blocked'}
            for k=1:I.n_frames
              Signal{i}{j,k}=mean(mean(squeeze(Y(k,:,:)).*...
                reshape(G.G(:,j),G.ny,G.nx)));
            end
          otherwise, error(I.stimulus_regime);
        end
      end
    end
    out.Signal=zeros(S.number_repetitions,G.dimension_X*G.number_locations,I.n_frames);
    for i=S.number_repetitions:-1:1
      %if SAVE_Y, out.Y(i,:,:,:)=Y{i}; end
      out.X(i,:,:)  =aux_X{i};
      out.G(i,:,:,:)=aux_G{i};
      out.O(i,:,:)  =aux_O{i};
      out.L(i,:,:)  =aux_L{i};
      out.T(i,:,:)  =aux_T{i};
      out.S(i,:)    =aux_S{i};
      for j=G.dimension_X*G.number_locations:-1:1
        switch I.stimulus_regime
          case {'static','blank'}
            out.Signal(i,j)=Signal{i}{j};
          case {'static-delayed','dynamic-delayed','dynamic-switching-signal','dynamic-switching-signal-blocked'}
            for k=I.n_frames:-1:1
              out.Signal(i,j,k)=Signal{i}{j,k};
            end
          otherwise, error(I.stimulus_regime);
        end
      end
    end
    
    %Convert back to single variable P for output, thne perform backwards
    %compatibility handling 
    
    P.number_orientations = G.number_orientations;
    P.prior_task = G.prior_task;
    P.number_locations = G.number_locations
    P.dimension_X = G.dimension_X;
    P.dimension_G = G.dimension_G;
    P.ny = G.ny;
    P.nx = G.nx;
    P.kappa_O = G.kappa_O;
    P.kappa_G = G.kappa_G;
    P.stimulus_regime = I.stimulus_regime;
    P.stimulus_contrast = I.stimulus_contrast;
    P.number_repetitions = S.number_repetitions;
    P.n_zero_signal = S.zero_signal;
    P.number_burn_in = S.number_burn_in;
    P.number_samples_to_use = S.number_samples_to_use;
    P.number_samples_per_evidence = S.number_samples_per_evidence;
    P.phi_O = S.phi_O;
    P.pO = S.pO;
    P.pL= S.pL;
    P.nT = S.nT;
    P.tauStyle = S.tauStyle;
    P.sigmaStyle = S.sigmaStyle;
    P.delta = G.delta;
    P.phi_x = G.phi_x;
    P.task = G.task;
    P.odds_inc = S.odds_inc;
    P.phi_g = G.phi_g;
    P.G = G.G;
    P.x = I.x;
    P.y = I.y;
    P.n_samples = S.n_samples;
    
    out.Projection=P; out.InputImage=I; out.Sampling=S;

    out = BackwardsComp(out);

end
