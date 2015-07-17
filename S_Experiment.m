function out = S_Experiment()
P = S_Exp_Para('paper-2AFC-corr');
% 2x2: two locations, two orientations (NIPS)
% 1xN: 1 location, many orientations
% nxN: n locations, many orientations (Current!)
if P.G.number_orientations<0, 
    task='detection'; 
    P.G.number_orientations=2; 
else task='discrimination'; end

    % Generative model
    P.G.nx=Get_Varargin;
    if isempty(P.G.nx),       P.G.nx=32; end
    p=C_Projection('nxN',P.G.nx,P.G.dimension_X,P.G.dimension_G,P.G.number_locations);
    P.fct='t-l-op-g-s';
    P.G.task=task;
    P.G.phi_x = p.phi_x;
    P.G.phi_g = p.phi_g;
    P.G.G = p.G;
    P.I.x = p.x;
    P.I.y=p.y;
    P.G.ny=p.ny;
    
    % Sampling setup
    P.S.n_samples=P.S.number_burn_in+P.S.number_samples_to_use;
    

    
    
    
    % Input
    P.I.fct='nx2';
    switch P.I.stimulus_regime
      case {'static','blank'}
        P.I.n_frames=1;
        P.S.access=ones(1,P.S.n_samples);
      case 'static-delayed'
        P.I.n_frames=2;
        P.S.access=ones(1,P.S.n_samples);
        P.S.access(P.I.n_zero_signal+1:end)=2;
      case 'dynamic-delayed'
        P.I.n_frames=P.S.n_samples;
        if P.I.n_zero_signal>=P.I.n_frames
          error('n_zero_signal>=n_frames!');
        end
        P.S.access=1:P.S.n_samples;
      case 'dynamic-switching-signal'
        P.I.n_frames=P.S.n_samples;
        if P.I.n_zero_signal>=P.I.n_frames
          error('n_zero_signal>=n_frames!');
        end
        P.S.access=1:P.S.n_samples;
      case 'dynamic-switching-signal-blocked'
        zs=P.I.n_zero_signal+1; ns=P.S.n_samples; spe=P.G.number_samples_per_evidence;
        P.S.access=[1:zs zs+ceil((1:ns-zs)/spe)];
        P.I.n_frames=numel(unique(P.S.access));
        if P.I.n_zero_signal>=P.I.n_frames
          error('n_zero_signal>=n_frames!');
        end
    end
    %zs=P.I.n_zero_signal+1;
    %ns=P.n_samples;
    %spe=P.number_samples_per_evidence;
    P.I.n_frames=numel(unique(P.S.access));
    if max(abs(P.I.stimulus_contrast))>0 && P.G.prior_task(1)<P.G.prior_task(2)
      warning('non-zero stim contrasts imply Task=1, not 2!!');
    end
    % Sanity checks
    if P.G.kappa_O(2)>=P.G.kappa_O(1), warning('are you sure about P.kappa?'); end

    % Repeating Gibbs sampling P.number_repetitions times
    out.X=zeros(P.S.number_repetitions,P.G.number_locations*P.G.dimension_X, P.S.n_samples);

    % the below two lines are my way of dealing with Matlab's lack of
    % macros. use comments to switch between serial and parallel processing

    %PARALLEL SUPPORT
    %parfor i=1:P.number_repetitions,
    warning('serial!!'); for i=1:P.S.number_repetitions, ProgressReport(10,i-1,P.S.number_repetitions);

    %parfor i=1:P.number_repetitions
    %warning('serial!!'); for i=1:P.number_repetitions, ProgressReport(10,i-1,P.number_repetitions);

      if mod(i,20)==0
        disp(['Computing Repetition ' num2str(i) ' / ' num2str(P.S.number_repetitions)]);
      end
      switch P.I.stimulus_regime
        case 'static'
          signal=P.I.stimulus_contrast;
        case 'static-delayed'
          signal(1,:,:)=zeros(size(P.I.stimulus_contrast));
          signal(2,:,:)=P.I.stimulus_contrast;
        case 'dynamic-delayed'
          if P.G.number_locations>1, error('not implemented'); end
          signal=zeros(P.I.n_frames,1,2);
          for k=P.I.n_zero_signal+1:P.I.n_frames
            signal(k,1,:)=P.I.stimulus_contrast;
          end
        case {'dynamic-switching-signal','dynamic-switching-signal-blocked'}
          if P.G.number_locations>1, error('not implemented'); end
          signal=zeros(P.I.n_frames,P.G.number_locations,2);
          for k=P.I.n_zero_signal+1:P.I.n_frames
            on=1+binornd(1,0.5);
            signal(k,1,on)=P.I.stimulus_contrast(on);
          end
      end
      switch P.I.stimulus_regime
        case 'blank'
          Y=zeros(P.G.ny,P.G.number_locations*P.G.ny);
        otherwise
          Y=InputImage(P.I.fct,P.G.number_locations,P.G.ny,signal);
      end
      
      %Perform a trial 
      [aux_X{i} aux_G{i} aux_O{i} aux_L{i} aux_S{i} aux_T{i}]=...
        Sampling_Gibbs_InPlace_Fast(P.G, P.S, P.I ,Y);

     for j=1:P.G.dimension_X*P.G.number_locations
        switch P.I.stimulus_regime
          case {'static','blank'}
            Signal{i}{j}=mean(mean(Y.*reshape(P.G.G(:,j),P.G.ny,P.G.nx)));
          case {'static-delayed','dynamic-delayed','dynamic-switching-signal','dynamic-switching-signal-blocked'}
            for k=1:P.I.n_frames
              Signal{i}{j,k}=mean(mean(squeeze(Y(k,:,:)).*...
                reshape(P.G.G(:,j),P.G.ny,P.G.nx)));
            end
          otherwise, error(P.I.stimulus_regime);
        end
      end
    end
    out.Signal=zeros(P.S.number_repetitions,P.G.dimension_X*P.G.number_locations,P.I.n_frames);
    for i=P.S.number_repetitions:-1:1
      %if SAVE_Y, out.Y(i,:,:,:)=Y{i}; end
      out.X(i,:,:)  =aux_X{i};
      out.G(i,:,:,:)=aux_G{i};
      out.O(i,:,:)  =aux_O{i};
      out.L(i,:,:)  =aux_L{i};
      out.T(i,:,:)  =aux_T{i};
      out.S(i,:)    =aux_S{i};
      for j=P.G.dimension_X*P.G.number_locations:-1:1
        switch P.I.stimulus_regime
          case {'static','blank'}
            out.Signal(i,j)=Signal{i}{j};
          case {'static-delayed','dynamic-delayed','dynamic-switching-signal','dynamic-switching-signal-blocked'}
            for k=P.I.n_frames:-1:1
              out.Signal(i,j,k)=Signal{i}{j,k};
            end
          otherwise, error(P.I.stimulus_regime);
        end
      end
    end
    
    %Convert back to single variable P for output, thne perform backwards
    %compatibility handling 

    
    out.Projection=P; out.InputImage=P.I; out.Sampling=P.S;

    out = BackwardsComp(out);

end
