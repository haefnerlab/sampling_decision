function out = S_Experiment(P)

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
    n_samples=P.number_burn_in+P.number_samples_to_use;
    
    % Input
    I.fct='nx2';
    I.n_zero_signal=20; % sometimes 50. Why?
    switch P.stimulus_regime
      case {'static','blank'}
        I.n_frames=1;
        S.access=ones(1,n_samples);
      case 'static-delayed'
        I.n_frames=2;
        S.access=ones(1,n_samples);
        S.access(I.n_zero_signal+1:end)=2;
      case 'dynamic-delayed'
        I.n_frames=n_samples;
        if I.n_zero_signal>=I.n_frames
          error('n_zero_signal>=n_frames!');
        end
        S.access=1:n_samples;
      case 'dynamic-switching-signal'
        I.n_frames=n_samples;
        if I.n_zero_signal>=I.n_frames
          error('n_zero_signal>=n_frames!');
        end
        S.access=1:n_samples;
      case 'dynamic-switching-signal-blocked'
        zs=I.n_zero_signal+1; ns=n_samples; spe=P.number_samples_per_evidence;
        S.access=[1:zs zs+ceil((1:ns-zs)/spe)];
        I.n_frames=numel(unique(S.access));
        if I.n_zero_signal>=I.n_frames
          error('n_zero_signal>=n_frames!');
        end
    end
    %zs=I.n_zero_signal+1;
    %ns=n_samples;
    %spe=P.number_samples_per_evidence;
    I.n_frames=numel(unique(S.access));
    if max(abs(P.stimulus_contrast))>0 && P.P.prior_task(1)<P.P.prior_task(2)
      warning('non-zero stim contrasts imply Task=1, not 2!!');
    end
    % Sanity checks
    if P.kappa_O(2)>=P.kappa_O(1), warning('are you sure about P.kappa?'); end
    % Repeating Gibbs sampling P.number_repetitions times
    out.X=zeros(P.number_repetitions,P.number_locations*P.dimension_X, n_samples);
    % the below two lines are my way of dealing with Matlab's lack of
    % macros. use comments to switch between serial and parallel processing
    parfor i=1:P.number_repetitions,
    %warning('serial!!'); for i=1:P.number_repetitions, ProgressReport(10,i-1,P.number_repetitions);
      if mod(i,20)==0
        disp(['Computing Repetition ' num2str(i) ' / ' num2str(P.number_repetitions)]);
      end
      switch P.stimulus_regime
        case 'static'
          signal=P.stimulus_contrast;
        case 'static-delayed'
          signal(1,:,:)=zeros(size(P.stimulus_contrast));
          signal(2,:,:)=P.stimulus_contrast;
        case 'dynamic-delayed'
          if P.number_locations>1, error('not implemented'); end
          signal=zeros(I.n_frames,1,2);
          for k=I.n_zero_signal+1:I.n_frames
            signal(k,1,:)=P.stimulus_contrast;
          end
        case {'dynamic-switching-signal','dynamic-switching-signal-blocked'}
          if P.number_locations>1, error('not implemented'); end
          signal=zeros(I.n_frames,P.number_locations,2);
          for k=I.n_zero_signal+1:I.n_frames
            on=1+binornd(1,0.5);
            signal(k,1,on)=P.stimulus_contrast(on);
          end
      end
      switch P.stimulus_regime
        case 'blank'
          Y=zeros(P.ny,P.number_locations*P.ny);
        otherwise
          Y=InputImage(I.fct,P.number_locations,P.ny,signal);
      end
      
      %Perform a trial 
      [aux_X{i} aux_G{i} aux_O{i} aux_L{i} aux_S{i} aux_T{i}]=...
        Sampling_Gibbs_InPlace_Fast(P,Y,n_samples,S.access,I.n_zero_signal);
     for j=1:P.dimension_X*P.number_locations
        switch P.stimulus_regime
          case {'static','blank'}
            Signal{i}{j}=mean(mean(Y.*reshape(P.G(:,j),P.ny,P.nx)));
          case {'static-delayed','dynamic-delayed','dynamic-switching-signal','dynamic-switching-signal-blocked'}
            for k=1:I.n_frames
              Signal{i}{j,k}=mean(mean(squeeze(Y(k,:,:)).*...
                reshape(P.G(:,j),P.ny,P.nx)));
            end
          otherwise, error(P.stimulus_regime);
        end
      end
    end
    out.Signal=zeros(P.number_repetitions,P.dimension_X*P.number_locations,I.n_frames);
    for i=P.number_repetitions:-1:1
      %if SAVE_Y, out.Y(i,:,:,:)=Y{i}; end
      out.X(i,:,:)  =aux_X{i};
      out.G(i,:,:,:)=aux_G{i};
      out.O(i,:,:)  =aux_O{i};
      out.L(i,:,:)  =aux_L{i};
      out.T(i,:,:)  =aux_T{i};
      out.S(i,:)    =aux_S{i};
      for j=P.dimension_X*P.number_locations:-1:1
        switch P.stimulus_regime
          case {'static','blank'}
            out.Signal(i,j)=Signal{i}{j};
          case {'static-delayed','dynamic-delayed','dynamic-switching-signal','dynamic-switching-signal-blocked'}
            for k=I.n_frames:-1:1
              out.Signal(i,j,k)=Signal{i}{j,k};
            end
          otherwise, error(P.stimulus_regime);
        end
      end
    end
    out.Projection=P; out.InputImage=I; out.Sampling=S;
end
