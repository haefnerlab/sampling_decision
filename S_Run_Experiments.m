function varargout = S_Run_Experiments(fct,varargin)

% create simulations for coupling strength dependency figure
switch fct
  case 'contrast'
    stepsize=3; nsteps=6; % computing 2*nsteps+1 contrast conditions!
    for i=1:nsteps, cond{nsteps+1-i}=stepsize*[i 0]; end
    cond{nsteps+1}=[0 0];
    for i=1:nsteps, cond{nsteps+1+i}=stepsize*[0 i]; end
    para=S_Exp_Para('paper-corr-performance')
    
    delta=[0 0.005 0.02 0.08];
    
    if nargin>1, e=varargin{1};
    else
      disp('initializing a completely new array of experiments');
      e=cell(length(delta),length(cond));
    end
    
    added_something=0;
    for i=1:length(delta), disp(['i, delta=' num2str([i delta(i)])]);
      para.G.delta=delta(i);
      for j=1:length(cond), disp(['j, c=' num2str([j cond{j}])]);
        if size(e,1)<i || size(e,2)<j || isempty(e{i,j})
          disp('computing since entry empty');
          para.I.stimulus_contrast=cond{j};
          if max(abs(cond{j}))==0
            para.G.prior_task=[1 0]; % not oblique for corr matrix
          else
            para.G.prior_task=[1 0]; % cardinal for non-zero constrast
          end
          disp(['task prior: ' num2str(para.G.prior_task)]);
          e{i,j}=S_Experiment(para);
          added_something=1;
        else
          disp('skipping since entry exists already');
        end
      end
      if added_something
        disp('saving intermediate results...');
        save('zwischenstand_S_Run_Experiment.mat','e','-v7.3');
      end
    end
    varargout{1}=e;
    
  otherwise
    warning(fct);
end

