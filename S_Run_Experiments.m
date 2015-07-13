function varargout = S_Run_Experiments(fct)

% create simulations for coupling strength dependency figure
switch fct
  case 'contrast'
    for i=1:6, cond{7-i}=3*[i 0]; end
    cond{7}=[0 0];
    for i=1:6, cond{7+i}=2*[0 i]; end
    para=S_Exp_Para('paper-corr-performance')
    
    delta=[0 0.005 0.02 0.08];
    
    for i=1:length(delta), disp(['i, delta=' num2str([i delta(i)])]);
      para.delta=delta(i);
      for j=1:length(cond), disp(['j, c=' num2str([j cond{j}])]);
        para.stimulus_contrast=cond{j};
        if max(abs(cond{j}))==0
          para.prior_task=[0 1]; % oblique for corr matrix
        else
          para.prior_task=[1 0]; % cardinal for non-zero constrast
        end
        disp(['task prior: ' num2str(para.prior_task)]);
        e{i,j}=S_Experiment(para);
      end
      save('zwischenstand_S_Run_Experiment.mat','e','-v7.3');
    end
    varargout{1}=e;
    
  otherwise
    warning(fct);
end

