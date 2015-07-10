function varargout = S_Run_Experiments(fct)

% create simulations for coupling strength dependency figure
switch fct
  case 'contrast'
    for i=1:11, cond{i}=[i-1 0]; end
    para=S_Exp_Para('test-contrast');
    delta=linspace(0,para.delta,4);
    
    for i=1:length(delta), disp(['i=' num2str(i)]);
      para.delta=delta(i);
      for j=1:length(cond), disp(['j=' num2str(j)]);
        para.stimulus_contrast=cond{j};
        e{i,j}=S_Experiment(para);
      end
    end
    varargout{i}=e;
    
  otherwise
    warning(fct);
end

