function varargout = S_Run_Experiments(fct,varargin)
%S_RUN_EXPERIMENTS Wrapper to run and save results of S_Experiment at
%   multiple contrast levels
%
%       e = S_RUN_EXPERIMENTS('contrast', [para, [e]]) Run S_Experiment
%       at 2*para.nsteps+1 contrast values (symmetric around zero). If some
%       results have been previously computed, they may be passed in as
%       'e'. The return value 'e' is a single struct containing results.
%
%       para is a struct from S_Exp_Para, and defaults to the
%       'paper-corr-performance' condition
%
%       para.nsteps (default 6) and para.stepsize (default 3) determine the
%       contrast values used.
%
%       Currently 'contrast' is the only implemented experiment.

% create simulations for coupling strength dependency figure
if nargin > 1
    para = varargin{1};
else
    para = S_Exp_Para('paper-corr-performance');
end

switch fct
    
    case 'contrast'
        if ~isfield(para, 'nsteps'), para.nsteps = 6; end
        if ~isfield(para, 'stepsize'), para.stepsize = 3; end
        % computing (2 * para.nsteps + 1) contrast conditions
        range = para.stepsize * para.nsteps;
        nconditions = 2 * para.nsteps + 1;
        condmat = [linspace(range, -range, nconditions);
            linspace(-range, range, nconditions)]';
        condmat(condmat < 0) = 0;
        cond = mat2cell(condmat, ones(nconditions,1), 2);
        
        delta=[0 0.005 0.02 0.08];
        
        if nargin>2, e=varargin{2};
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
end

