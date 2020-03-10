function [lshc_results, hslc_results] = S_Run_CB_Experiments(varargin)
%S_RUN_CB_EXPERIMENTS Wrapper to run and save CB results

% create simulations for coupling strength dependency figure
para = S_Exp_Para('cb');

savedir = pwd;
if isfield(para, 'savedir'), savedir = para.savedir; end

%% RUN LSHC
para.S.number_repetitions = 1000;
para.G.delta = .08;
para.I.stimulus_contrast = [+4 +4];
para.I.signal_match_probability = .9;
lshc_results = runWithParams(para, savedir);
CB_Diagnostics(lshc_results, 'LSHC Condition');

%% RUN HSLC
para.S.number_repetitions = 1000;
para.G.delta = .01;
para.I.stimulus_contrast = [+18 +18];
para.I.signal_match_probability = .6;
hslc_results = runWithParams(para, savedir);
CB_Diagnostics(hslc_results, 'HSLC Condition');
end

function results = runWithParams(para, savedir)
p_match = max(para.I.signal_match_probability, 1-para.I.signal_match_probability);

savename = sprintf('sampling_results_tr%d_c%.2f_pm%.2f_d%.3f.mat', ...
    para.S.number_repetitions, para.I.stimulus_contrast(1), p_match, para.G.delta);
savefile = fullfile(savedir, savename);

if exist(savefile, 'file')
    ld = load(savefile);
    results = ld.results;
else
    % Run once with the given parameters
    para.I.signal_match_probability = p_match;
    disp('=== Running Model 1/2 ===');
    tic;
    results{1} = S_Experiment(para);
    disp('=== Finished Model 1/2 ===');
    toc;
    % Run a second time, flipping which category is 'correct'
    para.I.signal_match_probability = 1 - p_match;
    disp('=== Running Model 2/2 ===');
    tic;
    results{2} = S_Experiment(para);
    disp('=== Finished Model 2/2 ===');
    toc;
    % Save results to skip computation later
    save(savefile, 'results', 'para');
end
end