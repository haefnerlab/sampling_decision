function [lshc_results, hslc_results] = S_Run_CB_Experiments(varargin)
%S_RUN_CB_EXPERIMENTS Wrapper to run and save CB results

% create simulations for coupling strength dependency figure
para = S_Exp_Para('cb');
cs_para = para; % Make this copy before restricting # of seeds so that CS-space plot uses iid trials

% correlations analysis requires many repeats from small number of seeds
nSeeds = 20;
rng('shuffle');
para.I.frame_seq_seeds = randi(1e9, 1, nSeeds);

savedir = pwd;
if isfield(para, 'savedir'), savedir = para.savedir; end

%% Run Category-Sensory Space

contrasts = 0.5:0.5:10;
match_probs = .51:.02:.99;
cs_para.S.number_repetitions = 200; % n trials per point in C-S space. may be small because we know the space is smooth
p_correct = LoadOrRun(@plotSensoryCategorySpace, ...
    {contrasts, match_probs, cs_para, 'defaults'}, ...
    fullfile(savedir, 'cache-SC-space.mat'));

smooth_p_correct = smoothn(p_correct);
cs_fig = figure;
imagesc('XData', contrasts, 'YData', match_probs, 'CData', smooth_p_correct, [.5 1]);
axis tight; axis square; set(gca, 'YDir', 'normal');
[cc, pp] = meshgrid(contrasts, match_probs);
hold on;
contour(cc, pp, smooth_p_correct, [0.7 0.7], '-w', 'LineWidth', 2);
xlabel('Sensory Information (contrast)');
ylabel('Category Information (frame category match)');

%% Set LSHC and HSLC parameters

% [contrast, p_match]
lshc_si_ci = [4, .90];
hslc_si_ci = [9, .61];

% Convention is to adjust the "true" kappa for a given CI s.t. kappa = 1 at CI = 0.9 (see
% @plotSensoryCategorySpace)
lshc_kappa = ci_to_kappa(lshc_si_ci(2)) / ci_to_kappa(.9);
hslc_kappa = ci_to_kappa(hslc_si_ci(2)) / ci_to_kappa(.9);

% Overlay these points on the performance curve above
plot(lshc_si_ci(1), lshc_si_ci(2), 'o', 'Color', [.6 0 0], 'MarkerFaceColor', [.6 0 0]);
plot(hslc_si_ci(1), hslc_si_ci(2), 'o', 'Color', [0 0 .6], 'MarkerFaceColor', [0 0 .6]);

%% RUN LSHC TASK
para.S.number_repetitions = 2500;
para.G.kappa_O = [lshc_kappa 0]; % attended and unattended
para.I.stimulus_contrast = lshc_si_ci(1)*[1 1];
para.I.signal_match_probability = lshc_si_ci(2);
lshc_results = runWithParams(para, savedir);
CB_Diagnostics(lshc_results, 'LSHC Condition');

%% RUN HSLC TASK
para.S.number_repetitions = 2500;
para.G.kappa_O = [hslc_kappa 0]; % attended and unattended
para.I.stimulus_contrast = hslc_si_ci(1)*[1 1];
para.I.signal_match_probability = hslc_si_ci(2);
hslc_results = runWithParams(para, savedir);
CB_Diagnostics(hslc_results, 'HSLC Condition');
end

function results = runWithParams(para, savedir)
savename = sprintf('sampling_results_tr%d_c%.2f_pm%.2f_k%s.mat', ...
    para.S.number_repetitions, para.I.stimulus_contrast(1), para.I.signal_match_probability, num2str(para.G.kappa_O(1)));
savefile = fullfile(savedir, savename);

if exist(savefile, 'file')
    ld = load(savefile);
    results = ld.results;
else
    % Run once with the given parameters
    p_positive = para.I.signal_match_probability;
    p_negative = 1-para.I.signal_match_probability;
    para.I.signal_match_probability = p_positive;
    disp('=== Running Model 1/2 ===');
    tic;
    results{1} = S_Experiment(para);
    disp('=== Finished Model 1/2 ===');
    toc;
    % Run a second time, flipping which category is 'correct'
    para.I.signal_match_probability = p_negative;
    disp('=== Running Model 2/2 ===');
    tic;
    results{2} = S_Experiment(para);
    disp('=== Finished Model 2/2 ===');
    toc;
    % Save results to skip computation later
    save(savefile, 'results', 'para');
end
end