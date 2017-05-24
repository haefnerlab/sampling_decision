function [optim_correct, optim_kappas, percent_correct, contrasts, match_probs, kappas] = ...
    plotSensoryCategorySpace(nGrid, params, postfix)

if ~exist('memo', 'dir'), mkdir('memo'); end
if nargin < 2
    params = S_Exp_Para('cb');
    postfix = 'defaults';
elseif nargin == 2
    error('if supplying params, also give a 3rd argument identifying it for memoization');
end

nProbs = min(nGrid, params.I.n_frames);
contrasts = linspace(0, 20, nGrid);
match_probs = linspace(0.5, 1, nProbs);
kappas = linspace(0, 3, nGrid);
[cc, pp, kk] = meshgrid(contrasts, match_probs, kappas);

%% Load or Run sampling model on each combination of contrast/prob/kappa

percent_correct = zeros(size(cc));

parfor i=1:numel(cc)
    fileName = ['CB-c' num2str(cc(i)) '-p' num2str(pp(i)) '-k' num2str(kk(i)) '-' postfix '.mat'];
    result = LoadOrRun(@S_Experiment, {params}, fullfile('memo', fileName));
    nTrials = size(result.O, 1);
    choices = arrayfun(@(t) mode(result.O(t, 1, :)), 1:nTrials);
    percent_correct(i) = sum(choices == 1) / nTrials;
end

% Maximize percent-correct across values of kappa.
% TODO - interpn before smoothing?
if exist('smoothn', 'file')
    percent_correct_smooth = smoothn(percent_correct);
else
    percent_correct_smooth = percent_correct;
    warning('Cannot find smoothn - %Correct results will be noisier');
end

optim_correct = zeros(nGrid, nProbs);
optim_kappas = zeros(nGrid, nProbs);

for i=1:nGrid
    for j=1:nProbs
        [max_pc, max_idx] = max(percent_correct_smooth(i, j, :));
        optim_correct(i, j) = max_pc;
        optim_kappas(i, j) = kappas(max_idx);
    end
end

%% Plot results

figure;
subplot(1,2,1);
imagesc(optim_correct, [0.5 1]);
axis image;
colorbar;
set(gca, 'YDir', 'Normal');
xlabel('contrast');
ylabel('p_{match}');
title('Percent Correct');

subplot(1,2,2);
imagesc(optim_kappas, [min(kappas) max(kappas)]);
axis image;
colorbar;
set(gca, 'YDir', 'Normal');
xlabel('contrast');
ylabel('p_{match}');
title('Optimized \kappa');

end