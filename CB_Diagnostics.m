function CB_Diagnostics(results, name)

results_pos = results{1};
results_neg = results{2};
% Assuming that the two structures really are comparable -- see @S_Run_CB_Experiments

zs = results_pos.InputImage.n_zero_signal;
nX = results_pos.Projection.dimension_X;

choice_pos = results_pos.O(:,3,end) > .5;
choice_neg = results_neg.O(:,3,end) > .5;
n_pos = length(choice_pos);
n_neg = length(choice_neg);
percent_correct = (mean(choice_pos == true)*n_pos + mean(choice_neg == false)*n_neg) / (n_pos + n_neg);

fprintf('=== CB Diagnostics: %s ===\n', name);
fprintf('\t%d trials\n', n_pos);
fprintf('\tcontrast set to %.1f, %.1f\n', results_pos.InputImage.c(1), results_pos.InputImage.c(2));
fprintf('\tkappa set to %.3f\n', results_pos.Projection.kappa_O(1));
fprintf('\tdelta set to %.3f\n', results_pos.Projection.delta);
fprintf('\t%.1f%% correct overall\n', 100*percent_correct);

figure;

%% Extract 'true' frame categories per frame
pool_projection = -cos(linspace(0, 2*pi, nX+1));
pool_projection = pool_projection(1:end-1);
sig_frames = zs+2:size(results_pos.Signal, 3);
signal_pos = squeeze(sum(results_pos.Signal(:,:,sig_frames) .* reshape(pool_projection, 1, nX, 1), 2));
signal_neg = squeeze(sum(results_neg.Signal(:,:,sig_frames) .* reshape(pool_projection, 1, nX, 1), 2));

% Get true frame categories, which is an index, i.e. one of [1 2]. Convert to [-1 +1].
true_cat = [results_pos.FrameCategory(:, zs+1:end); results_neg.FrameCategory(:, zs+1:end)];

%% PK analysis

w = glmfit([signal_pos; signal_neg], [choice_pos; choice_neg], 'binomial');
subplot(1,2,1);
w_norm = w(2:end) / mean(w(2:end));
plot(w_norm);
ylim([0 max(w_norm)]);
title('PK');

%% Correlations analysis
minTrialsPer = 10;
allX = cat(1, results_pos.X, results_neg.X);
allX = sum(allX(:, :, zs+2:end), 3);
% TODO: Select top 2ish repeates
% Consider trials 'frozen' post-hoc when they had the same sequence of frame categories.
[uSigs, ~, idxExpand] = unique(true_cat, 'rows');
fprintf('\t%d unique combinations of categories\n', size(uSigs, 1));
nRepeats = arrayfun(@(i) sum(idxExpand == i), 1:size(uSigs, 1));
fprintf('\t%d sequences with > %d repeats\n', sum(nRepeats > minTrialsPer), minTrialsPer);
zscoreX = nan(size(allX));
for iSig=size(uSigs, 1):-1:1
    trialsAt = idxExpand == iSig;
    % Don't bother with cases where there were fewer than 5 trials with the same frame pattern
    if nRepeats(iSig) > minTrialsPer
        zscoreX(trialsAt, :) = zscore(allX(trialsAt, :));
    end
end
corrX = nancov(zscoreX);
subplot(1,2,2);
imagesc(corrX - eye(size(corrX)));
axis image;
colorbar;
title('Correlations');

end