function CB_Diagnostics(results)

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

fprintf('=== CB Diagnostics ===\n');
fprintf('\t%d trials\n', n_pos);
fprintf('\tcontrast set to %.1f, %.1f\n', results_pos.InputImage.c(1), results_pos.InputImage.c(2));
fprintf('\tdelta set to %.3f\n', results_pos.Projection.delta);
fprintf('\t%.1f%% correct overall\n', 100*percent_correct);

figure;

%% Extract 'true' frame categories per frame
pool_projection = -cos(linspace(0, 2*pi, nX+1));
pool_projection = pool_projection(1:end-1);
sig_frames = zs+2:size(results_pos.Signal, 3);
signal_pos = squeeze(sum(results_pos.Signal(:,:,sig_frames) .* reshape(pool_projection, 1, nX, 1), 2));
signal_neg = squeeze(sum(results_neg.Signal(:,:,sig_frames) .* reshape(pool_projection, 1, nX, 1), 2));

true_cat = sign([signal_pos; signal_neg]);

%% PK analysis

w = glmfit([signal_pos; signal_neg], [choice_pos; choice_neg], 'binomial');
subplot(1,2,1);
plot(w(2:end));
title('PK');

%% Correlations analysis
allX = cat(1, results_pos.X, results_neg.X);
allX = sum(allX(:, :, zs+2:end), 3);
% Consider trials 'frozen' post-hoc when they had the same sequence of frame categories.
[uSigs, ~, idxExpand] = unique(true_cat, 'rows');
muXPerCond = nan(size(uSigs, 1), nX);
for iSig=length(uSigs):-1:1
    trialsAt = idxExpand == iSig;
    % Don't bother with cases where there were fewer than 5 trials with the same frame pattern
    if sum(trialsAt) > 5
        muXPerCond(iSig, :) = mean(allX(trialsAt, :), 1);
    end
end
zeroX = allX - muXPerCond(idxExpand, :);
validX = ~any(isnan(zeroX), 2);
zeroX = zeroX(validX, :);
moment2X = zeroX'*zeroX;
stdX = sqrt(diag(moment2X));
corrX = moment2X ./ (stdX .* stdX');
subplot(1,2,2);
imagesc(corrX - eye(size(moment2X)));
axis image;
colorbar;
title('Correlations');

end