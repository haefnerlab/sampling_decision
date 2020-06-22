function [percent_correct, fig] = plotSensoryCategorySpace(contrasts, match_probs, params, postfix)

if ~exist('memo', 'dir') && ~exist('memo', 'file'), mkdir('memo'); end
if nargin < 3
    params = S_Exp_Para('cb');
    postfix = 'defaults';
    
    % Large number of parameters to consider.. run fewer trials per parameter
    params.S.number_repetitions = 200;
elseif nargin == 3
    error('if supplying params, also give a 3rd argument identifying it for memoization');
end

[cc, pp] = meshgrid(contrasts, match_probs);

%% Load or Run sampling model on each combination of contrast/prob/kappa

percent_correct = zeros(size(cc));

for i=1:numel(cc)
    % Dividing by ci_to_kappa(.9) ensures a match to other simulations, i.e. that at CI=0.9, kappa
    % is 1
    kappa = ci_to_kappa(pp(i)) / ci_to_kappa(.9);
    fileName = ['CB-c' num2str(cc(i)) '-p' num2str(pp(i)) '-k' num2str(kappa) '-' postfix '.mat'];
    disp(fileName);
    params_copy = params;
    params_copy.G.kappa_O = [kappa 0];
    params_copy.I.stimulus_contrast = cc(i)*[+1 +1];
    params_copy.I.signal_match_probability = pp(i);
    result = LoadOrRun(@S_Experiment, {params_copy}, fullfile('memo', fileName));
    choices = result.O(:,2,end) < .5;
    ambivalent = result.O(:,2,end) == .5;
    choices(ambivalent) = rand(sum(ambivalent), 1) < .5;
    percent_correct(i) = mean(choices == 1);
end

%% Plot results

if nargout >= 2
    fig = figure;
    imagesc(percent_correct, [0.5 1]);
    axis image;
    colorbar;
    set(gca, 'XTick', 1:length(contrasts), 'XTickLabel', contrasts);
    set(gca, 'YTick', 1:length(match_probs), 'YTickLabel', match_probs, 'YDir', 'normal');
    xlabel('contrast');
    ylabel('p_{match}');
    title('Percent Correct');
end

end