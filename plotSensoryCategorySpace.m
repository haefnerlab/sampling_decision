function [percent_correct, contrasts, match_probs] = plotSensoryCategorySpace(nGrid, params, postfix)

if ~exist('memo', 'dir'), mkdir('memo'); end
if nargin < 2
    params = S_Exp_Para('cb');
    postfix = 'defaults';
elseif nargin == 2
    error('if supplying params, also give a 3rd argument identifying it for memoization');
end

% Large number of parameters to consider.. run only 100 trials each
params.S.number_repetitions = 100;

contrasts = linspace(0, 10, nGrid);
match_probs = linspace(0.5, 1, nGrid);
[cc, pp] = meshgrid(contrasts, match_probs);

%% Load or Run sampling model on each combination of contrast/prob/kappa

percent_correct = zeros(size(cc));

parfor i=1:numel(cc)
    kappa = ci_to_kappa(cc(i));
    fileName = ['CB-c' num2str(cc(i)) '-p' num2str(pp(i)) '-k' num2str(kappa) '-' postfix '.mat'];
    params_copy = params;
    params_copy.G.kappa_O = [kappa 0];
    params_copy.I.stimulus_contrast = cc(i)*[+1 +1];
    params_copy.I.signal_match_probability = pp(i);
    result = LoadOrRun(@S_Experiment, {params_copy}, fullfile('memo', fileName));
    choices = result.O(:,3,end) > .5;
    percent_correct(i) = mean(choices == 1);
end

%% Plot results

figure;
imagesc('XData', contrasts, 'YData', match_probs, 'CData', percent_correct, [0.5 1]);
axis image;
colorbar;
set(gca, 'YDir', 'Normal');
xlabel('contrast');
ylabel('p_{match}');
title('Percent Correct');

end

function kappa = ci_to_kappa(ci)
os = linspace(0,pi,1001);
os = os(1:end-1);
ks = linspace(0, 6);
for ik=length(ks):-1:1
    p1 = exp(ks(ik)*cos(os*2)); p1 = p1 ./ sum(p1);
    p2 = exp(-ks(ik)*cos(os*2)); p2 = p2 ./ sum(p2);
    ratio = p1 ./ (p1 + p2);
    emp_ci(ik) = dot(p1, ratio);
end

ks = [0 ks 7];
emp_ci = [.5 emp_ci 1];

kappa = interp1(emp_ci, ks, ci);
end