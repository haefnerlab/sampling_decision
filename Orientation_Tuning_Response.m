function [ orient_responses ] = Orientation_Tuning_Response( Ge, S, n_orient, trials )
%ORIENTATION_TUNING_RESPONSE given model params Ge, presents it with 0:180
%orientations ('n_orient' of them), each 'trials' times
%
% returns a struct array that has the scalar field 'orientation' plus one
% matrix for each variable in the model (dim x samples x trials)

assert(Ge.number_locations == 1, 'Tuning Response only works for a single location at the moment');

orientations = linspace(0,180,n_orient)';

orient_responses = struct('orientation', mat2cell(orientations, ones(n_orient,1)), ...
    'X', zeros(Ge.dimension_X, S.n_samples, trials), ...
    'G', zeros(Ge.dimension_G, S.n_samples, trials), ...
    'O', zeros(1, S.n_samples, trials), ...
    'L', zeros(1, S.n_samples, trials), ...
    'S', zeros(1, S.n_samples, trials), ...
    'T', zeros(1, S.n_samples, trials));

slices = eye(n_orient);
I = struct('n_zero_signal', 0);

for o_idx=1:n_orient
    o = orientations(o_idx);
    orient_responses(o_idx).orientation = o;
    img = InputImage('1xN', Ge.ny, n_orient, slices(o_idx,:));
    for t=1:trials
        fprintf('O %d/%d\tT %d/%d\n', o_idx, n_orient, t, trials);
        [X, G, O, L, St, T] = Sampling_Gibbs_InPlace_Fast(Ge, S, I, img);
        orient_responses(o_idx).X(:,:,t) = X;
        orient_responses(o_idx).G(:,:,t) = G;
        orient_responses(o_idx).O(:,:,t) = O(1,:);
        orient_responses(o_idx).L(:,:,t) = L(1,:);
        orient_responses(o_idx).S(:,:,t) = St(1,:);
        orient_responses(o_idx).T(:,:,t) = T(1,:);
    end
end

end

