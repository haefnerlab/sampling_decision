function [ models_tuning ] = Orientation_Tuning_Like( e, n_orient, trials )
%ORIENTATION_TUNING_LIKE given results of running experiments 'e', which is
%a cell array of results for various model and image paramters, this
%function makes a cell array of the same size as 'e' that contains
%orientation tuning results for each of the corresponding models.

models_tuning = cell(size(e));

for i=1:numel(e)
    model = e{i}
    models_tuning{i} = Orientation_Tuning_Response(model.Projection, model.Sampling, n_orient, trials);
end

end
