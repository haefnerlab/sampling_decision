function out = S_Experiment(params, varargin)
%S_EXPERIMENT runs the sampling decision model
%
%   out = S_Experiment(params) where params has been created by a call to
%   S_Exp_Para, runs the specified experiment and returns information in
%   the struct 'out'.

P = params;

% 2x2: two locations, two orientations (NIPS)
% 1xN: 1 location, many orientations
% nx2: n locations, two orientations
% nxN: n locations, many orientations (Current)

% Generative model
projective_fields = C_Projection('nxN', P.G.nx, P.G.dimension_X, P.G.dimension_G, P.G.number_locations);
P.fct = 't-l-op-g-s';
P.G.G     = projective_fields.G;
P.G.phi_x = projective_fields.phi_x;
P.G.phi_g = projective_fields.phi_g;
P.G.ny    = projective_fields.ny;
P.I.x     = projective_fields.x;
P.I.y     = projective_fields.y;

% Repeating Gibbs sampling P.number_repetitions times
out.X = zeros(P.S.number_repetitions,P.G.number_locations*P.G.dimension_X, P.S.n_samples);

% the below two lines are my way of dealing with Matlab's lack of
% macros. use comments to switch between serial and parallel processing

%PARALLEL SUPPORT
parfor i = 1:P.S.number_repetitions,
    %warning('serial!!'); for i = 1:P.S.number_repetitions, ProgressReport(10,i-1,P.S.number_repetitions);
    
    if mod(i,20) == 0
        disp(['Computing Repetition ' num2str(i) ' / ' num2str(P.S.number_repetitions)]);
    end
    switch P.I.stimulus_regime
        case 'static'
            signal = P.I.stimulus_contrast;
        case 'static-delayed'
            signal(1,:,:) = zeros(size(P.I.stimulus_contrast));
            signal(2,:,:) = P.I.stimulus_contrast;
        case 'dynamic-delayed'
            if P.G.number_locations>1, error('not implemented'); end
            signal = zeros(P.I.n_frames,1,2);
            for k = P.I.n_zero_signal+1:P.I.n_frames
                signal(k,1,:) = P.I.stimulus_contrast;
            end
        case {'dynamic-switching-signal','dynamic-switching-signal-blocked'}
            if P.G.number_locations>1, error('not implemented'); end
            signal = zeros(P.I.n_frames,P.G.number_locations,2);
            for k = P.I.n_zero_signal+1:P.I.n_frames
                on = 1+binornd(1,0.5);
                signal(k,1,on) = P.I.stimulus_contrast(on);
            end
    end
    switch P.I.stimulus_regime
        case 'blank'
            Y = zeros(P.G.ny,P.G.number_locations*P.G.ny);
        otherwise
            Y = InputImage(P.I.fct,P.G.number_locations,P.G.ny,signal);
    end
    
    %Perform a trial
    [aux_X{i} aux_G{i} aux_O{i} aux_L{i} aux_S{i} aux_T{i}] = ...
        Sampling_Gibbs_InPlace_Fast(P.G, P.S, P.I ,Y);
    
    for j = 1:P.G.dimension_X*P.G.number_locations
        switch P.I.stimulus_regime
            case {'static','blank'}
                Signal{i}{j} = mean(mean(Y.*reshape(P.G.G(:,j),P.G.ny,P.G.nx)));
            case {'static-delayed','dynamic-delayed','dynamic-switching-signal','dynamic-switching-signal-blocked'}
                for k = 1:P.I.n_frames
                    Signal{i}{j,k} = mean(mean(squeeze(Y(k,:,:)).*...
                        reshape(P.G.G(:,j),P.G.ny,P.G.nx)));
                end
            otherwise, error(P.I.stimulus_regime);
        end
    end
end
out.Signal = zeros(P.S.number_repetitions, P.G.dimension_X*P.G.number_locations, P.I.n_frames);
for i = P.S.number_repetitions:-1:1
    %if SAVE_Y, out.Y(i,:,:,:) = Y{i}; end
    out.X(i,:,:)  =aux_X{i};
    out.G(i,:,:,:) = aux_G{i};
    out.O(i,:,:)  =aux_O{i};
    out.L(i,:,:)  =aux_L{i};
    out.T(i,:,:)  =aux_T{i};
    out.S(i,:)    =aux_S{i};
    for j = P.G.dimension_X*P.G.number_locations:-1:1
        switch P.I.stimulus_regime
            case {'static','blank'}
                out.Signal(i,j) = Signal{i}{j};
            case {'static-delayed','dynamic-delayed','dynamic-switching-signal','dynamic-switching-signal-blocked'}
                for k = P.I.n_frames:-1:1
                    out.Signal(i,j,k) = Signal{i}{j,k};
                end
            otherwise, error(P.I.stimulus_regime);
        end
    end
end

%Convert back to single variable P for output, thne perform backwards
%compatibility handling


out.Projection = P; out.InputImage = P.I; out.Sampling = P.S;

out = BackwardsComp(out);

end

function out = BackwardsComp(out)
%Manually handle variables whose names have changed
out.Projection.nO = out.Projection.G.number_orientations;
out.Projection.nL = out.Projection.G.number_locations;
out.Projection.nX = out.Projection.G.dimension_X;
out.Projection.nG = out.Projection.G.dimension_G;
out.Projection.pT = out.Projection.G.prior_task;
out.Projection.phi_x = out.Projection.G.phi_x;
out.Projection.phi_g = out.Projection.G.phi_g;
out.InputImage.dyn = out.Projection.I.stimulus_regime;
out.InputImage.c = out.Projection.I.stimulus_contrast;
out.Sampling.n_burn = out.Projection.S.number_burn_in;
out.Sampling.n_use = out.Projection.S.number_samples_to_use;
out.Sampling.nsamples_per_evidence = out.Projection.G.number_samples_per_evidence;
out.Projection.delta = out.Projection.G.delta;
out.Projection.kappa_O = out.Projection.G.kappa_O;
out.Projection.kappa_G = out.Projection.G.kappa_G;


%Automatically update non-name changing variables
names = fieldnames(out.Projection.G);
data = struct2cell(out.Projection.G);
l = length(names);
for i = 1:l;
    out.Projection.(names{i}) = data{i};
end
names = fieldnames(out.Projection.S);
data = struct2cell(out.Projection.S);
l = length(names);
for i = 1:l;
    out.Projection.(names{i}) = data{i};
end
names = fieldnames(out.Projection.I);
data = struct2cell(out.Projection.I);
l = length(names);
for i = 1:l;
    out.Projection.(names{i}) = data{i};
end
end