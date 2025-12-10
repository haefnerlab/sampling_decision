function out = S_Experiment(params)
%S_EXPERIMENT runs the sampling decision model
%
%   out = S_Experiment(params) where params has been created by a call to
%   S_Exp_Para, runs the specified experiment and returns information in
%   the struct 'out'.

P = params;
Check_Parameter_Sanity(P);

%% local copies of variables with shorter or more descriptive names
n_trials = P.S.number_repetitions;
regime = P.I.stimulus_regime;
contrast = P.I.stimulus_contrast;
n_frames = P.I.n_frames;
im_type = P.I.fct;
n_locs = P.G.number_locations;
im_height = P.G.ny;
n_zero_sig = P.I.n_zero_signal;
n_neurons = P.G.dimension_X * P.G.number_locations;
n_pixels = size(P.G.G,1);
%%%%% By Shizhao Liu 03/05/2025. Adding a parameter (image_task) that controls whether
%%%%% the input images are cardinal (0/90 degrees, default) or oblique (45/135)
%%%%% This parameter is also passed to create_trial_stimulus and then
%%%%% Input_image
image_task = P.I.image_task;
% function handle to make a new stimulus
create_stimulus_handle = @() create_trial_stimulus(regime, im_type, contrast, im_height, n_locs, n_frames, n_zero_sig, image_task);

% pre-allocate the return variables
Signal = zeros(n_trials, n_neurons, n_frames);

%%%%% By Shizhao Liu 12/04/2025. Pre-allocate orientation signal
run_ori_energy = P.I.run_ori_energy;

if run_ori_energy
    n_ori_bin       = P.I.n_ori_bin;
    n_frame_energy  = n_frames - n_zero_sig;
    orientation_energy = zeros(n_trials, n_ori_bin,n_ori_bin);
end
%% loop over trials, parallelized over multiple cores if possible
% (if not, parfor defaults to a backwards for loop)
for i = 1:n_trials
    
    if mod(i,20) == 0
        disp(['Computing Repetition ' num2str(i) ' / ' num2str(n_trials)]);
    end
    
    Y = create_stimulus_handle();
    
    % Perform a trial
    % (suppressing warning that P is broadcast.. it's unavoidable)
    switch P.G.switching_mode
        case 'single' 
            [aux_X{i}, aux_G{i}, aux_O{i}, aux_L{i}, aux_S{i}, aux_T{i}] = ...
                Sampling_Gibbs_InPlace_Fast(P.G, P.S, P.I ,Y); %#ok<PFBNS>
        case 'dual'
             [aux_X{i}, aux_G{i}, aux_O{i}, aux_L{i}, aux_S{i}, aux_T{i}, aux_T_cardinal{i}, aux_T_oblique{i}] = ...
                Sampling_Gibbs_InPlace_Fast(P.G, P.S, P.I ,Y); %#ok<PFBNS>
    end
    
    % Signal at trial i, neuron j, frame k is mean convolution of the image
    % with the neuron's projective field
    switch regime
        case {'static', 'blank'}
            Signal(i,:) = (Y(:)' * P.G.G) / n_pixels;
        case {'static-delayed','dynamic-delayed','dynamic-switching-signal','dynamic-switching-signal-blocked'}
            % do convolution at each frame (tmp variable to assist in
            % parfor slicing)
            tmp = zeros(n_neurons, n_frames);
            for k = 1:n_frames
                tmp(:, k) = (Y(k,:) * P.G.G) / n_pixels;
            end
            Signal(i,:) = tmp(:)';
        otherwise, error(regime);
    end

    %%% also save orientation energy of images
    if run_ori_energy
        
        for im = 1:n_frames
            [orientation_bins_center, orientation_energy(i,:,im)] = get_ori_energy_model(squeeze(Y(im,:,:)), n_ori_bin);
        end
    end
    
end

%% Copy results to output struct
out.Signal = Signal;
if run_ori_energy
    out.oriEnergy = orientation_energy;
    out.orientation_bins_center  = orientation_bins_center;
end
for i = n_trials:-1:1
    out.X(i,:,:)   = aux_X{i};
    out.G(i,:,:,:) = aux_G{i};
    out.O(i,:,:)   = aux_O{i};
    out.L(i,:,:)   = aux_L{i};
    out.T(i,:,:)   = aux_T{i};
    if strcmp(P.G.switching_mode,'dual')
        out.T_cardinal(i,:,:)   = aux_T_cardinal{i};
        out.T_oblique(i,:,:)    = aux_T_oblique{i};
    end
    out.S(i,:)     = aux_S{i};
    
end

% Rename/alias variables for backwards compatibility with old analysis code
out.Projection = P;
out.InputImage = P.I;
out.Sampling = P.S;
out = BackwardsComp(out);

end

function stim = create_trial_stimulus(regime, im_type, contrast, im_height, n_locs, n_frames, n_zero_sig, image_task)
% helper function to create the stimulus for each trial
%%%%% By Shizhao Liu 03/05/2025. Adding a parameter (image_task) that controls whether
%%%%% the input images are cardinal (0/90 degrees, default) or oblique (45/135)
%%%%% This parameter is then passed to Input_image.

if strcmp(regime, 'blank')
    stim = zeros(im_height, n_locs*im_height);
    return
end
% As long as regime is not 'blank', the image itself changes, drawn
% from some distribution based on the 'signal', i.e. contrast at
% each location, which varies sample to sample in the
% dynamic-switching-* regimes
switch regime
    case 'static'
        signal = contrast;
    case 'static-delayed'
        % signal starts as zeros then becomes the stimulus
        signal(1,:,:) = zeros(size(contrast));
        signal(2,:,:) = contrast;
    case 'dynamic-delayed'
        % starting at n_zero_signal+1, signal starts being dynamic
        signal = zeros(n_frames,1,2);
        for frame = n_zero_sig+1:n_frames
            signal(frame,1,:) = contrast;
        end
    case {'dynamic-switching-signal','dynamic-switching-signal-blocked'}
        signal = zeros(n_frames,1,2);
        for frame = n_zero_sig+1:n_frames
            on = 1+binornd(1,0.5);
            signal(frame,1,on) = contrast(on);
        end
end
%stim = InputImage(im_type, n_locs, im_height, signal);
stim = InputImage(im_type, n_locs, im_height, signal, image_task);
end

function out = BackwardsComp(out)
%Manually handle variables whose names have changed
out.Projection.nO = out.Projection.G.number_orientations;
out.Projection.nL = out.Projection.G.number_locations;
out.Projection.nX = out.Projection.G.dimension_X;
out.Projection.nG = out.Projection.G.dimension_G;
out.Projection.pT = out.Projection.G.prior_task;

if isfield(out.Projection.G, 'priot_task_cardinal')
    out.Projection.pT_cardinal = out.Projection.G.prior_task_cardinal;
    out.Projection.pT_oblique  = out.Projection.G.prior_task_oblique;
end
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
for i = 1:l
    out.Projection.(names{i}) = data{i};
end
names = fieldnames(out.Projection.S);
data = struct2cell(out.Projection.S);
l = length(names);
for i = 1:l
    out.Projection.(names{i}) = data{i};
end
names = fieldnames(out.Projection.I);
data = struct2cell(out.Projection.I);
l = length(names);
for i = 1:l
    out.Projection.(names{i}) = data{i};
end
end