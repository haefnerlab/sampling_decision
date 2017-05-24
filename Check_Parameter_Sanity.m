function Check_Parameter_Sanity(P)
%CHECK_PARAMETER_SANITY warn or error certain cases of parameter values

if ~any(strcmp(P.I.stimulus_regime, {'blank', 'static', 'static-delayed', 'dynamic-delayed', 'dynamic-switching-signal','dynamic-switching-signal-blocked', 'dynamic-shuffled-signal', 'dynamic-shuffled-signal-blocked'}))
    error('Unknown stimulus regime: %s', P.I.stimulus_regime);
end

if ~any(strcmp(P.I.stimulus_regime, {'blank', 'static'})) ...
        && P.I.n_zero_signal >= P.I.n_frames
    error('n_zero_signal (%d) >= n_frames (%d). Nothing left to be signal!', P.I.n_zero_signal, P.I.n_frames);
end

if any(strcmp(P.I.stimulus_regime, {'dynamic-delayed', 'dynamic-switching-signal','dynamic-switching-signal-blocked','dynamic-shuffled-signal','dynamic-shuffled-signal-blocked'})) ...
        && P.G.number_locations > 1,
    error('dynamic signal at multiple locations not implemented!');
end

if strcmp(P.G.task, 'detection') && P.G.number_orientations ~= 2
    warning('detection task expects number_orientations to be 2');
end

if P.G.kappa_O(2) >= P.G.kappa_O(1)
    warning('are you sure about P.kappa?'); 
end

if max(abs(P.I.stimulus_contrast)) > 0 && P.G.prior_task(1) < P.G.prior_task(2)
    warning('non-zero stim contrasts imply Task = 1, not 2!!');
end

if sum(abs(diag(P.G.R)+1) > 1e-5)
    error('R_ii has to be -1!');
end

if ~any(strcmp(P.G.prior_style, {'sparse', 'slow'}))
    error('unknown style prior: %s', P.G.prior_style);
end

if ~any(strcmp(P.G.task, {'discrimination', 'detection'}))
    error('unknown task: %s', P.G.task);
end

end