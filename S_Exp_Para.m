function P = S_Exp_Para(mode, varargin)
%S_EXP_PARA create parameters struct to be passed into S_Experiment
%
%   P = S_EXP_PARA(mode, ...), where mode can be any of 'test_2AFC_corr', 
%       'debugging', 'top_downs_disabled', 'test_contrast',  'paper_2AFC_PK',
%       'paper_2AFC_corr', or 'paper_corr_performance'
%
%       Further refinement of parameters can be set using varargin, for example:
%
%          S_EXP_PARA(mode, 'G.number_orientations', 3, ...)
%
%       but NOTE that some parameters depend on others, e.g. 
%
%          P.I.stimulus_contrast = zeros(1,P.G.number_orientations)
%
%       ...and these cannot be set directly via varargin (but feel free to
%       set G.number_orientations, and its dependencies will be re-adjusted
%       accordingly)
%
%   Parameters are grouped into three sub-structs:
%   - P.G sets model parameters
%   - P.S sets the sampler parameters
%   - P.I sets stimulus parameters

%% Begin with the 'default' setup - the most common cases, which will later
% be refined for each mode

% MODEL params
P.G.number_orientations = 2;
P.G.prior_task = [0 1]; % [cardinal, oblique]
P.G.number_locations = 1;
P.G.dimension_X = 256;
P.G.dimension_G = 64;
P.G.kappa_O = [1 0]; % attended and unattendeds
P.G.kappa_G = 3;
P.G.pO = ones(1,P.G.number_orientations)/P.G.number_orientations;
P.G.pL = ones(1,P.G.number_locations)/P.G.number_locations;
P.G.phi_O = [(0:2:2*(P.G.number_orientations-1));...
             (1:2:2*P.G.number_orientations-1)] *...
             pi/2/P.G.number_orientations;
P.G.sigmaStyle = 0.1;
P.G.odds_inc = 1/20; % with 80 samples/sec, every 250ms independent signal
P.G.delta = .08; % strength of X-G coupling for corr & CPs
P.G.nT = 2; % number of possible tasks
P.G.tauStyle = 0.1; % initial value for style
P.G.number_samples_per_evidence = 2; % for dynamic-switching-signal-blocked

% SAMPLER params
P.S.number_repetitions = 512;    % number of trials
P.S.number_burn_in = 0;          % number of burn-in samples
P.S.alpha = 1.0;                 % strength of top-down influence (0 to 1)
P.S.number_samples_to_use = 100; % number of non-burn samples to be used for evidence

% STIMULUS params
P.I.stimulus_regime = 'static';
P.I.n_zero_signal = 20; % number of frames before onset of stimulus
P.I.stimulus_contrast = zeros(1,P.G.number_orientations);

%% Specialize default values based on mode
switch mode
    case 'paper-2AFC-corr'
        P.G.dimension_X = 1024;
        P.G.dimension_G = 256;
        P.S.number_repetitions = 100; 
        
    case 'debugging'
        P.G.number_orientations = 7;
        P.G.dimension_X = 64;
        P.G.dimension_G = 16;
        P.S.number_repetitions = 10;
        P.S.alpha = 0.0;

    case 'top-downs-disabled'                
        P.S.number_repetitions = 2000;
        P.S.alpha = 0.5;        
        
    case 'paper-2AFC-PK'
        P.G.dimension_X = 1024;
        P.G.dimension_G = 256;
        P.G.delta = 80;
        P.G.tauStyle = 1;
        P.G.sigmaStyle = 1e-10;

        P.S.number_repetitions = 16;

        P.I.stimulus_regime='dynamic-switching-signal-blocked';
        P.I.n_zero_signal = 50;
        
    case 'test-2AFC-corr'
        P.G.dimension_X = 128;
        P.G.dimension_G = 16;
        P.G.delta = 0.01;

        P.S.number_repetitions = 128; % number of trials
        
    case 'paper-corr-performance'
        P.G.prior_task=[1 0]; % [cardinal, oblique]
        P.G.dimension_X = 512;        
        
    case 'test-contrast'
        P.G.prior_task=[1 0];
        P.G.delta = 80;
        
        P.S.number_repetitions=1;        
        
    otherwise
        warning('invalid option');
end

% add in any additional args that were specified in 'varargin'
for i=1:2:length(varargin)-1
    param_name = varargin{i};
    param_value = varargin{i+1};
    % strange but effective way of assigning into a nested struct
    % so that we can do S_Exp_Para('debugging', 'S.alpha', 0.9) sort of
    % thing (i.e. reference a struct in a struct)
    field_reference = strsplit(param_name, '.');
    P = subsasgn(P, struct('type', '.', 'subs', field_reference), param_value);
end

% recompute the variables that are dependent on others in case the dependencies changed
% via varargin (note this means varargin cannot be used to override any of these directly)
P.G.pO = ones(1,P.G.number_orientations)/P.G.number_orientations;
P.G.pL = ones(1,P.G.number_locations)/P.G.number_locations;
P.G.phi_O = [[0:2:2*(P.G.number_orientations-1)];[1:2:2*P.G.number_orientations-1]]*pi/2/P.G.number_orientations;
P.I.stimulus_contrast = zeros(1,P.G.number_orientations);

end