function [X, G, O, L, Style, Task, Task_cardinal, Task_oblique] = Sampling_Gibbs_InPlace_Fast(Ge, S, I, Image, DBG)
% SAMPLING_GIBBS_INPLACE_FAST performs gibbs sampling on the given
% generative model with G and X 'layers' of variables
%
%   [X,G,O,L,Style,Task] = SAMPLING_GIBBS_INPLACE_FAST(Ge,S,I,Image,[debug])
%   Uses the generative model Ge, sampling parameters S, stimulus info I,
%   and stimulus Image.
%
%   Return values are:
%   - X: a (neurons x samples) array of spike rates of 'X' layer
%   - G: likewise for the 'G' layer
%   - O: a (orientations+1 x samples) array where the first row gives sampled values and
%        rows 2:orientations+1 give the posterior it was drawn from
%   - L: like 'O' but for locations. (locations+1 x samples),
%        first row is values rest are posterior
%   - Task: like ?O? and ?L? for the currently inferred task
%   - Style: ??? (TODO)

%%%%% By Shizhao Liu, 11/19/2025:
%%%%% Changed the code to accommodate a mode of dual-task-switching, where
%%%%% there are two T variables, T_cardinal and T_oblique. They are
%%%%% independent and can be 0 or 1.

if nargin < 5, DBG = 0; end

%% convenient aliases

% standard deviation of noise at each pixel (iid)
sigy = 1;

nL = Ge.number_locations;
dimX = Ge.dimension_X;
dimG = Ge.dimension_G;
% there is one of each X and G at each location
nX = nL * dimX;
nG = nL * dimG;

%% make input image into column vector(s)
% (one per frame if multiple frames)

nI = size(Image);
switch length(nI)
    case 2
        % [width x height] reshape to [npixels x 1]
        Image = reshape(Image, [numel(Image) 1]);
    case 3
        % [frames x width x height] reshaped to [npixels x frames]
        Image = reshape(permute(Image, [2,3,1]), [nI(2)*nI(3), nI(1)]);
    otherwise
        error(size(Image));
end

%% Kernel initialization
% (kernels are storage for intermediate computations, or lookup tables of
% precomputed probabilities)

% kernel_G(g,x) is the precomputed circular von Mises functions that give
% (part of) the relationship E[x|g]:
% E[x|g] = 1 + delta sum_{k=1}^{n_g} g_k exp(lambda cos 2(phi^x_i-phi^g_k))
%
% hence (in mixed notation), E[x|g] = 1 + sum(kernel_G(g=='on', :))
kernel_G = zeros(dimG, dimX);
norm_kernel_G = 1 / besseli(0, Ge.kappa_G);
for k = 1:dimG
    for l = 1:dimX
        % note that kappa_G is referred to as lambda in the text
        % (TODO be more consistent about naming)
        kernel_G(k,l) = norm_kernel_G * Ge.delta * ...
            exp(Ge.kappa_G * cos(2*(Ge.phi_x(l) - Ge.phi_g(k))));
    end
end

% kernel_O(attended, which_task, which_O, g) is the precomputed circular
% von Mises functions which give the probability that 'g' is 1 given its
% graphical-model-parents (whether it's in the attended location, which
% task is being performed, and the current Orientation belief)

%%%%% By Shizhao Liu 11/19/2025. Added a mode that uses two Ts, one for
%%%%% cardinal, one for oblique. Need to design variable dependency for
%%%%% this mode ("Dual"). "Single" represents Ralf's original mode, where there is only one T. 
switch Ge.switching_mode
    case 'single' 
        % dimension 1 is [attended, unattended], the rest are self-explanatory
        kernel_O = zeros([2 Ge.nT, Ge.number_orientations dimG]);
        % normalization such that sum over g of p(g | O) is 1 (on average one g is
        % 'on' at a time)
        norm_kernel_O = 1 ./ besseli(0, Ge.kappa_O) / dimG;
        for T = 1:Ge.nT
            for O = 1:Ge.number_orientations
                if O == 1 && strcmp(Ge.task, 'detection')
                    delta_cos = 0; % stimulus absent
                else
                    % same for 'discrimination' and 'detection' tasks here
                    delta_cos = cos(2 * (Ge.phi_O(T,O) - Ge.phi_g));
                end
                
                % the difference between attended and unattended is in the width of
                % the von Mises, controlled by kappa_O (see S_Exp_Para)
                kernel_O(1,T,O,:) = norm_kernel_O(1) * exp(Ge.kappa_O(1)*delta_cos);
                kernel_O(2,T,O,:) = norm_kernel_O(2) * exp(Ge.kappa_O(2)*delta_cos);
            end
        end
        kernel_O(kernel_O > 1) = 1;
    case 'dual'
        %%%% Need to design a kernel_O for dual mode
        % Define two separate kernels. One for dependency between O and G
        % in cardinal task, one for dependency in oblique task

        % dimension 1 is [attended, unattended]
        % the rest are self-explanatory
        kernel_O_cardinal = zeros([2, 1, Ge.number_orientations, dimG]);
        kernel_O_oblique  = zeros([2, 1, Ge.number_orientations, dimG]);
        kernel_O_uniform  = ones([2, 1, Ge.number_orientations, dimG]);
         % normalization such that sum over g of p(g | O) is 1 (on average one g is
        % 'on' at a time)
        norm_kernel_O = 1 ./ besseli(0, Ge.kappa_O) / dimG;

        %%%%% define cardinal kernel
        for O = 1:Ge.number_orientations
            delta_cos = cos(2 * (Ge.phi_O(1,O) - Ge.phi_g)); % the first row in phi_O corresponds to the cardinal task
                
            kernel_O_cardinal(1,1,O,:) = norm_kernel_O(1) * exp(Ge.kappa_O(1)*delta_cos);
            kernel_O_cardinal(2,1,O,:) = norm_kernel_O(2) * exp(Ge.kappa_O(2)*delta_cos);
        end

        %%%% define oblique kernel
        for O = 1:Ge.number_orientations
            delta_cos = cos(2 * (Ge.phi_O(2,O) - Ge.phi_g)); % the second row in phi_O corresponds to the oblique task
                
            kernel_O_oblique(1,1,O,:) = norm_kernel_O(1) * exp(Ge.kappa_O(1)*delta_cos);
            kernel_O_oblique(2,1,O,:) = norm_kernel_O(2) * exp(Ge.kappa_O(2)*delta_cos);
        end
        
        %%%% define a uniform kernel (no task is on)
        kernel_O_uniform             = kernel_O_uniform ./ besseli(0, 0) / dimG / 2;


end

%% initialization and allocation of variables

% pO, pL, prior_task are our priors
pO = Ge.pO;
pL = Ge.pL;


X     = zeros(nX, S.n_samples); % X is a flattened [nL, dimX] matrix
G     = zeros(nL, dimG, S.n_samples);
Style = zeros(1, S.n_samples);
% Note: the three _Posterior variables are concatenated (along
% dimension 1) into Task, O, and L after all the sampling is done
% (i.e. in the return values, O(1,:) is the sampled orientation and
% O(2:end,:) is pO_Posterior, likewise for task and L)
L     = zeros(1, S.n_samples);
O     = zeros(1, S.n_samples);
pL_Posterior = zeros(nL, S.n_samples);
pO_Posterior = zeros(Ge.number_orientations, S.n_samples);

switch Ge.switching_mode
    case 'single' 
        prior_task = Ge.prior_task;
        Task  = zeros(1, S.n_samples);
        pT_Posterior = zeros(Ge.nT, S.n_samples);
    case 'dual'
        prior_task_cardinal     = Ge.prior_task_cardinal;
        prior_task_oblique      = Ge.prior_task_oblique;
        Task_cardinal           = zeros(1, S.n_samples);
        Task_oblique            = zeros(1, S.n_samples);
        pT_cardinal_Posterior   = zeros(1, S.n_samples);
        pT_oblique_Posterior    = zeros(1, S.n_samples); 
    
        %%%% create this variable just for the sake of returning
        Task    = nan * ones(1, S.n_samples); 
end



%% Draw initial values from priors
% starting at the top with L, O, Task.. then G.. then X

Style(1) = Ge.tauStyle; % start with prior mean

% Samplel initial values of L, O, and T from their prior
L(1)    = randidx(pL);
O(1)    = randidx(pO);
pL_Posterior(:,1) = pL;
pO_Posterior(:,1) = pO;


switch Ge.switching_mode
    case 'single' 
        Task(1)             = randidx(prior_task);
        pT_Posterior(:,1)   = prior_task;
        % Sample initial G conditioned on T, O, L
        G(:,:,1) = init_G(kernel_O, Task(1), O(1), L(1), nG, nL, dimG);
    case 'dual'
        Task_cardinal(1)             = randidx([1-prior_task_cardinal, prior_task_cardinal]) - 1;
        Task_oblique(1)              = randidx([1-prior_task_oblique, prior_task_oblique]) - 1;
        pT_cardinal_Posterior(:,1)   = prior_task_cardinal;
        pT_oblique_Posterior(:,1)    = prior_task_oblique;

        % Sample initial G conditioned on T, O, L
        %%%% First, get the kernel_O to use based on T variables
        kernel_O_use = get_kernel_O_use_dual(kernel_O_cardinal, kernel_O_oblique, kernel_O_uniform, Task_cardinal(1),Task_oblique(1));
        %%% The second input is task, which should always be 1 because kernel_O_use already accounts for task 
        T_fix = 1;
        
        G(:,:,1) = init_G(kernel_O_use, T_fix, O(1), L(1), nG, nL, dimG);
end


% Sample initial X conditioned on G
X(:,1)   = init_X(kernel_G, G(:,:,1), Ge.delta, nX, dimX);



%% Gibbs sampling

for samp = 2:S.n_samples
    % Copy entire current state forward by 1 step
    X(:,samp) = X(:,samp-1);
    G(:,:,samp) = G(:,:,samp-1);
    O(:,samp) = O(:,samp-1);
    L(:,samp) = L(:,samp-1);
    Style(:,samp) = Style(:,samp-1);
    pO_Posterior(:,samp) = pO_Posterior(:,samp-1);
    pL_Posterior(:,samp) = pL_Posterior(:,samp-1);
    switch Ge.switching_mode
        case 'single'
            Task(:,samp) = Task(:,samp-1);
            pT_Posterior(:,samp) = pT_Posterior(:,samp-1);
        case 'dual'
            Task_cardinal(:,samp)           = Task_cardinal(:,samp-1);
            Task_oblique(:,samp)            = Task_oblique(:, samp-1);
            pT_cardinal_Posterior(:,samp)   = pT_cardinal_Posterior(:,samp-1);
            pT_oblique_Posterior(:,samp)    = pT_oblique_Posterior(:,samp-1);
    end
    

   
    
    % Get Image for this sample
    cur_Image = Image(:, S.access(samp));
    
    % Update priors (accumulated evidence after signal starts)
    if samp > I.n_zero_signal 
        % previous posterior is current prior
        pO = pO_Posterior(:,samp)';
        pL = pL_Posterior(:,samp)';
        if Ge.clamp_prior % do not update prior for task
            switch Ge.switching_mode
                case 'single'
                    prior_task = Ge.prior_task;
                case 'dual'
                    prior_task_cardinal = Ge.prior_task_cardinal;
                    prior_task_oblique  = Ge.prior_task_oblique;
            end
        else
            switch Ge.switching_mode
                case 'single'
                    prior_task = pT_Posterior(:,samp)';
                case 'dual'
                    prior_task_cardinal = pT_cardinal_Posterior(:,samp)';
                    prior_task_oblique  = pT_oblique_Posterior(:,samp)';
            end
        end
    else
        % if no signal, just use initial prior
        pO = Ge.pO;
        pL = Ge.pL;

        switch Ge.switching_mode
            case 'single'
            prior_task = Ge.prior_task;
            case 'dual'
            prior_task_cardinal = Ge.prior_task_cardinal;
            prior_task_oblique  = Ge.prior_task_oblique;
        end

    end
    
    %%%%%%%%%% Update variables
    switch Ge.switching_mode
        case 'single'
            % Update X
            X(:,samp) = sample_X(cur_Image, Ge, kernel_G, X(:,samp), G(:,:,samp), ...
                Style(1,samp), sigy, S.alpha, nX, dimX);
            
            % Update G
            G(:,:,samp) = sample_G(kernel_O, kernel_G, G(:,:,samp), X(:,samp), ...
                Task(1,samp), O(1,samp), L(1,samp), S.alpha, nL, dimG, dimX);
            
            % Update O
            [O(1,samp), pO_Posterior(:,samp)] = ...
                sample_O(kernel_O, Ge, G(:,:,samp), L(1,samp), Task(1,samp), pO, nL);
            if DBG, disp(['i pO: ' num2str([samp pO_Posterior(:,samp)])]); end
            
            % Update L
            [L(1,samp), pL_Posterior(:,samp)] = ...
                sample_L(kernel_O, Ge, G(:,:,samp), Task(1,samp), O(1,samp), pL, nL);
            
            % Update Task
            [Task(1,samp), pT_Posterior(:,samp)] = ...
                sample_Task(kernel_O, Ge, G(:,:,samp), O(1,samp), L(1,samp), prior_task, nL);
            if DBG, disp(['i prior_task: ' num2str([samp pT_Posterior(:,samp)])]); end
            
            % Update Style
            Style(:,samp) = sample_Style(Ge, Style(:,samp), X(:,samp), cur_Image, sigy);
        case 'dual'
            %%%% Instead of creating new sampling functions, use the
            %%%% existing ones, but change input to accommodate
            %%%% dual-task-switching
            %%%% first, get the kernel_O to use for this sample based on
            %%%% Task_cardinal and T_oblique, which are not updated 
            kernel_O_use = get_kernel_O_use_dual(kernel_O_cardinal, kernel_O_oblique, kernel_O_uniform, Task_cardinal(1,samp), Task_oblique(1, samp));
            %%% The task input should always be 1 because kernel_O_use already accounts for task 
            T_fix = 1;
            

            % Update X
            X(:,samp) = sample_X(cur_Image, Ge, kernel_G, X(:,samp), G(:,:,samp), ...
                Style(1,samp), sigy, S.alpha, nX, dimX);
            
            % Update G
            G(:,:,samp) = sample_G(kernel_O_use, kernel_G, G(:,:,samp), X(:,samp), ...
                T_fix, O(1,samp), L(1,samp), S.alpha, nL, dimG, dimX);
            
            % Update O
            [O(1,samp), pO_Posterior(:,samp)] = ...
                sample_O(kernel_O_use, Ge, G(:,:,samp), L(1,samp), T_fix, pO, nL);
            if DBG, disp(['i pO: ' num2str([samp pO_Posterior(:,samp)])]); end
            
            % Update L
            [L(1,samp), pL_Posterior(:,samp)] = ...
                sample_L(kernel_O_use, Ge, G(:,:,samp), T_fix, O(1,samp), pL, nL);
            
            % Update Task_cardinal
            [Task_cardinal(1,samp), pT_cardinal_Posterior(:,samp)] = ...
                sample_task_dual(kernel_O_cardinal,kernel_O_oblique,kernel_O_uniform, Ge, G(:,:,samp), O(1,samp), L(1,samp), Task_oblique(1,samp), prior_task_cardinal, nL, 'cardinal');
            
            % Update Task_oblique
            [Task_oblique(1,samp), pT_oblique_Posterior(:,samp)] = ...
                sample_task_dual(kernel_O_cardinal,kernel_O_oblique,kernel_O_uniform, Ge, G(:,:,samp), O(1,samp), L(1,samp), Task_cardinal(1,samp), prior_task_oblique, nL, 'oblique');
            
            if DBG, disp(['i prior_task_cardinal: ' num2str([samp pT_cardinal_Posterior(:,samp)])]); end
            if DBG, disp(['i prior_task_oblique: ' num2str([samp pT_oblique_Posterior(:,samp)])]); end
            
            % Update Style
            Style(:,samp) = sample_Style(Ge, Style(:,samp), X(:,samp), cur_Image, sigy);
    end
            
end

if DBG
    Get_Figure('Post-debug'); plot(pO_Posterior(2,:),'*-');
end
if DBG
    Get_Figure('1'); Subplot(2,1,2,1);
    Histogram(TAU);
    Subplot(2); Histogram(MUK);
end

if DBG, disp(['size O: ' num2str(size(O))]); end

O   (2:1+Ge.number_orientations,:) = pO_Posterior;
L   (2:1+nL,:) = pL_Posterior;
switch Ge.switching_mode
    case 'single'
        Task(2:1+Ge.nT,:) = pT_Posterior;
    case 'dual'
        Task_cardinal(2,:) = pT_cardinal_Posterior;
        Task_oblique(2,:)  = pT_oblique_Posterior; 
end
end

function newX = sample_X(img, Ge, kernel_G, X, G, s, sigy, alpha, nX, dimX)
%SAMPLE_X do a Gibbs sampling step on a random order of all X variables

newX = X;
% randomize update order for X
order = randperm(nX);
% X_k to be sampled
for i = order
    % x is drawn from a 'cut gaussian' with mean mu and variance sig^2
    idx_no_i = [1:i-1 i+1:nX]; % all but i index
    sig = sigy / s;
    kO  = mod(i-1, dimX)+1;  % index of X_i's Gabor orientation
    kL  = 1 + (i-kO) / dimX; % index of X_i's location
    % 0<=alpha<=1 controls the strength of top-down influence (of G on X)
    tau = 1 + alpha * G(kL,:) * kernel_G(:,kO); % see def of kernel_G
    if tau <= 0, error('in sampling x_i, got tau_i <= 0'); end
    mu = s * Ge.G(:,i)' * img + ...
        s^2 * Ge.R(i,idx_no_i) * newX(idx_no_i) - ...
        sigy^2 / tau;
    mu = mu / s^2;
    newX(i) = Cut_Gaussian('random', mu, sig);
    if ~isreal(newX(i)), error(['X NaN, mu sigma: ' num2str([mu sig])]); end
end
end

function newG = sample_G(kernel_O, kernel_G, G, X, T, O, L, alpha, nL, dimG, dimX)
%SAMPLE_G do a Gibbs sampling step on a random order of all g variables

newG = G;

% unflatten X for easy slicing later
% (X is stored with the orientations dimension contiguous)
X = reshape(X, [dimX, nL]);

% randomize update order for G
order = randperm(nL * dimG);
for jk = order, [l, k] = ind2sub([nL dimG], jk);
    if L == l, attn = 1; else attn = 2; end
    % p(g|...) factorizes to p(g|D,T,O)p(x|g)
    % p_g_topdown gets us p(g|D,T,O) with 'alpha controlling the strength
    % of this top-down influence
    %
    % note that (1-alpha)*(sum over all orientations except O) is
    % implemented as [(1-alpha)*(sum over all) - (1-alpha)(just O)]
    p_g_topdown = kernel_O(attn, T, O, k) + ...
        (1.0 - alpha) * sum(kernel_O(attn, T, :, k)) - ...
        (1.0 - alpha) * kernel_O(attn, T, O, k);
    
    % now we get the 'bottom up' effect of x on g. To do this, we consider
    % both cases, G(j,k) = 0 or 1, in terms of their relative likelihood
    % of producing X
    log_p_each_case = zeros(1,2); % log probability when G(j,k) is 0,1 respectively
    for Gjk_val = [0 1], case_idx = Gjk_val+1;
        
        % set G to this case to see what likelihoods we get
        newG(l,k) = Gjk_val;
        
        % if G(g,k) is 'off', the top-down likelihood is 1-p_g_topdown
        switch Gjk_val
            case 0, log_p_each_case(1) = log(1-p_g_topdown);
            case 1, log_p_each_case(2) = log(p_g_topdown);
        end
        
        % 'tau' is the expected value of x given g, defined mathematically
        %   E[x|g] = 1 + delta sum_{k=1}^{n_g} g_k exp(lambda cos 2(phi^x_i-phi^g_k))
        % where kernel_G(k,i) contains
        %   delta * exp(lambda cos 2(phi^x_i-phi^g_k))
        %
        % so tau = E[x|g] is just 1 + (dot product of current G with kernel)
        tau = 1 + (newG(l,:) * kernel_G)';
        
        % x|g is drawn from an exponential distribution with mean 'tau', so
        % its log-likelihood is -log(tau) - x/tau
        log_p_each_case(case_idx) = ...
            log_p_each_case(case_idx) + ...
            sum(-log(tau) - X(:,l)./tau);
    end
    
    % return from log-probability to probability space in a relatively
    % numerically stable way (avoiding exp(very_very_negative_number))
    pG = exp(log_p_each_case-max(log_p_each_case));
    pG = pG / sum(pG);
    
    %ERROR HANDLING
    if any(~isreal(pG))
        disp(['tau : ' num2str(tau')]);
        disp(['aux : ' num2str(log_p_each_case)]);
        error('NaN probability p(G|...)');
    end
    
    if pG(2)>1, error('p(g=1|...) > 1'); end
    
    % new weighted coin flip for G (pG(2) contains prob G is 1)
    newG(l,k) = binornd(1, pG(2), 1);
end
end

function [newO, pO] = sample_O(kernel_O, Ge, G, L, T, pO, nL)
%UPDATE_O compute new posterior over O and take a sample
%
% P(O|G,..) \propto P(G|O,..)P(O) = pO * product_k P(g_k|O,..)
%
% where kernel_O has much of the joint P(attending, T, O, G) precomputed

% TODO - O, T, and L are interdependent - do we really want to stop updating when pO reaches 1? Will it ever reach 1?
if max(pO) < 1 % there is some uncertainty left about O
    log_like_O = zeros(size(pO));
    for l = 1:nL
        if L == l, attn = 1; else attn = 2; end
        
        % log(P(O|G)) sum over G dimension of log kernels
        log_like_O = log_like_O ...
            +reshape(sum(log(1-kernel_O(attn, T, :, G(l,:)==0)),4), size(pO)) ...
            +reshape(sum(log(  kernel_O(attn, T, :, G(l,:)==1)),4), size(pO));
        
        % take a step moving log_pO towards log_like_O
        log_pO = Ge.odds_inc * log_like_O + log(pO);
        % return from log-probability to probability space in a relatively
        % numerically stable way (avoiding exp(very_very_negative_number))
        pO = exp(log_pO-max(log_pO));
        pO = pO / sum(pO);
    end
end
if any(isnan(pO)), error(['NaN value in posterior over O.. ' num2str(pO)]); end
% sample guessed orientation: a multinomial random draw
newO = randidx(pO);
end

function [newL, pL] = sample_L(kernel_O, Ge, G, T, O, pL, nL)
%UPDATE_L compute new posterior over L and take a sample

if max(pL) < 1 % there is some uncertainty left about L
    log_like_L = zeros(size(pL));
    % compute posterior at each location L
    for L = 1:nL
        % requires a nested loop because in order to check each candidate
        % value for L we must see how it affects likelhiood G at all other
        % locations
        for l = 1:nL
            if L == l, attn = 1; else attn = 2; end
            log_like_L(L) = log_like_L(L)...
                +reshape(sum(log(1-kernel_O(attn, T, O, G(l,:) == 0)),4), size(pL))...
                +reshape(sum(log(  kernel_O(attn, T, O, G(l,:) == 1)),4), size(pL));
        end
    end
    
    % take a step moving log_pL towards log_like_L
    log_pL = Ge.odds_inc * log_like_L + log(pL);
    % return from log-probability to probability space in a relatively
    % numerically stable way (avoiding exp(very_very_negative_number))
    pL = exp(log_pL - max(log_pL));
    pL = pL / sum(pL);
end
if any(isnan(pL)), error(['NaN value in posterior over L.. ' num2str(pL)]); end
% sample a guessed location according to weighted posterior
newL = randidx(pL);
end

function [newT, pT] = sample_Task(kernel_O, Ge, G, O, L, pT, nL)
%SAMPLE_TASK compute new posterior over Task and take a sample

if max(pT) < 1 % there is some uncertainty left about Task
    log_like_T = zeros(size(pT));
    
    for l = 1:nL
        % L is spatial location of attention, so attn state is 1 ('attended') when L==l, and 2
        % ('unattended') otherwise
        if L == l, attn = 1; else attn = 2; end
        % kernel_O has shape (2 (attention states), #tasks, #orientations-per-task, #G). Compute
        % p(Task|G) by summing log-evidence from all G==0 and G==1 states. kernel_O is a lookup
        % table of these Task|G conditional probabilities
        log_like_T = log_like_T...
            +reshape(sum(log(  kernel_O(attn,:,O,G(l,:) == 1)),4), size(pT))...
            +reshape(sum(log(1-kernel_O(attn,:,O,G(l,:) == 0)),4), size(pT));
    end
    
    % take a step moving log_pT towards log_like_T
    log_prior_task = Ge.odds_inc * log_like_T + log(pT);
    % return from log-probability to probability space in a relatively
    % numerically stable way (avoiding exp(very_very_negative_number))
    pT = exp(log_prior_task-max(log_prior_task));
    pT = pT/sum(pT);
end
if any(isnan(pT)), error(['NaN value in posterior over T.. ' num2str(pT)]); end
newT = randidx(pT);
end

function newStyle = sample_Style(Ge, Style, X, img, sigy)
% SAMPLE_STYLE sample a new value for the 'style' or 's' variable

xRx = X' * Ge.R * X;
if xRx > 0, error('xRx must be <0!'); end

yGx = img' * Ge.G * X;
% TODO - avoid switching every loop iteration. Not like prior_style is changing..
switch Ge.prior_style
    case 'sparse'
        sigs = sigy / sqrt(-xRx);
        mus = (yGx/sigy^2 - 1/Ge.tauStyle) * sigy^2/(-xRx);
    case 'slow'
        sigs = sqrt(1 / (1/Ge.sigmaStyle^2 - xRx/sigy^2));
        mus = sigs^2 * (yGx/sigy^2 + Style/Ge.sigmaStyle^2);
end
newStyle = Cut_Gaussian('random', mus, sigs);
end

function newG = init_G(kernel_O, Task, O, L, nG, nL, dimG)
%INIT_G draw initial values for G based on inital values of O and L

newG = zeros(nL, dimG);
% order doesn't matter for initialization
for jk = 1:nG
    [l, k] = ind2sub([nL, dimG], jk);
    % if the location matches, use 'attended' kernel (slice 1) otherwise
    % the 'unattended' kernel (slice 2)
    if L == l, attn = 1; else attn = 2; end
    prob = kernel_O(attn, Task, O, k);
    newG(l,k) = binornd(1,prob,1); % 0 or 1
end
end

function newX = init_X(kernel_G, G, delta, nX, dimX)
%INIT_X draw initial values for X based on inital values of G

% order doesn't matter for initialization
for i = nX:-1:1
    iO = mod(i-1, dimX) + 1; % unflatten index (col, orientation)
    iL = 1 + (i-iO) / dimX;  % unflatten index (row, location)
    % tau(i) = 1 + G(kL,:) * kernel_G(:, kO); % OLD version
    tau(i) = delta + G(iL,:) * kernel_G(:, iO);
end
newX = exprnd(tau);
end

function idx = randidx(probs)
probs = probs(:) / sum(probs);
idx = find(rand < cumsum(probs), 1);
end



function [newT,pT] = sample_task_dual(kernel_O_cardinal, kernel_O_oblique, kernel_O_uniform, Ge, G, O, L, T_other, pT, nL, option)

%%%% sample one task variable T, conditioned on the other task variable
%%%% use option to specific which task T to sample

if pT < 1 & pT > 0 % there is some uncertainty left about Task
    log_like_T = zeros(1,2);
    
    for l = 1:nL
        % L is spatial location of attention, so attn state is 1 ('attended') when L==l, and 2
        % ('unattended') otherwise
        if L == l, attn = 1; else attn = 2; end
        % kernel_O has shape (2 (attention states), #tasks, #orientations-per-task, #G). Compute
        % p(Task|G) by summing log-evidence from all G==0 and G==1 states. kernel_O is a lookup
        % table of these Task|G conditional probabilities

        switch option
            case 'cardinal'
                %%%%% Compute the kernel_O when T_cardinal  = 0 or 1 (on or off)
                %%%% T_other is T_oblique.
                
                kernel_O_use_off = get_kernel_O_use_dual(kernel_O_cardinal, kernel_O_oblique, kernel_O_uniform, 0, T_other);
                kernel_O_use_on  = get_kernel_O_use_dual(kernel_O_cardinal, kernel_O_oblique, kernel_O_uniform, 1, T_other);
            case 'oblique'
                %%%%% Compute the kernel_O when T_oblique  = 0 or 1 (on or off)
                %%%% T_other is T_cardinal.
                
                kernel_O_use_off = get_kernel_O_use_dual(kernel_O_cardinal, kernel_O_oblique, kernel_O_uniform, T_other, 0);
                kernel_O_use_on  = get_kernel_O_use_dual(kernel_O_cardinal, kernel_O_oblique, kernel_O_uniform, T_other, 1);

        end
        
        %%%% likelihood of this task is off
        log_like_T(1) = log_like_T(1)...
            +reshape(sum(log(  kernel_O_use_off(attn,1,O,G(l,:) == 1)),4), 1,1)...
            +reshape(sum(log(1-kernel_O_use_off(attn,1,O,G(l,:) == 0)),4), 1,1);

        %%%% likelihood of this task is on
        log_like_T(2) = log_like_T(2)...
            +reshape(sum(log(  kernel_O_use_on(attn,1,O,G(l,:) == 1)),4), 1,1)...
            +reshape(sum(log(1-kernel_O_use_on(attn,1,O,G(l,:) == 0)),4), 1,1);

       
    end
    
    % take a step moving log_pT towards log_like_T
    % pT is a scalar, the prior prob that T is on
    log_prior_task = Ge.odds_inc * log_like_T + [log(1-pT), log(pT)] ;
    % return from log-probability to probability space in a relatively
    % numerically stable way (avoiding exp(very_very_negative_number))
    posT    = exp(log_prior_task-max(log_prior_task));
    posT    = posT/sum(posT);
    pT      = posT(2); % pT is a scalar that this task is on. Need to return pT
end
if any(isnan(pT)), error(['NaN value in posterior over T.. ' num2str(pT)]); end

newT = randidx([1-pT, pT]) - 1; %%%% sample T. [p_task_off, p_task_on]. randidx returns 1 or 2, need to minus 1
end



function  kernel_O_use = get_kernel_O_use_dual(kernel_O_cardinal, kernel_O_oblique,kernel_O_uniform, T_cardinal, T_oblique)

%%%%% By Shizhao Liu 11/20/2025. Get the kernel_O to use under
%%%%% dual-task-switching mode, based on pre-defined kernel_O_cardinal,
%%%%% kernel_O_oblique, and the samples of two task variable T.

%%%%% This line of code is equivalent to:
% % if T_cardinal == 1 & T_oblique == 0
% %     % Just like cardinal
% %     kernel_O_use = kernel_O_cardinal + kernel_O_uniform;
% % elseif T_cardinal == 0 & T_oblique == 1
% %     % Just like oblique
% %     kernel_O_use = kernel_O_oblique + kernel_O_uniform;
% % elseif T_cardinal == 0 & T_oblique == 0
% %     % Uniform profile
% %     kernel_O = kernel_O_uniform + kernel_O_uniform;
% %     
% % elseif T_cardinal == 1 & T_oblique == 1
% %     % Sum of two kernels
% %     kernel_O_use = kernel_O_cardinal + kernel_O_oblique;
% % end
kernel_O_use = kernel_O_cardinal * T_cardinal + kernel_O_oblique * T_oblique + ...
            (1 - T_cardinal) * (1 - T_oblique) * kernel_O_uniform;
end