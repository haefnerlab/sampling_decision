function [X, G, O, L, Style, Task] = Sampling_Gibbs_InPlace_Fast(Ge, S, I, Image, debug)
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


if nargin < 5, debug = 0; end

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
    for i = 1:dimX
        % note that kappa_G is referred to as lambda in the text
        % (TODO be more consistent about naming)
        kernel_G(k,i) = norm_kernel_G * Ge.delta * ...
            exp(Ge.kappa_G * cos(2*(Ge.phi_x(i) - Ge.phi_g(k))));
    end
end

% kernel_O(attended, which_task, which_O, g) is the precomputed circular
% von Mises functions which give the probability that 'g' is 1 given its
% graphical-model-parents (whether it's in the attended location, which
% task is being performed, and the current Orientation belief)

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

%% initialization and allocation of variables

% pO, pL, prior_task are our priors
pO = Ge.pO;
pL = Ge.pL;
prior_task = Ge.prior_task;

X     = zeros(nX, S.n_samples); % X is a flattened [nL, dimX] matrix
G     = zeros(nL, dimG, S.n_samples);
Style = zeros(1, S.n_samples);
% Note: the three _Posterior variables are concatenated (along
% dimension 1) into Task, O, and L after all the sampling is done
% (i.e. in the return values, O(1,:) is the sampled orientation and
% O(2:end,:) is pO_Posterior, likewise for task and L)
L     = zeros(1, S.n_samples);
O     = zeros(1, S.n_samples);
Task  = zeros(1, S.n_samples);
pL_Posterior = zeros(nL, S.n_samples);
pO_Posterior = zeros(Ge.number_orientations, S.n_samples);
pT_Posterior = zeros(Ge.nT, S.n_samples);

%% Draw initial values from priors
% starting at the top with L, O, Task.. then G.. then X

Style(1) = Ge.tauStyle; % start with prior mean
% find() of cumsum()>rand is the discrete version of inverting the cdf to
% draw a sample from the prior
L(1)    = find(cumsum(pL) > rand(1), 1, 'first');
O(1)    = find(mnrnd(1,pO) == 1, 1);
Task(1) = find(mnrnd(1,prior_task) == 1, 1);
pL_Posterior(:,1) = pL;
pO_Posterior(:,1) = pO;
pT_Posterior(:,1) = prior_task;

G(:,:,1) = init_G(kernel_O, Task(1), O(1), L(1), nG, nL, dimG);
X(:,1)   = init_X(kernel_G, G(:,:,1), Ge.delta, nX, dimX);

%% Gibbs sampling

for samp = 2:S.n_samples
    % Copy entire current state forward by 1 step
    X(:,samp) = X(:,samp-1);
    G(:,:,samp) = G(:,:,samp-1);
    O(:,samp) = O(:,samp-1);
    L(:,samp) = L(:,samp-1);
    Task(:,samp) = Task(:,samp-1);
    Style(:,samp) = Style(:,samp-1);
    pO_Posterior(:,samp) = pO_Posterior(:,samp-1);
    pL_Posterior(:,samp) = pL_Posterior(:,samp-1);
    pT_Posterior(:,samp) = pT_Posterior(:,samp-1);
    
    % Get Image for this sample
    cur_Image = Image(:, S.access(samp));
    
    % Update X
    X(:,samp) = sample_X(cur_Image, Ge, kernel_G, X(:,samp), G(:,:,samp), ...
        Style(1,samp), sigy, S.alpha, nX, dimX);
    
    % Update G
    G(:,:,samp) = sample_G(kernel_O, kernel_G, G(:,:,samp), X(:,samp), ...
        Task(1,samp), O(1,samp), L(1,samp), S.alpha, nL, dimG, dimX);
    
    % UPDATE O ------------------------------------------
    % first prior (posterior from last time step)
    if samp>I.n_zero_signal % accumulate evidence after signal starts
        pO = pO_Posterior(:,samp)'; % use exact posterior
        pL = pL_Posterior(:,samp)';
        prior_task = pT_Posterior(:,samp)';
    else
        pO = Ge.pO; % initial prior (same one, no accumulation)
        pL = Ge.pL;
        prior_task = Ge.prior_task;
    end
    if max(pO)<1 % uncertainty left about O
        log_like_O = zeros(Ge.number_orientations,1);
        for i = 1:nL
            if L(samp) ==i, jL = 1; else jL = 2; end
            log_like_O = log_like_O...
                +squeeze(sum(log(1-kernel_O(jL,Task(samp),:,G(i,:,samp) ==0)),4))... % grating off
                +squeeze(sum(log(  kernel_O(jL,Task(samp),:,G(i,:,samp) ==1)),4));   % grating on
            log_pO = Ge.odds_inc*log_like_O'+log(pO);
            pO = exp(log_pO-max(log_pO)); pO = pO/sum(pO);
        end
    end
    pO_Posterior(:,samp) = pO;
    if debug, disp(['i pO: ' num2str([samp pO])]); end
    if sum(isnan(pO))>0, error(num2str(pO)); end
    O(1,samp) = find(mnrnd(1,pO) ==1,1); %1+binornd(1,pO(2),1);
    % UPDATE L -------------------------------------------
    if max(pL)<1 % uncertainty left
        log_like_L = [0 0];
        for LL = 1:nL
            for i = 1:nL
                if LL ==i, jL = 1; else jL = 2; end
                log_like_L(LL) = log_like_L(LL)...
                    +sum(log(1-kernel_O(jL,Task(samp),O(samp),G(i,:,samp) ==0)))...
                    +sum(log(  kernel_O(jL,Task(samp),O(samp),G(i,:,samp) ==1)));
            end
        end
        log_pL = Ge.odds_inc*log_like_L+log(pL);
        pL = exp(log_pL-max(log_pL)); pL = pL/sum(pL);
    end
    pL_Posterior(:,samp) = pL;
    L(samp) = find(cumsum(pL)>rand(1),1,'first');
    % UPDATE Task ----------------------------------------
    if max(prior_task)<1 % uncertainty left about Task
        %log_like_T = S_ConditionalProbs('log_prior_task_gLO',G(:,:,i),L(i),O(1,i),phi_O,phi_g,P);
        log_like_T = [0 0];
        for i = 1:nL
            if L ==i, jL = 1; else jL = 2; end
            log_like_T = log_like_T...
                +squeeze(sum(log(1-kernel_O(jL,:,O(samp),G(i,:,samp) ==0)),4))...
                +squeeze(sum(log(  kernel_O(jL,:,O(samp),G(i,:,samp) ==1)),4));
        end
        log_prior_task = Ge.odds_inc*log_like_T+log(prior_task);
        prior_task = exp(log_prior_task-max(log_prior_task)); prior_task = prior_task/sum(prior_task);
    end
    pT_Posterior(:,samp) = prior_task;
    if debug, disp(['i prior_task: ' num2str([samp prior_task])]); end
    Task(samp) = find(mnrnd(1,prior_task) ==1,1); %1+binornd(1,prior_task(2),1); % 1 or 2
    % UPDATE Style ---------------------------------------
    xRx = X(:,samp)'*Ge.R*X(:,samp); if xRx>0, error('xRx must be <0!'); end
    yGx = cur_Image' * Ge.G * X(:,samp);
    if 0 % sparse prior
        sigs = sigy/sqrt(-xRx);
        mus = (yGx/sigy^2-1/Ge.tauStyle)*sigy^2/(-xRx);
    else % slow prior
        sigs = sqrt(1/(1/Ge.sigmaStyle^2-xRx/sigy^2));
        mus = sigs^2*(yGx/sigy^2+Style(samp)/Ge.sigmaStyle^2);
    end
    Style(samp) = Cut_Gaussian('random',mus,sigs,1);
end

if debug
    Get_Figure('Post-debug'); plot(pO_Posterior(2,:),'*-');
end
if debug
    Get_Figure('1'); Subplot(2,1,2,1);
    Histogram(TAU);
    Subplot(2); Histogram(MUK);
end

disp(size(O)) % 1xn_samples

O   (2:1+Ge.number_orientations,:) = pO_Posterior;
Task(2:1+Ge.nT,:) = pT_Posterior;
L   (2:1+nL,:) = pL_Posterior;

end

function newX = sample_X(img, Ge, kernel_G, X, G, s, sigy, alpha, nX, dimX)
%SAMPLE_X do a Gibbs sampling step on a random order of all X variables

newX = X;
% randomize update order for X
order = randperm(nX);
% X_k to be sampled
for i = order
    % x is drawn from a 'cut gaussian' with mean mu and variance sig^2
    % TODO - where did R_ii go (eqns 18-20 in 'ff vs fb' document)? Did we
    % ensure that each projective field has an L2 norm of 1?
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
    newX(i) = Cut_Gaussian('random', mu, sig, 1);
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
for jk = order, [j, k] = ind2sub([nL dimG], jk);
    if L == j, attn = 1; else attn = 2; end
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
        newG(j,k) = Gjk_val;
        
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
        tau = 1 + (newG(j,:) * kernel_G)';
        
        % x|g is drawn from an exponential distribution with mean 'tau', so
        % its log-likelihood is -log(tau) - x/tau
        log_p_each_case(case_idx) = ...
            log_p_each_case(case_idx) + ...
            sum(-log(tau) - X(:,j)./tau);
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
    newG(j,k) = binornd(1, pG(2), 1);
end
end

function newG = init_G(kernel_O, Task, O, L, nG, nL, dimG)
%INIT_G draw initial values for G based on inital values of O and L

newG = zeros(nL, dimG);
% order doesn't matter for initialization
for jk = 1:nG
    [j, k] = ind2sub([nL, dimG], jk);
    % if the location matches, use 'attended' kernel (slice 1) otherwise
    % the 'unattended' kernel (slice 2)
    if L == j, attn = 1; else attn = 2; end
    prob = kernel_O(attn, Task, O, k);
    newG(j,k) = binornd(1,prob,1); % 0 or 1
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