%Inputs are Ge: Generative Model, S for sampling information, I for Image
%related information, Input (the image) access and zero signal.n

function [X G O L Style Task] = Sampling_Gibbs_InPlace_Fast(Ge, S, I, Input,access, n_zero_signal)

% created as fast version of original Sampling_Gibbs_InPlace
% incorporates S_ConditionalProbs (or at least 'tau' and 'pO_...'

DEBUG=0;

sigy=1;


% access determines which frames to use in the input


nI=size(Input);
switch length(nI)
  case 2
    Input=reshape(Input,[numel(Input) 1]);
  case 3
    aux=zeros(nI(2)*nI(3),nI(1));
    for i=1:nI(1)
      aux(:,i)=reshape(Input(i,:,:),[nI(2)*nI(3) 1]);
    end
    Input=aux;
  otherwise
    error(size(Input));
end

NX=Ge.number_locations*Ge.dimension_X; % Number of X variables, <>P.nx which is the x-size of the image

R=-Ge.G'*Ge.G; % Ge.G contains projective fields of X
if sum(abs(diag(R)+1)>1e-5), error('R_ii has to be -1!'); end

% initialization of kernels
kernel_G=zeros(Ge.dimension_G,Ge.dimension_X);
C=1/besseli(0,Ge.kappa_G);
for i=1:Ge.dimension_G
  for j=1:Ge.dimension_X
    %kernel_G(i,j)=Ge.delta/G.dimension_X*(C*exp(G.kappa_G*cos(2*(G.phi_x(j)-G.phi_g(i))))); % OLD
    kernel_G(i,j)=Ge.delta*(C*exp(Ge.kappa_G*cos(2*(Ge.phi_x(j)-Ge.phi_g(i))))); % July 9, 2015
  end
end
% check for worst case scenario
%aux=zeros(1,G.dimension_G); aux(kernel_G(:,1)<0)=1;
%if aux*kernel_G(:,1)<1, error('kernel too negative!'); end
%
kernel_O=zeros([2 G.nT, Ge.number_orientations Ge.dimension_G]);
% 1st dim: 1 at attended location, 2 different location
C1=1/besseli(0,Ge.kappa_O(1))/Ge.dimension_G; 
C2=1/besseli(0,Ge.kappa_O(2))/Ge.dimension_G;
for T=1:G.nT
  for O=1:Ge.number_orientations
    switch Ge.task
      case 'discrimination', delta_cos=cos(2*(G.phi_O(T,O)-Ge.phi_g));
      case 'detection'
        if O==1,             delta_cos=0; % stimulus absent
        else                 delta_cos=cos(2*(G.phi_O(T,O)-Ge.phi_g));
        end
    end
    kernel_O(1,T,O,:)=1*C1*exp(Ge.kappa_O(1)*delta_cos);
    kernel_O(2,T,O,:)=1*C2*exp(Ge.kappa_O(2)*delta_cos);
  end
end
kernel_O(kernel_O>1)=1;

% initialization of samples
pO=G.pO; pL=G.pL; prior_task=Ge.prior_task;

X=zeros(NX,S.n_samples); %
G=zeros(Ge.number_locations,Ge.dimension_G,S.n_samples);
Style=zeros(1,S.n_samples);
L=zeros(   1,S.n_samples);
Style(1)=1; % start with Olshausen & Field
X=zeros(NX, S.n_samples); %
G=zeros(Ge.number_locations,Ge.dimension_G, S.n_samples);
Style=zeros(1,S.n_samples);
L=zeros(   1,S.n_samples);
Style(1)=G.tauStyle; % start with prior mean



L(1)=find(cumsum(pL)>rand(1),1,'first');
O=zeros(1,S.n_samples); % 1st: sample, 2nd: ratio 2/1
O(1)=find(mnrnd(1,pO)==1,1); %1+binornd(1,pO(2));
pO_Posterior=zeros(Ge.number_orientations,S.n_samples);
pO_Posterior(:,1)=pO;
Task =zeros(1,S.n_samples); % 1st: sample, 2nd: ratio 2/1
prior_task_Posterior=zeros(G.nT,S.n_samples);
prior_task_Posterior(:,1)=prior_task;
pL_Posterior=zeros(Ge.number_locations,S.n_samples);
pL_Posterior(:,1)=pL;
Task(1)=find(mnrnd(1,prior_task)==1,1); %1+binornd(1,prior_task(2));


i=1; % first sample
% --------- sampling G's---------------
order=randperm(Ge.number_locations*Ge.dimension_G); % randomize update order for G
for jk=order, [j k]=ind2sub([Ge.number_locations Ge.dimension_G],jk);
  if L(i)==j, iL=1; else iL=2; end
  prob=kernel_O(iL,Task(1,i),O(1,i),k);
  G(j,k,i)=binornd(1,prob,1); % 0 or 1
end
% --------- sampling X's---------------
for k=1:NX % current X to be sampled, order doesn't make a difference!
  kO=mod(k-1,Ge.dimension_X)+1; % index of current Gabor orientation
  kL=1+(k-kO)/Ge.dimension_X; % index of current location
  tau=1+G(kL,:,i)*kernel_G(:,kO); % simplification July 2013
  %tau=Ge.delta+G(kL,:,i)*kernel_G(:,kO); % modification July 2015
  X(k,i)=exprnd(tau);
end
% Gibbs sampling ---------------------------------------------------------

for i=2:S.n_samples
  % Copy entire current state forward by 1 step
  X(:,i)=X(:,i-1);
  G(:,:,i)=G(:,:,i-1);
  O(:,i)=O(:,i-1);
  L(:,i)=L(:,i-1);
  Task(:,i)=Task(:,i-1);
  Style(:,i)=Style(:,i-1);
  pO_Posterior(:,i)=pO_Posterior(:,i-1);
  pL_Posterior(:,i)=pL_Posterior(:,i-1);
  prior_task_Posterior(:,i)=prior_task_Posterior(:,i-1);


  % UPDATE X
  order=randperm(NX); % randomize update order for X
  for k=order % current X to be sampled
    idx_no_k=[1:k-1 k+1:NX]; % all but k index
    sigk=sigy/Style(i); % realized July 2013
    kO=mod(k-1,Ge.dimension_X)+1; % index of current Gabor orientation
    kL=1+(k-kO)/Ge.dimension_X; % index of current location
    tau=1+G(kL,:,i)*kernel_G(:,kO); % simplification July 2013
    %tau=Ge.delta+G(kL,:,i)*kernel_G(:,kO); % modification July 2015
    if tau<=0, error(['tau<=0']); end
    muk=Style(i)*Ge.G(:,k)'*Input(:,access(i))+Style(i)^2*R(k,idx_no_k)*X(idx_no_k,i)-sigy^2/tau; % my calc - confirmed by PB
    muk=muk/Style(i)^2; % realized July 2013
    X(k,i)=Cut_Gaussian('random',muk,sigk,1);
    if ~isreal(X(k,i)), error(['X NaN, mu sigma: ' num2str([muk sigk])]); end
  end
  % UPDATE everything else


  % UPDATE G ------------------------------------------
  order=randperm(Ge.number_locations*Ge.dimension_G); % randomize update order for G
  for jk=order, [j k]=ind2sub([Ge.number_locations Ge.dimension_G],jk);
    %pG=S_ConditionalProbs('pG_kgTLOx',[j k],G(:,:,i),Task(1,i),L(i),O(1,i),X(:,i),phi_O,phi_g,P);
    if L(i)==j, iL=1; else iL=2; end
    prob=kernel_O(iL,Task(1,i),O(1,i),k);
    GG=G(:,:,i);
    aux=zeros(1,2);
    for glk=[0 1]
      GG(j,k)=glk;
      switch glk
        case 0, aux(1)=log(1-prob); % grating off
        case 1, aux(2)=log(prob); % grating is on
      end
      tau=1+(GG(j,:)*kernel_G)';
      %tau=Ge.delta+(GG(j,:)*kernel_G)'; % modification July 2015
      aux(1+glk)=aux(1+glk)-sum(log(tau)+X((j-1)*Ge.dimension_X+(1:Ge.dimension_X),i)./tau);

      %ERROR HANDLING
      if ~isreal(aux(1+glk))
        disp(['tau : ' num2str(tau')]);
        disp(['aux : ' num2str(aux)]);
        disp(['sum(.): ' num2str((log(tau)+X((j-1)*Ge.dimension_X+(1:G.dimension_X),i)./tau)')]);
        error('aux NaN');
      end

    end


    aux=exp(aux-max(aux));


    pG=aux/sum(aux);

    %ERROR HANDLING
    if sum(~isreal(pG))>0,
      disp(['tau : ' num2str(tau')]);
      disp(['aux : ' num2str(aux)]);
      error('pG NaN'); 
    end

    if pG(2)>1, error('pG(2)>1'); end



    G(j,k,i)=binornd(1,pG(2),1); % 0 or 1
  end


  %G(1,:,i)=[0 0 0 0]; % useful for debugging! seems to work!
  % UPDATE O ------------------------------------------
  % first prior (posterior from last time step)
  if i>n_zero_signal % accumulate evidence after signal starts
    pO=pO_Posterior(:,i)'; % use exact posterior
    pL=pL_Posterior(:,i)';
    prior_task=prior_task_Posterior(:,i)';
  else
    pO=G.pO; % initial prior (same one, no accumulation)
    pL=G.pL;
    prior_task=Ge.prior_task;
  end
  if max(pO)<1 % uncertainty left about O
    log_like_O=zeros(Ge.number_orientations,1);
    for j=1:Ge.number_locations
      if L(i)==j, jL=1; else jL=2; end
      log_like_O=log_like_O...
        +squeeze(sum(log(1-kernel_O(jL,Task(i),:,G(j,:,i)==0)),4))... % grating off
        +squeeze(sum(log(  kernel_O(jL,Task(i),:,G(j,:,i)==1)),4));   % grating on
      log_pO=G.odds_inc*log_like_O'+log(pO);
      pO=exp(log_pO-max(log_pO)); pO=pO/sum(pO);
    end
  end  
  pO_Posterior(:,i)=pO;
  if DEBUG, disp(['i pO: ' num2str([i pO])]); end
  if sum(isnan(pO))>0, error(num2str(pO)); end
  O(1,i)=find(mnrnd(1,pO)==1,1); %1+binornd(1,pO(2),1);
  % UPDATE L -------------------------------------------
  if max(pL)<1 % uncertainty left
    log_like_L=[0 0];
    for LL=1:Ge.number_locations
      for j=1:Ge.number_locations
        if LL==j, jL=1; else jL=2; end
        log_like_L(LL)=log_like_L(LL)...
          +sum(log(1-kernel_O(jL,Task(i),O(i),G(j,:,i)==0)))...
          +sum(log(  kernel_O(jL,Task(i),O(i),G(j,:,i)==1)));
      end
    end
    log_pL=G.odds_inc*log_like_L+log(pL);
    pL=exp(log_pL-max(log_pL)); pL=pL/sum(pL);
  end
  pL_Posterior(:,i)=pL;
  L(i)=find(cumsum(pL)>rand(1),1,'first');
  % UPDATE Task ----------------------------------------
  if max(prior_task)<1 % uncertainty left about Task
    %log_like_T=S_ConditionalProbs('log_prior_task_gLO',G(:,:,i),L(i),O(1,i),phi_O,phi_g,P);
    log_like_T=[0 0];
    for j=1:Ge.number_locations
      if L==j, jL=1; else jL=2; end
      log_like_T=log_like_T...
        +squeeze(sum(log(1-kernel_O(jL,:,O(i),G(j,:,i)==0)),4))...
        +squeeze(sum(log(  kernel_O(jL,:,O(i),G(j,:,i)==1)),4));
    end
    log_prior_task=G.odds_inc*log_like_T+log(prior_task);
    prior_task=exp(log_prior_task-max(log_prior_task)); prior_task=prior_task/sum(prior_task);
  end
  prior_task_Posterior(:,i)=prior_task;
  if DEBUG, disp(['i prior_task: ' num2str([i prior_task])]); end
  Task(i)=find(mnrnd(1,prior_task)==1,1); %1+binornd(1,prior_task(2),1); % 1 or 2
  % UPDATE Style ---------------------------------------
  xRx=X(:,i)'*R*X(:,i); if xRx>0, error('xRx must be <0!'); end
  yGx=Input(:,access(i))'*Ge.G*X(:,i);
  if 0 % sparse prior
    sigs=sigy/sqrt(-xRx);
    mus=(yGx/sigy^2-1/G.tauStyle)*sigy^2/(-xRx);
  else % slow prior
    sigs=sqrt(1/(1/G.sigmaStyle^2-xRx/sigy^2));
    mus=sigs^2*(yGx/sigy^2+Style(i)/G.sigmaStyle^2);
  end
  Style(i)=Cut_Gaussian('random',mus,sigs,1);
end

if DEBUG
  Get_Figure('Post-debug'); plot(pO_Posterior(2,:),'*-'); 
end
if DEBUG
  Get_Figure('1'); Subplot(2,1,2,1);
  Histogram(TAU);
  Subplot(2); Histogram(MUK);
end

O   (2:1+Ge.number_orientations,:)=pO_Posterior;
Task(2:1+G.nT,:)=prior_task_Posterior;
L   (2:1+Ge.number_locations,:)=pL_Posterior;
%disp('end Gibbs sampling');
