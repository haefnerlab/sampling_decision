function [X G O L Style] = Sampling_Gibbs(P,Prior,Input,NSamples,varargin)

sigy=1;

if nargin<5 || isempty(varargin{1}), access=ones(1,NSamples); else access=varargin{1}; end
% access determines which frames to use in the input
% ones means always the first frame
% 1:NSamples means a different one every time (dynamic random dots)

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

nx=P.nL*P.nX; % P.nx ist the x-size of the image

R=-P.G'*P.G;
if sum(abs(diag(R)+1)>1e-5), error('R_ii has to be -1!'); end

% initialization of samples

switch Prior
  case 'iid-exp'
    X=zeros(nx,NSamples);
    X(:,1)=random('exp',1,[nx 1]);
  case 'l-o' % location - orientation
    if nargin<6, pL=[0.5 0.5]; else pL=varargin{2}; end
    if nargin<7, pO=[0.5 0.5]; else pO=varargin{3}; end
    X=zeros(nx+2,NSamples); % + [L O]
    X(1:nx,1)=random('exp',1,[nx 1]);
    [~, X(nx+1,1)]=find(cumsum(pL)>rand(1),1,'first');
    [~, X(nx+2,1)]=find(cumsum(pO)>rand(1),1,'first');
    kLO=LO2K(X(nx+1),X(nx+2));
  case 'o-g' % orientation - grating
    S=varargin{2};
    phi_O=[pi/4 3/4*pi];
    ng=P.nG; phi_g=[0:ng-1]/ng*pi;
    pO=P.pO;
    X=zeros(nx,NSamples); %
    G=zeros(ng,NSamples);
    O=zeros( 1,NSamples);
    O(1)=1+binornd(1,pO(1)); % 1 and 2
    for k=1:ng
      G(k,1)=binornd(1,P.a+P.b*cos(2*(phi_O(O(1))-phi_g(k))),1); % 0 and 1
    end
    X(1:nx,1)=random('exp',1,[nx 1]); % not quite correct but easy and irrelevant after burn-in
  case 'l-o-g' % location-orientation - grating
    S=varargin{2};
    phi_O=[pi/4 3/4*pi];
    ng=P.nG; phi_g=[0:ng-1]/ng*pi;
    nL=P.nL; %phi_g=ones(nL,1)*phi_g;
    pO=P.pO; pL=P.pL;
    X=zeros(nx,NSamples); %
    G=zeros(nL,ng,NSamples);
    L=zeros( 1,NSamples);
    O=zeros( 1,NSamples);
    L(1)=find(cumsum(pL)>rand(1),1,'first');
    O(1)=1+binornd(1,pO(1)); % 1 and 2
    for j=1:nL
      if L(1)==j, beta=P.beta_hit; else beta=P.beta_miss; end
      for k=1:ng
        prob=beta*exp(P.kappa_O*(cos(2*(phi_O(O(1))-phi_g(k)))-1));
        G(j,k,1)=binornd(1,prob,1); % 0 and 1
      end
    end
    X(1:nx,1)=random('exp',1,[nx 1]); % not quite correct but easy and irrelevant after burn-in
  case 'l-o-g-s' % location-orientation - grating - style
    S=varargin{2};
    phi_O=P.phi_O;
    ng=P.nG; phi_g=[0:ng-1]/ng*pi;
    nL=P.nL; %phi_g=ones(nL,1)*phi_g;
    nX=P.nX; % number of X's per location
    pO=P.pO; pL=P.pL;
    X=zeros(nx,NSamples); %
    G=zeros(nL,ng,NSamples);
    Style=zeros( 1,NSamples);
    L=zeros( 1,NSamples);
    O=zeros( 1,NSamples);
    %Style(1)=gamrnd(P.paraStyle(1),P.paraStyle(2),1); disp('Style ~ Gamma!');
    Style(1)=exprnd(P.tauStyle,1); %disp('Style ~ Exp!');
    L(1)=find(cumsum(pL)>rand(1),1,'first');
    O(1)=1+binornd(1,pO(1)); % 1 and 2
    for j=1:nL
      if L(1)==j, beta=P.beta(2); else beta=0; end
      for k=1:ng
        prob=P.beta(1)+beta*exp(P.kappa_O*(cos(2*(phi_O(O(1))-phi_g(k)))-1));
        PROB(j,k)=prob;
        prob=min(0.9,prob);
        G(j,k,1)=binornd(1,prob,1); % 0 and 1
      end
    end
    for j=1:nL
      for i=1:P.nX
        tau=P.tau;
        for t2=1:P.nG
          tau=tau+(1-P.tau)*G(j,t2)*exp(P.kappa_G*(cos(2*(P.phi_x(i)-phi_g(t2)))-1));
        end
        X((j-1)*nX+i,1)=exprnd(tau,1);
      end
    end
    
  case 'd-o' % don't remember
    if nargin<6, kernelD='sin'; else kernelD=varargin{2}; end
    if nargin<7, kernelO='cos'; else kernelO=varargin{3}; end
    if nargin<8, pD=[0.5 0.5]; else pD=varargin{4}; end
    if nargin<9, pO=ones(1,nx)/nx; else pO=varargin{5}; end
    X=zeros(nx+2,NSamples);
    X(1:nx,1)=random('exp',1,[nx 1]);
    [~, X(nx+1,1)]=find(cumsum(pD)>rand(1),1,'first');
    [~, X(nx+2,1)]=find(cumsum(pO)>rand(1),1,'first');
    %kDO=DO2K(X(nx+1),X(nx+2));
    ks=0:nx-1;
    switch kernelD
      case 'sin', tauD_pre=P.tauD_amp*sin(2*(ks+0.5)/nx*pi)+1;
      otherwise, error(kernelD);
    end
    switch kernelO
      case 'cos', tauO_pre=P.tauO_amp*cos(2*ks/nx*pi)+1;
      otherwise, error(kernelO);
    end

  otherwise
    error(Prior);
end

% Gibbs sampling
for i=2:NSamples
  X(:,i)=X(:,i-1); % update to
  switch Prior
    case 'o-g',     G(:,i)=G(:,i-1); O(:,i)=O(:,i-1);
    case 'l-o-g',   
      error('correct?'); % not missing location dimension?
      G(:,i)=G(:,i-1); 
      O(:,i)=O(:,i-1); L(:,i)=L(:,i-1);
    case 'l-o-g-s'
      G(:,:,i)=G(:,:,i-1); 
      O(:,i)=O(:,i-1); L(:,i)=L(:,i-1); Style(i)=Style(i-1);
  end
  order=randperm(nx); % randomize update order
  for k=order % current source to be sampled
    idx_no_k=[1:k-1 k+1:nx]; % all but k index
    switch Prior
      case 'iid-exp'
        sigk=sigy;
        muk=P.G(:,k)'*Input(:,access(i))+R(k,idx_no_k)*X(idx_no_k,i)-2*sigk^2; % my calc - confirmed
        X(k,i)=Cut_Gaussian('random',muk,sigk,1);
        while isnan(X(k,i)) || abs(X(k,i))>1e10
          disp(['resampling X(' num2str([k i]) ')']);
          X(k,i)=Cut_Gaussian('random',muk,sigk,1);
        end
      case 'l-o'
        if k<=nx
          sigk=sigy;
          if k==kLO
            muk=P.G(:,k)'*Input(:,access(i))+R(k,idx_no_k)*X(idx_no_k,i)-2*sigk^2; % my calc - confirmed by PB
          else
            muk=P.G(:,k)'*Input(:,access(i))+R(k,idx_no_k)*X(idx_no_k,i)-2*sigk^2/P.tau; % my calc - confirmed by PB
          end
          X(k,i)=Cut_Gaussian('random',muk,sigk,1);
          while isnan(X(k,i)) || abs(X(k,i))>1e10
            disp(['resampling X(' num2str([k i]) ')']);
            X(k,i)=Cut_Gaussian('random',muk,sigk,1);
          end
        elseif k==nx+1
          p=Aux_PL(X(1:nx,i),X(nx+2,i),P.tau);
          X(k,i)=find(cumsum(p)>rand(1),1,'first');
          kLO=LO2K(X(nx+1,i),X(nx+2,i));
        else
          p=Aux_PO(X(1:nx,i),X(nx+1,i),P.tau);
          X(k,i)=find(cumsum(p)>rand(1),1,'first');
          kLO=LO2K(X(nx+1,i),X(nx+2,i));
        end
        
      case 'o-g'
        sigk=sigy;
        tau=S_ConditionalProbs('tau-kg',k,G(:,i),P.tau0,P.tau_amp,P.kappa,P.phi_x,phi_g);
        muk=P.G(:,k)'*Input(:,access(i))+R(k,idx_no_k)*X(idx_no_k,i)-2*sigk^2/tau; % my calc - confirmed by PB
        X(k,i)=Cut_Gaussian('random',muk,sigk,1);
        while isnan(X(k,i)) || abs(X(k,i))>1e10
          disp(['resampling X(' num2str([k i]) ')']);
          X(k,i)=Cut_Gaussian('random',muk,sigk,1);
        end
        
      case 'l-o-g'
        sigk=sigy;
        kk=mod(k-1,P.nL)+1;
        tau=S_ConditionalProbs('tau-kg',kk,G(:,:,i),P.tau0,P.tau_amp,P.kappa_O,P.phi_x,phi_g);
        muk=P.G(:,k)'*Input(:,access(i))+R(k,idx_no_k)*X(idx_no_k,i)-2*sigk^2/tau; % my calc - confirmed by PB
        %Dtau(k)=tau; Dmuk(k)=muk;
        %Dgi(k)=P.G(:,k)'*Input(:,access(i));
        %Drx(k)=R(k,idx_no_k)*X(idx_no_k,i);
        X(k,i)=Cut_Gaussian('random',muk,sigk,1);
        while isnan(X(k,i)) || abs(X(k,i))>1e10
          disp(['resampling X(' num2str([k i]) ')']);
          X(k,i)=Cut_Gaussian('random',muk,sigk,1);
        end
        
      case 'l-o-g-s'
        sigk=sigy;
        kO=mod(k-1,P.nX)+1; % index of current grating orientation
        kL=1+(k-kO)/P.nX; % index of current location
        tau=S_ConditionalProbs('tau-kg',kO,G(kL,:,i),P.tau,P.kappa_O,P.phi_x,phi_g);
        muk=Style(i)*P.G(:,k)'*Input(:,access(i))+Style(i)^2*R(k,idx_no_k)*X(idx_no_k,i)-sigk^2/tau; % my calc - confirmed by PB
        TAU(k)=tau; MUK(k)=muk;
        X(k,i)=Cut_Gaussian('random',muk,sigk,1);
        if X(k,i)>100 || isnan(X(k,i))
          disp(['X=' num2str(X(k,i))]);
          X(k,i)=100;
        end
        while isnan(X(k,i)) || abs(X(k,i))>1e10
          disp(['resampling X(' num2str([k i]) ')=' num2str(X(k,i))]);
          X(k,i)=Cut_Gaussian('random',muk,sigk,1);
        end
        
      case 'd-o'
        sigk=sigy;
        if k==1
          dd=X(nx+1,i); oo=X(nx+2,i);
          tauO=circshift(tauO_pre,[0 oo]);
          if dd==1, tau=tauO.*tauD_pre; else tau=tauO.*(2-tauD_pre); end
        end
        if k<=nx
          muk=P.G(:,k)'*Input(:,access(i))+R(k,idx_no_k)*X(idx_no_k,i)-2*sigk^2/tau(k);
          X(k,i)=Cut_Gaussian('random',muk,sigk,1);
          while isnan(X(k,i)) || abs(X(k,i))>1e10
            disp(['resampling X(' num2str([k i]) ')']);
            X(k,i)=Cut_Gaussian('random',muk,sigk,1);
          end
        elseif k==nx+1
          p=Aux_DO_PD(X(1:nx,i));
          X(k,i)=find(cumsum(p)>rand(1),1,'first');
        else
          p=Aux_DO_PO(X(1:nx,i),X(nx+1,i));
          X(k,i)=find(cumsum(p)>rand(1),1,'first');
        end
        
      otherwise
        error(Prior);
    end
  end
  switch Prior
    case 'o-g'
      for k=1:ng
        pG=S_ConditionalProbs('pG_kgOx',k,G(:,i),O(i),X(:,i),phi_O,phi_g,P.a,P.b);
        G(k,i)=binornd(1,pG(2),1); % 0 and 1
      end
      pO=S_ConditionalProbs('pO_g',G(:,i),phi_O,phi_g,P.a,P.b,P.pO);
      O(i)=1+binornd(1,pO(2),1); % 1 and 2
    case 'l-o-g'
      for j=1:nL
        for k=1:ng
          pG=S_ConditionalProbs('pG_kgLOx',[j k],G(:,:,i),L(i),O(i),X(:,i),phi_O,phi_g,P);
          G(j,k,i)=binornd(1,pG(2),1); % 0 and 1
        end
        pO=S_ConditionalProbs('pO_gL',G(:,:,i),L(i),phi_O,phi_g,P);
        O(i)=1+binornd(1,pO(2),1); % 1 and 2
        pL=S_ConditionalProbs('pL_gO',G(:,:,i),O(i),phi_O,phi_g,P);
        L(i)=find(cumsum(pL)>rand(1),1,'first');
      end
    case 'l-o-g-s'
      for j=1:nL
        for k=1:ng
          pG=S_ConditionalProbs('pG_kgLOx',[j k],G(:,:,i),L(i),O(i),X(:,i),phi_O,phi_g,P);
          PG(j,k)=pG(2);
          if sum(isnan(pG))>0
            error('NaN');
          end
          G(j,k,i)=binornd(1,pG(2),1); % 0 or 1
        end
      end
      pO=S_ConditionalProbs('pO_gL',G(:,:,i),L(i),phi_O,phi_g,P);
      O(i)=1+binornd(1,pO(2),1); % 1 and 2
      pL=S_ConditionalProbs('pL_gO',G(:,:,i),O(i),phi_O,phi_g,P);
      L(i)=find(cumsum(pL)>rand(1),1,'first');
      Gx2=-X(:,i)'*R*X(:,i); if Gx2<0, error('Gx2 must be >0!'); end
      sigs=sigy/sqrt(Gx2);
      mus=(Input(:,access(i))'*P.G*X(:,i)... % y^t G x
        -sigy^2/P.tauStyle)/Gx2;          %
      Style(i)=Cut_Gaussian('random',mus,sigs,1);
      if Style(i)>1e2
        disp(['Style=' num2str(Style(i))]);
        Style(i)=1e2;
      end
      while isnan(Style(i)) || abs(Style(i))>1e10
        disp(['resampling Style(' num2str(i) ')=' num2str(s(i))]);
        s(i)=Cut_Gaussian('random',muk,sigy,1);
      end
    case 'l-o-g-s-vec'
      sigk=sigy;
      % replace k by vectorization!
      kO=mod(k-1,P.nX)+1; % index of current grating orientation
      kL=1+(k-kO)/P.nX; % index of current location
      %tau=S_ConditionalProbs('tau-kg',kO,G(kL,:,i),P.tau,P.kappa_O,P.phi_x,phi_g);
      tau=S_ConditionalProbs('tau-kg',kO,G(:,:,i),P.tau,P.kappa_O,P.phi_x,phi_g);
      %muk=Style(i)*P.G(:,k)'*Input(:,access(i))+Style(i)^2*R(k,idx_no_k)*X(idx_no_k,i)-sigk^2/tau; % my calc - confirmed by PB
      muk=Style(i)*P.G(:,:)'*Input(:,access(i))+Style(i)^2*R(:,:)*X(:,i)-sigk^2/tau; % my calc - confirmed by PB
      X(:,i)=Cut_Gaussian('random',muk,sigk,[nx 1]);
      %if X(k,i)>100 || isnan(X(k,i))
      %  disp(['X=' num2str(X(k,i))]);
      %  X(k,i)=100;
      %end
      X(X(:,i)>100,i)=100;
      X(isnan(X(:,i)),i)=100;
      % end of vectorization of X
      for j=1:nL
        for k=1:ng
          pG=S_ConditionalProbs('pG_kgLOx',[j k],G(:,:,i),L(i),O(i),X(:,i),phi_O,phi_g,P);
          if sum(isnan(pG))>0
            error('NaN');
          end
          G(j,k,i)=binornd(1,pG(2),1); % 0 or 1
        end
      end
      pO=S_ConditionalProbs('pO_gL',G(:,:,i),L(i),phi_O,phi_g,P);
      O(i)=1+binornd(1,pO(2),1); % 1 and 2
      pL=S_ConditionalProbs('pL_gO',G(:,:,i),O(i),phi_O,phi_g,P);
      L(i)=find(cumsum(pL)>rand(1),1,'first');
      Gx2=-X(:,i)'*R*X(:,i); if Gx2<0, error('Gx2 must be >0!'); end
      sigs=sigy/sqrt(Gx2);
      mus=(Input(:,access(i))'*P.G*X(:,i)... % y^t G x
        -sigy^2/P.tauStyle)/Gx2;          %
      Style(i)=Cut_Gaussian('random',mus,sigs,1);
      if Style(i)>1e2
        disp(['Style=' num2str(Style(i))]);
        Style(i)=1e2;
      end
      while isnan(Style(i)) || abs(Style(i))>1e10
        disp(['resampling Style(' num2str(i) ')=' num2str(s(i))]);
        s(i)=Cut_Gaussian('random',muk,sigy,1);
      end
  end
  if isnan(Style(i)) sum(isnan(X(:,i)))>0 || sum(sum(isnan(G(:,:,i))))>0 || sum(isnan([L(i) O(i)]))>0
    error(i);
  end
end

%disp('end Gibbs sampling');

  function kk = LO2K(L,O)
    kk=(L-1)*2+O; % 
  end

  function aux = Aux_PL(x,O,tau)
    xx=reshape(x,[2 2])';
    aux=zeros(1,2);
    s=sum(xx(:));
    for LL=[1 2]
      aux(LL)=1/tau^3*exp(-xx(LL,O))*exp(-(s-xx(LL,O))/tau)*pL(LL);
    end
    aux=aux/sum(aux);
  end
  function aux = Aux_PO(x,L,tau)
    xx=reshape(x,[2 2])';
    aux=zeros(1,2);
    s=sum(xx(:));
    for OO=[1 2]
      aux(OO)=1/tau^3*exp(-xx(L,OO))*exp(-(s-xx(L,OO))/tau)*pO(OO);
    end
    aux=aux/sum(aux);
  end

  function aux = Aux_DO_PD(x)
    aux=zeros(1,2);
    for D=[1 2]
      if D==1, taux=tauO.*tauD_pre; else taux=tauO.*(2-tauD_pre); end
      aux(D)=1/prod(taux)*exp(-sum(x./taux'))*pD(D);
    end
    aux=aux/sum(aux);
  end
  function aux = Aux_DO_PO(x,D)
    if D==1, tauD=tauD_pre; else tauD=2-tauD_pre; end
    aux=zeros(1,nx);
    for OO=1:nx
      taux=tauD.*circshift(tauO_pre,[0 OO]);
      aux(OO)=1/prod(taux)*exp(-sum(x./taux'))*pO(OO);
    end
    aux=aux/sum(aux);
  end

end