function [Y noiseless]= InputImage(fct,varargin)

DEBUG = 0;

switch fct
    case 'binary'
        nx=2;
        c=varargin{1};
        if nargin<3, nr=1; else nr=varargin{2}; end
        Y=randn(nr,nx); % white noise background
        if numel(c)==2
            Y=Y+ones(nr,1)*reshape(c,1,2);
            noiseless=reshape(c,1,2);
        else
            Y=Y+c;
            noiseless=[max(c(:)) min(c(:))];
        end
        if DEBUG, debug_image; end
        
    case '1x2' % 1 location, 2 orientations
        nx=varargin{1}; ny=nx;
        c =varargin{2}; % contrast of each patch [LV LH]
        if nargin<4, nr=1; else nr=varargin{3}; end
        %Y=1-2*round(rand(nr,ny,nx)); % RDS background
        Y=squeeze(randn(nr,ny,nx)); % white noise background
        x=linspace(-1,1,nx); y=linspace(-1,1,ny);
        [xx yy]=meshgrid(x,y);
        ph=0;
        LV=Gabor_neu([0 1 0.5 ph 0 0.2],xx,'orig').*normpdf(yy, 0,0.6);
        LH=Gabor_neu([0 1 0.5 ph 0 0.2],yy,'orig').*normpdf(xx, 0,0.6);
        if numel(c)==2
            aux=c(1)*LV/norm(LV)+c(2)*LH/norm(LH);
            Y=Y+shiftdim(reshape(repmat(aux,[1 nr]),[nx ny nr]),2);
            noiseless=aux;
        else
            for i=1:nr
                aux=c(i,1)*LV/norm(LV)+c(i,2)*LH/norm(LH);
                Y(i,:,:)=squeeze(Y(i,:,:))+aux;
            end
            noiseless=max(c(1,:))*LV/norm(LV);
        end
        if DEBUG, debug_image; end
        Y=squeeze(Y);
        
    case '2x2'
        nx=varargin{1}; ny=nx;
        c =varargin{2}; % contrast of each patch [LV LH RV RH]
        if nargin<4, nr=1; else nr=varargin{3}; end
        %Y=1-2*round(rand(nr,ny,nx)); % RDS background
        Y=randn(nr,ny,nx); % white noise background
        x=linspace(-2,2,nx); y=linspace(-2,2,ny);
        [xx yy]=meshgrid(x,y);
        LV=Gabor_neu([0 1 1 pi/2 -1 0.2],xx,'orig').*normpdf(yy, 0,0.3);
        RV=Gabor_neu([0 1 1 pi/2  1 0.2],xx,'orig').*normpdf(yy, 0,0.3);
        LH=Gabor_neu([0 1 1 pi/2  0 0.2],yy,'orig').*normpdf(xx,-1,0.3);
        RH=Gabor_neu([0 1 1 pi/2  0 0.2],yy,'orig').*normpdf(xx, 1,0.3);
        aux=c(1)*LV/norm(LV)+c(2)*LH/norm(LH)+c(3)*RV/norm(RV)+c(4)*RH/norm(RH);
        Y=squeeze(Y)+shiftdim(reshape(repmat(aux,[1 nr]),[nx ny nr]),2);
        if DEBUG, debug_image; end
        %Y=squeeze(Y);
        noiseless=aux;
        
    case 'nx2' % for the time being there is no nxN
        nL = varargin{1};
        ny = varargin{2}; nx=nL*ny;
        c  = varargin{3}; % contrast of each patch [time,location,orientation]
        %%% By Shizhao Liu 03/05/2025:
        %%% Add another input argument "image_task", which controls whether
        %%% the images are cardinal (0/90 degrees, originally implemented
        %%% by this code) or oblique (45/135 degrees. implement this by
        %%% rotating the coordinates by 45 degrees)
        %%% Note this was only done for this nx2 option. 
        image_task = varargin{4};
        switch length(size(c))
            case 2,    nr=1; c=reshape(c,[1 size(c)]);
            case 3,    nr=size(c,1);
            otherwise, error('wrong size contrast dimensions');
        end
        Y=randn(nr,ny,nx); % white noise background
        %Y=0*Y; warning('zero noise!!!');
        x=linspace(-1/2,nL-1/2,nx); y=linspace(-1/2,1/2,ny);
        [xx yy]=meshgrid(x,y);
        switch image_task 
            case 'cardinal'
                % no need to rotate, just use the orginal  coordinates
                xxr = xx;
                yyr = yy;
            case 'oblique'
                % rotate the pixel coordinates we got from meshgrid
                rotate_phi = pi/4;
                rot = [[cos(rotate_phi) sin(rotate_phi)]; [-sin(rotate_phi) cos(rotate_phi)]];
                zza = rot * [xx(:)'; yy(:)'];
                xxr = reshape(zza(1,:), ny,nx);
                yyr = reshape(zza(2,:), ny,nx);
        end
        for i=1:nL
            %auxV=Gabor_neu([0 1 2 0  i-1 0.1],xx,'orig').*normpdf(yy,  0,0.2);
            auxV=Gabor_neu([0 1 2 0  i-1 0.1],xxr,'orig'); % new July 2015
            auxV=auxV/norm(auxV);
            %auxH=Gabor_neu([0 1 2 0    0 0.1],yy,'orig').*normpdf(xx,i-1,0.2);
            auxH=Gabor_neu([0 1 2 0    0 0.1],yyr,'orig'); % new July 2015
            auxH=auxH/norm(auxH);
            for k=1:nr
                Y(k,:,:)=squeeze(Y(k,:,:))+c(k,i,1)*auxV+c(k,i,2)*auxH;
            end
        end
        Y=squeeze(Y);
        if DEBUG, debug_image; end
        % if nargout<1, noiseless=aux; end % don't understand this statement!
        
    case 'nx2-OLD' % for the time being there is no nxN
        nL=varargin{1};
        ny=varargin{2}; nx=nL*ny;
        c =varargin{3}; % contrast of each patch [location,orientation]
        if nargin<5, nr=1; else nr=varargin{4}; end
        Y=zeros(ny,nx);
        Noise=randn(nr,ny,nx); % white noise background
        x=linspace(-1/2,nL-1/2,nx); y=linspace(-1/2,1/2,ny);
        [xx yy]=meshgrid(x,y);
        V=zeros(nL,nx,ny); H=zeros(size(V));
        for i=1:nL
            auxG=Gabor_neu([0 1 2 0  i-1 0.1],xx,'orig').*normpdf(yy,  0,0.2);
            Y=Y+c(i,1)*auxG/norm(auxG);
            auxG=Gabor_neu([0 1 2 0    0 0.1],yy,'orig').*normpdf(xx,i-1,0.2);
            Y=Y+c(i,2)*auxG/norm(auxG);
        end
        Y=squeeze(Noise)+shiftdim(reshape(repmat(Y,[1 nr]),[ny nx nr]),2);
        %Y=squeeze(Y);
        if DEBUG, debug_image; end
        if nargout<1, noiseless=aux; end
        
    case 'logs-gen' % same as nx2 but based on true generative model
        ny=varargin{1};
        S=varargin{2};
        phi_O=[pi/4 3/4*pi];
        ng=S.nG; phi_g=[0:ng-1]/ng*pi;
        phi=(0:S.nX-1)/S.nX*pi;
        nL=S.nL; %phi_g=ones(nL,1)*phi_g;
        pO=G.pO; pL=G.pL;
        %s=gamrnd(S.paraStyle(1),S.paraStyle(2),1); disp('s ~ Gamma!');
        s=exprnd(S.tau_Style,1); disp('s ~ Exp!');
        L(1)=find(cumsum(G.pL)>rand(1),1,'first');
        O(1)=1+binornd(1,G.pO(1)); % 1 and 2
        for j=1:S.nL
            if L(1)==j, beta=S.beta(2); else beta=S.beta(1); end
            for k=1:ng
                prob=beta*exp(S.kappa_O*(cos(2*(phi_O(O(1))-phi_g(k)))-1));
                G(j,k)=binornd(1,prob,1); % 0 and 1
            end
        end
        nx=nL*ny;
        Y=randn(ny,nx); % white noise background
        x=linspace(-1/2,nL-1/2,nx); y=linspace(-1/2,1/2,ny);
        [xx yy]=meshgrid(x,y);
        for j=1:nL
            cx=j-1; % center
            [xx yy]=meshgrid(x-cx,y);
            for i=1:S.nX
                tau=S.tau;
                for t2=1:S.nG
                    tau=tau+(1-S.tau)*G(j,t2)*exp(S.kappa_G*(cos(2*(phi(i)-phi_g(t2)))-1));
                end
                X(j,i)=exprnd(tau,1);
                %X(j,i)=tau;
                rot=[[cos(phi(i)) sin(phi(i))];[-sin(phi(i)) cos(phi(i))]];
                zza=rot*[xx(:)'; yy(:)'];
                xxs=zza(1,:); yys=zza(2,:);
                auxG=Gabor_neu([0 1 2 0 0 0.1],xxs,'orig').*normpdf(yys,0,0.2);
                Y=Y+s*X(j,i)*reshape(auxG,[ny nx]);
            end
        end
        %Y=squeeze(Y)+shiftdim(reshape(repmat(aux,[1 nr]),[nx ny nr]),2);
        if DEBUG, debug_image; end
        %Y=squeeze(Y);
        if nargout>1, noiseless=G; end
        
    case '1xN'
        nx=varargin{1}; ny=nx;
        no=varargin{2};
        c =varargin{3}; % contrast of each patch
        if nargin<5, nr=1; else nr=varargin{4}; end
        Y=randn(nr,ny,nx); % white noise background
        x=linspace(-1/2,1/2,nx); y=linspace(-1/2,1/2,ny);
        [xx yy]=meshgrid(x,y);
        phi=(0:no-1)/no*pi;
        aux=zeros(1,nx*ny);
        for i=1:no
            rot=[[cos(phi(i)) sin(phi(i))];[-sin(phi(i)) cos(phi(i))]];
            zza=rot*[xx(:)'; yy(:)']; xxs=zza(1,:); yys=zza(2,:);
            auxG=Gabor_neu([0 1 2 0 0 0.1],xxs,'orig').*normpdf(yys,0,2);
            aux=aux+c(i)*auxG;
        end
        aux=reshape(aux,nx,ny);
        Y=squeeze(Y)+squeeze(shiftdim(reshape(repmat(aux,[1 nr]),[nx ny nr]),2));
        if DEBUG, debug_image; end
        Y=squeeze(Y);
        noiseless=aux;
    otherwise
        error(fct);
end





    function debug_image
        Get_Figure('II');
        switch fct
            case 'logs-gen'
                disp(['s=' num2str(s)]);
                disp(['L=' num2str(L)]);
                disp(['O=' num2str(O)]);
                disp(G);
                disp(X);
            otherwise
                title(['input: ' num2str(c(:)')]);
                if length(size(Y))==3
                    imagesc(x,y,squeeze(Y(1,:,:)));
                else
                    imagesc(x,y,Y);
                end
                axis image; colorbar;
                drawnow;
        end
        
    end
end