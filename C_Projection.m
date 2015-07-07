function P = C_Projection(fct,varargin)

DEBUG = 0;
% P.G projection matrix
% P.x x-axis
% P.y y-axis
% P.nx, P.ny size of image
% P.tau tau for prior of X not in agreement with l & o

nx=Get_Varargin(varargin);
if isempty(nx), nx=32; end

switch fct
  case '2x2'
    if nargin<3, P.tau=1; else P.tau=varargin{2}; end
    ny=nx;
    x=linspace(-2,2,nx); y=linspace(-2,2,ny);
    [xx yy]=meshgrid(x,y);
    LV=Gabor_neu([0 1 1 pi/2 -1 0.2],xx,'orig').*normpdf(yy, 0,0.3);
    RV=Gabor_neu([0 1 1 pi/2  1 0.2],xx,'orig').*normpdf(yy, 0,0.3);
    LH=Gabor_neu([0 1 1 pi/2  0 0.2],yy,'orig').*normpdf(xx,-1,0.3);
    RH=Gabor_neu([0 1 1 pi/2  0 0.2],yy,'orig').*normpdf(xx, 1,0.3);
    if DEBUG, debug_RF; end
    G=zeros(nx*ny,4);
    G(:,1)=LV(:)/norm(LV(:)); G(:,2)=LH(:)/norm(LH(:));
    G(:,3)=RV(:)/norm(RV(:)); G(:,4)=RH(:)/norm(RH(:));
    P.nx=nx; P.ny=ny; P.x=x; P.y=y; P.G=G;
    
  case 'nx2' % n locations, 2 orientations, NOT TESTED!
    if nargin<3, P.tau=1; else P.tau=varargin{2}; end
    if nargin<4, P.nL=1; else P.nL=varargin{3}; end
    ny=nx; nx=P.nL*ny;
    x=linspace(-1/2,1/2,nx); y=linspace(-1/2,1/2,ny);
    [xx yy]=meshgrid(x,y);
    G=zeros(nx*ny,2*nL);
    for i=1:nL
      aux=Gabor_neu([0 1 2 pi/2  i-1 0.1],xx,'orig').*normpdf(yy,  0,0.2);
      G(:,2*i-1)=aux(:)/norm(aux);
      aux=Gabor_neu([0 1 2 pi/2    0 0.1],yy,'orig').*normpdf(xx,i-1,0.2);
      G(:,2*i)=aux(:)/norm(aux);
    end
    P.nx=nx; P.ny=ny; P.x=x; P.y=y; P.G=G;
    
  case '1xN'
    if nargin<3, P.tau=1; else P.tau=varargin{2}; end
    if nargin<4, nX=8;  else nX=varargin{3}; end % number of orientations
    %if nargin<4, P.tauD_amp=0.5;  else P.tauD_amp=varargin{3}; end
    %if nargin<5, P.tauO_amp=0.1;  else P.tauO_amp=varargin{4}; end
    ny=nx;
    x=linspace(-1/2,1/2,nx); y=linspace(-1/2,1/2,ny);
    %InputImage anpassen!!
    [xx yy]=meshgrid(x,y);
    phi=(0:nX-1)/nX*pi;
    G=zeros(nx*ny,nX);
    for i=1:nX
      rot=[[cos(phi(i)) sin(phi(i))];[-sin(phi(i)) cos(phi(i))]];
      zza=rot*[xx(:)'; yy(:)']; xxs=zza(1,:); yys=zza(2,:);
      auxG=Gabor_neu([0 1 2 0 0 0.1],xxs,'orig').*normpdf(yys,0,0.2);
      G(:,i)=auxG/norm(auxG);
    end
    P.nx=nx; P.ny=ny; P.x=x; P.y=y; P.G=G; P.nX=nX; P.phi_x=phi;
    
  case 'nxN' % n locations, N orientations
    [P.nX P.nG P.nL]=Get_Varargin;
    if isempty(P.nX), P.nX=8; end % number of orientations
    if isempty(P.nG), P.nG=P.nX/2; end
    if isempty(P.nL), P.nL=1; end
    ny=nx; nx=P.nL*ny; % size of the image (Y)
    x=linspace(-1/2,P.nL-1/2,nx); y=linspace(-1/2,1/2,ny);
    phi_x=(0:P.nX-1)/P.nX*pi;
    phi_g=(0:P.nG-1)/P.nG*pi;
    G=zeros(nx*ny,P.nX*P.nL); % projection matrix
    for j=1:P.nL
      cx=j-1; % center
      [xx yy]=meshgrid(x-cx,y);
      for i=1:P.nX
        rot=[[cos(phi_x(i)) sin(phi_x(i))];[-sin(phi_x(i)) cos(phi_x(i))]];
        zza=rot*[xx(:)'; yy(:)'];
        xxs=zza(1,:); yys=zza(2,:);
        auxG=Gabor_neu([0 1 2 0 0 0.1],xxs,'orig').*normpdf(yys,0,0.2);
        G(:,(j-1)*P.nX+i)=auxG/norm(auxG); % division by norm!
      end
    end
    P.nx=nx; P.ny=ny; P.x=x; P.y=y; P.G=G; 
    P.phi_x=phi_x; P.phi_g=phi_g;
    
  otherwise
    error(fct);
end

if DEBUG, debug_RF; end


  function debug_RF
    Get_Figure(['CPRF:' fct]);
    clim=[-0.8 0.8];
    switch fct
      case '2x2'
        Subplot(4,1,2,2); imagesc(x,y,LV,clim); axis image; colorbar; title('LV 1');
        Subplot; imagesc(x,y,RV,clim); axis image; colorbar; title('RV 3');
        Subplot; imagesc(x,y,LH,clim); axis image; colorbar; title('LH 2');
        Subplot; imagesc(x,y,RH,clim); axis image; colorbar; title('RH 4');
      case '1xN'
        Subplot(nX,1,nX,1);
        for i=1:nX, Subplot(i);
          imagesc(reshape(G(:,i),nx,ny)); axis image;
        end
      case 'nxN'
        Subplot(P.nL*P.nX,0,P.nX,P.nL);
        for j=1:P.nL
          for i=1:P.nX, Subplot;
            imagesc(reshape(G(:,(j-1)*P.nX+i),ny,nx)); axis image;
          end
        end
    end
  end

end