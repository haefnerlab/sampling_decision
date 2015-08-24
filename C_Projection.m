function P = C_Projection(fct, varargin)
%C_PROJECTION creates gabor projective fields
%
%   P = C_PROJECTION(fct, ...) returns a struct with the following fields:
%       P.G projection matrices as (flattened 1D pixels)x(n_matrices)
%       P.x x-axis
%       P.y y-axis
%       P.nx, P.ny size of each projective field 
%                  (i.e. P.nx*P.ny == size(P.G(:,1))
%       P.tau controls prior of X not in agreement with l & o (not used in
%             nxN case)
%
%       fct is a string with '#locations x #orientations' and may be '2x2',
%       'nx2', '1xN', or 'nxN'
%
%   C_PROJECTION(fct, ..., 'debug') visualize the resulting projective fields
%
%   P = C_Projection('2x2', [nx, [tau]]) Creates two gabors in two locations
%       (Left, Right) x (Vertical, Horizontal). nx defaults to 32, tau to 1
%
%   P = C_Projection('nx2', [nx, [tau, [nL]]) Creates gabors in nL locations,
%       where each gabor is given a square patch. The resulting projective
%       fields are each (nx*nL) by (nx) pixels. nL defaults to 1
%       The field P.nL is also set.
%
%   P = C_Projection('1xN', [nx, [tau, [nX]]]) Creates gabors at nX orientations
%       The fields P.nX, P.phi_x are also set.
%
%   P = C_Projection('nxN', [nx, [nX, [nG, [nL]]]]]) Creates gabores at nL
%       locations, with nX orientations. (unclear what nG does)
%       P.nX, P.nG, P.nL, P.phi_x, P.phi_g all set

% find the string 'debug' anywhere in varargin then remove it so it doesn't
% interfere with any other varargin parameters
nargin_ = nargin;
dbstr = strcmp('debug', varargin);
varargin = varargin(~dbstr);
DEBUG = any(dbstr);
% decrement nargin_ if DEBUG is set (which is why we can't use the regular
% nargin here)
nargin_ = nargin_ - DEBUG;

if nargin_ < 2
    P.nx = 32;
else
    P.nx = varargin{1};
end

switch fct
    case '2x2'
        if nargin_ < 3
            P.tau = 1;
        else
            P.tau = varargin{2};
        end
        P.ny = P.nx;
        % create x, y spatial coordinates
        P.x = linspace(-2,2,P.nx);
        P.y = linspace(-2,2,P.ny);
        [xx, yy] = meshgrid(P.x,P.y);
        % create gabors (Left, Right)x(Vertical, Horizontal)
        LV = Gabor_neu([0 1 1 pi/2 -1 0.2], xx, 'orig') .* normpdf(yy, 0,0.3);
        RV = Gabor_neu([0 1 1 pi/2  1 0.2], xx, 'orig') .* normpdf(yy, 0,0.3);
        LH = Gabor_neu([0 1 1 pi/2  0 0.2], yy, 'orig') .* normpdf(xx,-1,0.3);
        RH = Gabor_neu([0 1 1 pi/2  0 0.2], yy, 'orig') .* normpdf(xx, 1,0.3);
        % create flattened matrix to store each projective field as a column vector
        P.G = zeros(P.nx*P.ny, 4);
        P.G(:,1) = LV(:) / norm(LV(:));
        P.G(:,2) = LH(:) / norm(LH(:));
        P.G(:,3) = RV(:) / norm(RV(:));
        P.G(:,4) = RH(:) / norm(RH(:));
        
    case 'nx2' % n locations, 2 orientations, NOT TESTED!
        if nargin_ < 3,
            P.tau = 1;
        else
            P.tau = varargin{2};
        end
        if nargin_ < 4,
            P.nL = 1;
        else
            P.nL = varargin{3};
        end
        % height=width to make each patch square
        P.ny = P.nx;
        % one square patch per location, so total width is nx * nL
        P.nx = P.nL * P.nx;
        P.x = linspace(-1/2, 1/2, P.nx); % TODO - should this be linspace(-P.nL/2, P.nL/2, P.nx*P.nL) ??
        P.y = linspace(-1/2, 1/2, P.ny);
        [xx, yy] = meshgrid(P.x, P.y);
        % create flattened matrix to store each projective field as a column vector
        P.G = zeros(P.nx*P.ny, 2 * P.nL);
        for i=1:P.nL
            % one horizontal and one vertical gabor per patch
            % TODO - these don't look right
            gaborV = Gabor_neu([0 1 2 pi/2 i-1 0.1], xx, 'orig') .* normpdf(yy, 0, 0.2);
            P.G(:,2*i-1) = gaborV(:) / norm(gaborV);
            gaborH = Gabor_neu([0 1 2 pi/2 0 0.1], yy, 'orig') .* normpdf(xx, i-1, 0.2);
            P.G(:,2*i)=gaborH(:)/norm(gaborH);
        end
        
    case '1xN'
        if nargin_ < 3
            P.tau = 1;
        else
            P.tau = varargin{2};
        end
        % 'nX' is number of orientations
        if nargin_ < 4
            P.nX = 8;
        else
            P.nX = varargin{3};
        end 
        % set up square image space
        P.ny = P.nx;
        P.x = linspace(-1/2, 1/2, P.nx);
        P.y = linspace(-1/2, 1/2, P.ny);
        [xx, yy] = meshgrid(P.x, P.y);
        % uniformly distribute orientations
        P.phi_x = (0:P.nX-1) / P.nX * pi;
        % create flattened matrix to store each projective field as a column vector
        P.G = zeros(P.nx*P.ny,P.nX);
        % create gabor at each orientation
        for i=1:P.nX
            % rotation matrix transform
            c = cos(P.phi_x(i));
            s = sin(P.phi_x(i));
            rot = [[c s]; [-s c]];
            % rotate the pixel coordinates we got from meshgrid
            zza = rot * [xx(:)'; yy(:)'];
            xxs = zza(1,:);
            yys = zza(2,:);
            % get Gabor on rotated coordinates
            gabor = Gabor_neu([0 1 2 0 0 0.1], xxs, 'orig') .* normpdf(yys,0,0.2);
            P.G(:,i) = gabor / norm(gabor);
        end
        
    case 'nxN'
        % 'nX' is number of orientations
        if nargin_ < 3
            P.nX = 8;
        else
            P.nX = varargin{2};
        end
        if nargin_ < 4
            P.nG = P.nX / 2;
        else
            P.nG = varargin{3};
        end
        if nargin_ < 5
            P.nL = 1;
        else
            P.nL = varargin{4};
        end
        % set up projective field sizes
        P.ny = P.nx;
        P.nx = P.nL * P.ny;
        P.x = linspace(-1/2, P.nL-1/2, P.nx); % TODO - make this and P.x definition in 'nx2' the same ?
        P.y = linspace(-1/2, 1/2, P.ny);
        P.phi_x = (0:P.nX-1) / P.nX * pi;
        P.phi_g = (0:P.nG-1) / P.nG * pi;
        % create flattened matrix to store each image as a column vector
        P.G = zeros(P.nx*P.ny,P.nX*P.nL);
        % loop over locations
        for j=1:P.nL
            cx = j-1; % center x location (see definition of P.x above)
            [xx, yy]=meshgrid(P.x-cx,P.y);
            % loop over orientations
            for i=1:P.nX
                % rotation matrix transform
                c = cos(P.phi_x(i));
                s = sin(P.phi_x(i));
                rot = [[c s]; [-s c]];
                % rotate the pixel coordinates we got from meshgrid
                zza = rot * [xx(:)'; yy(:)'];
                xxs = zza(1,:);
                yys = zza(2,:);
                % get Gabor on rotated coordinates
                gabor = Gabor_neu([0 1 2 0 0 0.1], xxs, 'orig') .* normpdf(yys,0,0.2);
                P.G(:,(j-1)*P.nX+i) = gabor / norm(gabor);
            end
        end
        
    otherwise
        error(fct);
end

if DEBUG
    debug_RF(fct, P);
end

end

function debug_RF(fct, P)
Get_Figure(['CPRF:' fct]);
clim = [-0.8 0.8];
switch fct
    case '2x2'
        Subplot(4,1,2,2); imagesc(P.x,P.y,LV,clim); axis image; colorbar; title('LV 1');
        Subplot; imagesc(P.x,P.y,RV,clim); axis image; colorbar; title('RV 3');
        Subplot; imagesc(P.x,P.y,LH,clim); axis image; colorbar; title('LH 2');
        Subplot; imagesc(P.x,P.y,RH,clim); axis image; colorbar; title('RH 4');
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