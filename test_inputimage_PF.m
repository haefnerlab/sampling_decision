clear all
clc
close all
%%
% Parameters
isize = 32;         % Size of the patch (pixels)
sigma = 0.2;         % Gaussian envelope standard deviation
f = 2;               % Spatial frequency 

phi = 0;            % Phase offset
gamma = 1;          % Aspect ratio
A = 1;              % Amplitude (contrast)

figure
theta_list = [0,90,45,135];
[orientation_bins, orientation_energy, gabor_image] = deal(cell(numel(theta_list), 1));

for n = 1:numel(theta_list)
    theta = theta_list(n);         % Orientation angle in degrees
    % Create a grid
    x=linspace(-1/2,1/2,isize); y=linspace(-1/2,1/2,isize);
    [xx, yy]=meshgrid(x,y);
    
    % Convert theta to radians
    theta_rad = deg2rad(theta);
    
    % Rotate the grid
    x_theta = xx * cos(theta_rad) + yy * sin(theta_rad);
    y_theta = -xx * sin(theta_rad) + yy * cos(theta_rad);
    
    % % Gabor function
    % gabor_1 = A * exp(- (x_theta.^2 + y_theta.^2) / (2 * sigma^2)) ...
    %             .* cos(2 * pi * f * x_theta + phi);
    % gabor_2 = A * exp(- (y_theta.^2 +  x_theta.^2) / (2 * sigma^2)) ...
    %             .* cos(2 * pi * f * y_theta + phi);

    gabor_1=Gabor_neu([0 1 f 0  0 0.1], x_theta, 'orig').* normpdf(y_theta,0,sigma);
    %gabor_1((x_theta.^2+y_theta.^2) > (1/2).^2) = 0;
    norm(gabor_1)
    norm(gabor_1(:))
    gabor_1 = gabor_1 / norm(gabor_1(:));
    
    gabor_image{n} = gabor_1;
    subplot(1,4,n);imagesc(gabor_1); %clim([-A,A])
   % [orientation_bins{n}, orientation_energy{n}] = get_ori_energy_model(gabor_1);
    
end
% figure
% for n  = 1:4
%     plot(orientation_bins{n}, orientation_energy{n}); hold on
% end

%% simulate projective field
P.nX = 64;
P.nx = isize;P.ny = isize;
P.x = linspace(-1/2, 1, P.nx); % TODO - make this and P.x definition in 'nx2' the same ?
P.y = linspace(-1/2, 1/2, P.ny);
P.phi_x = (0:P.nX-1) / P.nX * pi;
P.G     = zeros(P.nx*P.ny,P.nX);
for i=1:P.nX
    
    % rotation matrix transform
    c = cos(P.phi_x(i));
    s = sin(P.phi_x(i));
    rot = [[c s]; [-s c]];
    % rotate the pixel coordinates we got from meshgrid
    zza = rot * [xx(:)'; yy(:)'];
    xxs = zza(1,:);
    yys = zza(2,:);
    gabor = Gabor_neu([0 1 f 0 0 sigma], xxs, 'orig').* normpdf(yys,0,sigma);
    %gabor((xxs.^2+yys.^2) > (1/2).^2) = 0;
    gabor = gabor/norm(gabor);

    P.G(:,i) = gabor;
    P.G_norm(i) = norm(gabor);
end
%% 
figure
%%%%% cardinal signal
response = zeros(numel(theta_list), P.nX);
for i = 1:numel(theta_list)
    response(n,:) = arrayfun(@(n)sum(reshape(gabor_image{i}, size(P.G(:,1)))  .* P.G(:,n)), [1:P.nX]);

    plot(P.phi_x, response(n,:)); hold on
end

            
%%

