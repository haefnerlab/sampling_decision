clear all
clc
%close all
%%% This script tests image rotation
%%

ny = 32;
nx = 32;
nL = 1;

%Y=0*Y; warning('zero noise!!!');


c(1,1,1) = 0;
c(1,1,2) = 40;
c(2,1,1) = 40;
c(2,1,2) = 0;
nr = size(c,1);

x=linspace(-1/2,nL-1/2,nx); y=linspace(-1/2,1/2,ny);
[xx, yy] = meshgrid(x,y);

image_task = 'oblique';
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
Y = randn(nr,ny,nx); % white noise background


for i=1:nL
    %%% cardinal imagesc
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
Y = squeeze(Y);
figure
for k = 1:nr
    subplot(nr, 1, k)
    imagesc(squeeze(Y(k,:,:)))

end
