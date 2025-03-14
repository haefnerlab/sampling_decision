clear all
clc
close all
%%% This script tests image rotation
%% use ralf's code
figure
imSize = 32;
ny = imSize;nx = imSize;nL = 1;
nIm = 1;
contrast_list = [5,15,25,35];
for i_stim = 1:numel(contrast_list)
    contrast = contrast_list(i_stim);
    
    c(1,1,1) = 0;c(1,1,2) = contrast;
    c(2,1,1) = contrast;c(2,1,2) = 0;
    nr = size(c,1);
    x=linspace(-1/2,nL-1/2,nx); y=linspace(-1/2,1/2,ny);
    [xx, yy] = meshgrid(x,y);

    image_task_list = {'cardinal', 'oblique'};
    for itask = 1:numel(image_task_list)
        switch image_task_list{itask}
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
            auxV=Gabor_neu([0 1 2 0  i-1 0.1],xxr,'orig').*normpdf(yyr,  0,0.2);
            %auxV=Gabor_neu([0 1 2 0  i-1 0.1],xxr,'orig'); % new July 2015
            auxV=auxV/norm(auxV);
            auxH=Gabor_neu([0 1 2 0    0 0.1],yyr,'orig').*normpdf(xxr,i-1,0.2);
           % auxH=Gabor_neu([0 1 2 0    0 0.1],yyr,'orig'); % new July 2015
            auxH=auxH/norm(auxH);
            for k=1:nr
                im_signal = c(k,i,1)*auxV+c(k,i,2)*auxH;
                im_signal(xx.^2+yy.^2 > (1/2).^2) = 0;
                Y(k,:,:)=squeeze(Y(k,:,:)) + im_signal;
            end
        end
        Y = squeeze(Y);
    
        for k = 1:nr
            subplot(4, 4, nIm)
            imagesc(squeeze(Y(k,:,:)));
            title(sprintf('Contrast = %d',contrast))
            colormap gray
            colorbar
            nIm = nIm + 1;
        
        end
    end
end
sgtitle("Ralf code")
%% use adam's code
contrast_list = [0.5, 1];
noise_list    = [0, 0.2]; 
orientation_list = [0,90,45,135];
nIm = 1;
figure
for i = 1:numel(contrast_list)
    for j = 1:numel(noise_list)
        for n = 1:numel(orientation_list)
            gabor = make_gabor_acs('isize',imSize,'orient_deg',orientation_list(n),'contrast',contrast_list(i),'noise',noise_list(j));
            
            subplot(4,4,nIm)
            imagesc(gabor)
            colormap gray
            colorbar
            clim([0,1])
            title(sprintf('Contrast = %.1f, noise = %.1f',contrast_list(i), noise_list(j)))
            nIm = nIm + 1;
        end
    end
end
sgtitle("Adam code")
%% a combination between ralf's and adam's
%%% generate white noise then add output of make gabor
contrast_list = [5,15,25,35];
figure
nIm = 1;
for i = 1:numel(contrast_list)
    contrast = contrast_list(i);
    
    c(1,1,1) = 0;c(1,1,2) = contrast;
    c(2,1,1) = contrast;c(2,1,2) = 0;
    nr = size(c,1);
   

    for n = 1:numel(orientation_list)
        Y = randn(ny,nx); % white noise background
        gabor = make_gabor_acs('isize',imSize,'orient_deg',orientation_list(n),'contrast',1,'noise',0);
        Y = Y + gabor .* contrast;
        subplot(4,4,nIm)
        imagesc(Y)
        colormap gray
        colorbar
        %clim([0,1])
        title(sprintf('Contrast = %d',contrast_list(i)))
        nIm = nIm + 1;
    end
  
end
sgtitle("Ralf & adam code")