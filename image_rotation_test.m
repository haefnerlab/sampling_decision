clear all
clc
close all
%%% This script tests image rotation
%%
nr = 1;
ny = 32;
nx = 32;
nL = 1;
Y=randn(nr,ny,nx); % white noise background
%Y=0*Y; warning('zero noise!!!');
x=linspace(-1/2,nL-1/2,nx); y=linspace(-1/2,1/2,ny);
[xx yy]=meshgrid(x,y);

for i=1:nL
    %auxV=Gabor_neu([0 1 2 0  i-1 0.1],xx,'orig').*normpdf(yy,  0,0.2);
    auxV=Gabor_neu([0 1 2 0  i-1 0.1],xx,'orig'); % new July 2015
    auxV=auxV/norm(auxV);
    %auxH=Gabor_neu([0 1 2 0    0 0.1],yy,'orig').*normpdf(xx,i-1,0.2);
    auxH=Gabor_neu([0 1 2 0    0 0.1],yy,'orig'); % new July 2015
    auxH=auxH/norm(auxH);
    for k=1:nr
        Y(k,:,:)=squeeze(Y(k,:,:))+c(k,i,1)*auxV+c(k,i,2)*auxH;
    end
end
Y=squeeze(Y);