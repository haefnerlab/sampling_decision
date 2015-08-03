%calculating correlation matrix, eigenvalue and eigenvectors
%created by Shuchen Wu 
%06/2015
function experiment_PCA(e)
close all;

%load e_Detection_T0_1_X1024_G256_k10_3_t80_c0_rep1000_time100_b5_s001_001_i20.mat
X = e.X;
X(isnan(X))=0;

%%compute Corr_Matrix
%%each row specify neuron(repetition) and each 100 column represent time in
%%a trail
%% should sum over all the samples in one trail for each neuron
%% compute the covariance matrix between total responses in each trail
%%(trail * #of neuron *  time)
close all;
figure
subplot(2,2,1);
X_reshape = sum(X,3); %%sum the time trail of that data,1000*1024 array
Corr_Matrix = corr(X_reshape);
Corr_Matrix(isnan(Corr_Matrix)) =0;
X_range = e.Projection.phi_x;%obtain axis range
imagesc(X_range,X_range,Corr_Matrix);
colorbar;
caxis auto;
%%why is x and y labelled deltaphi?

xlabel('\Delta\phi'); 
ylabel('\Delta\phi');
set(gca,'xtick',[0 pi],'xticklabel',{'0','pi'},...
       'ytick',[0 pi],'yticklabel',{'0','pi'});
set(gca,'Ydir','Normal');
axis tight;
caxis([0,0.4]);



X1_reshape = sum(X,3);%%sum the time trail of that data,1000*1024 array
Corr_Matrix_1 = corr(X1_reshape);

X_range_1 = e.Projection.phi_x;%obtain axis range

imagesc(X_range_1,X_range_1,Corr_Matrix_1);
colorbar;
caxis auto;
%%why is x and y labelled deltaphi?

xlabel('\Delta\phi'); 
ylabel('\Delta\phi');
set(gca,'xtick',[0 pi],'xticklabel',{'0','pi'},...
       'ytick',[0 pi],'yticklabel',{'0','pi'});
set(gca,'Ydir','Normal');
axis tight;
caxis([0,0.4]);

%% eigenvector, eigenvalue
[V, D] = eig(Corr_Matrix);

Avg_response_var = mean(e.X,2);
size(Avg_response_var);
 %%plotting the average response of all neurons over time for the first
 %%trail
%  
% plot(squeeze(Avg_response_var(1,1,:)));
% title('Average resoponse of all neurons given the first stimulus over time')
% xlabel('time');
% ylabel('average spikes')
% figure;
% 
% plot(real(eig(Corr_Matrix)), '*')%%unnormalized eigenvalue
% title('plotting the unnormalized eigenvalue')
% figure;

%% the biggest eigenvectors and eigenvalues
%[V_max, D_max] = eigs(Corr_Matrix);

% figure;
% plot(real(Corr_Matrix(:,1)),real(Corr_Matrix(:,2)),'b.');
% title('distribution of correlation matrix')
% hold on;
% plot(3*[-V_max(2,1) V_max(2,1)],3*[-V_max(2,2) V_max(2,2)],'k')
% plot(3*[-V_max(1,1) V_max(1,1)],3*[-V_max(1,2) V_max(1,2)],'k')
    % use the following statement to retrieve X for ICA or so
    %varargout{1}=sum(X(:,:,n0S:end),3); return;
   
    
%%Compute CP
O = e.O;
TPR =  O(:,2,end)>0.5;%True positive rate drawn at the end of the trails?
FPR = O(:,3,end)>0.5;%false positive rate 
O_reshape = sum(O,3);%sum obver time domain
%  trying to construct a ROC curve 


% for k = 1:100
%   TPR =  O(:,2,k)>0.5;%True positive rate drawn at the end of the trails?
%   FPR = O(:,3,k)>0.5;
%   plot(length(O(TPR))/1000,length(O(FPR))/1000,'*')
%   hold on;
% end

[trails, neurons, time] = size(X);

CP = zeros(neurons,time);%%cp of neurons across time
%minlength = min(length(pref),length(anti))
 
for i = 1: neurons
  for j = 1 : time
    pref=X(TPR,i,j);
    anti=X(FPR,i,j);
    L_pref = length(pref);
    L_anti = length(anti);
%     nPref = 1:L_pref;
%     nAnti = 1:L_anti;
    
    minlength = min(L_pref,L_anti);

    %CP(i,j) = ROC(pref(nPref),anti(nAnti));
    %CP(i,j)=(sum(x<y)+0.5*sum(x==y))/length(x)
    CP(i,j) = (sum(pref(1:minlength)<anti(1:minlength))+0.5*sum(pref(1:minlength)==anti(1:minlength)))/minlength;

  end
end
CP = mean(CP,2);% average Cp over time
 subplot(2,2,2);
plot(e.Projection.phi_x, CP,'b');
title('simulated CP')
xlabel('orientation')
ylabel('CP')

%Corr_Inv = inv(Corr_Matrix);
%Weights = zeros(neurons,1); % construct 1024*1 matrix for weights 
Const =  pi/sqrt(2);
 subplot(2,2,3);

C = Corr_Matrix;
Gamma = (CP - 0.5).*Const;%%Sqrt(Ckk) is always 1, so doesnt matter
Weights = Gamma*1024\C;%Gamma = C*Weights
plot(e.Projection.phi_x, Weights,'r*');
title('modelled CP')
xlabel('orientation')
ylabel('CP')



% for k = 1: neurons
%   Sum = 0;
%   for l = 1: neurons
%     Sum = Sum + Corr_Inv(k,l)*sqrt(Corr_Matrix(l,l))*(CP(l)-0.5);
%   end
%    Weights(k,1) = Sum;
% end 
% Weights = Weights.*Const./(neurons);
% figure
% plot(e.Projection.phi_x,Weights(:,1), 'r*')
Gamma_re =  C*Weights';
CP_re = Gamma_re./Const+0.5;
subplot(2,2,4);

plot(e.Projection.phi_x,CP_re,'r*')
hold on
plot(e.Projection.phi_x,CP,'b*')
title('comparison between simulated CP and modelled CP')
xlabel('orientation')
ylabel('CP')



    
 
    
    


