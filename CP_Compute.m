function CP_Compute(e)
% CP_COMPUTE compute and plot choice probabilities from results of
% S_Run_Experiments
%
%       CP_COMPUTE(e) where 'e' is the output of S_Run_Experiments and
%       contains field X, a (Trials x Neurons x Time) matrix of spike rates

figure;

nconditions = length(e); % e is a cell array of structs, one per condition

% compute <|CP(0.5) - 0.5|> using the unbiased stimulus
% identify TPR and FPR for that particular struct
O = e.O;
X = e.X; % load X.

TPR = O(:,2,end) > 0.5; % True positive rate drawn at the end of the trails
FPR = O(:,3,end) > 0.5; % false positive rate

[trials, neurons, time] = size(X);

CP_half = zeros(neurons,time); % cp of neurons across time
% minlength = min(length(pref),length(anti))
choiceRatio_half = zeros(neurons,time);

%% Loop through neurons to generate half CP and half choice ratio
for i = 1: neurons
    for j = 1 : time
        pref = X(TPR,i,j);
        anti = X(FPR,i,j);
        L_pref = length(pref);
        L_anti = length(anti);
        minlength = min(L_pref,L_anti);
        % when length of pref and anti are different, roc function should be used
        %   CP(i,j) = ROC(pref(nPref),anti(nAnti));
        %   CP(i,j)=(sum(x<y)+0.5*sum(x == y))/length(x)
        % because of the length of pref and anti are the same, this method, instead of the ROC curve, can be used
        CP_half(i,j) = (sum(pref(1:minlength)<anti(1:minlength))+0.5*sum(pref(1:minlength)== anti(1:minlength)))/minlength;
        choiceRatio_half(i,j) = L_pref/(L_pref+L_anti);
    end
end

% plot <|CP(0.5) - 0.5|>
CP_half = mean(CP_half,2); % average CP_half over all time bin
choiceRatio_half = mean(choiceRatio_half,2);
subplot(2,2,1);
denom = abs(CP_half-0.5);
plot(choiceRatio_half, denom,'*');
title(' <|CP(0.5) - 0.5|> ');
% ratio_sum = zeros(1,structsize); % ratiosum sums all rations across different stimulus

%% 2 matrix to obtain all choice ratio, averages, and error value across each loop
all_choice_ratio = zeros(1, nconditions);
all_average = zeros(1,nconditions);
all_errorbar = zeros(1,nconditions);

ratio_matrix_x = zeros(1,nconditions*neurons);
ratio_matrix_y = zeros(1,nconditions*neurons);

for m = 1:nconditions
    O = e.O;
    TPR = O(:,2,end)>0.5; % True positive rate drawn at the end of the trails
    FPR = O(:,3,end)>0.5; % false positive rate
    X = e.X; % load X.
    
    % optional plot of ROC curve
    %
    % for k = 1:80 % for each time bin, calculate TPR & FPR for the time bin
    %   TPR =  O(:,2,k)>0.5; % True positive rate drawn at the end of the trails?
    %   FPR = O(:,3,k)>0.5;
    %    plot(length(O(TPR))/1000,length(O(FPR))/1000,'*')
    %   hold on;
    % end
    
    
    [trials, neurons, time] = size(X);
    
    CP = zeros(neurons,time); %% cp of neurons across time
    % minlength = min(length(pref),length(anti))
    choiceRatio = zeros(neurons,time);
    
    % in struct 7, has a TPR of 0
    for i = 1: neurons
        for j = 1 : time
            pref = X(TPR,i,j);
            anti = X(FPR,i,j);
            L_pref = length(pref);
            L_anti = length(anti);
            %     nPref = 1:L_pref;
            %     nAnti = 1:L_anti;
            minlength = min(L_pref,L_anti);
            % CP(i,j) = ROC(pref(nPref),anti(nAnti));
            % CP(i,j)=(sum(x<y)+0.5*sum(x == y))/length(x)
            CP(i,j) = (sum(pref(1:minlength)<anti(1:minlength))+0.5*sum(pref(1:minlength)== anti(1:minlength)))/minlength;
            choiceRatio(i,j) = L_pref/(L_pref+L_anti);
        end
    end
    
    CP = mean(CP,2); % average Cp over time
    choiceRatio = mean(choiceRatio,2);
    nom = abs(CP-0.5);
    % plot(e.Projection.phi_x, CP,'b');
    subplot(2,2,2);
    plot(choiceRatio, nom,'*');
    hold on
    title('<|CP - 0.5|>')
    ratio = nom./denom; %% num and denom are two arrays with the length of the number of neurons
    
    % Standard error of the <CP-.5>/<CP(0.5)-0.5>
    SE = std(ratio)/sqrt(length(ratio));
    
    % % plot all data distribution
    % subplot(2,2,3);
    % plot(choiceRatio, ratio ,'*');
    % title('<|CP - 0.5|>/<|CP(0.5) - 0.5|> ')
    % hold on
    
    
    %% load data and ratio into the grand ratio matrix to plot later and divide
    %% into groups
    start = (m-1)*256+1;
    End = m*256;
    ratio_matrix_x(start:End) = choiceRatio(1,1); %% since choiceratio doesnt change across trails
    ratioT = ratio';
    ratio_matrix_y(start:End) = ratioT(1:256)';
    
    % calculate mean
    ratio_sum = sum(ratio);
    average = ratio_sum./neurons;
    
    
    % plot error bar with mean
    subplot(2,2,4)
    plot(choiceRatio,average,'b')
    choiceRatio = choiceRatio(1); % choice ratio does not change across all neurons
    all_choice_ratio(1, m) = choiceRatio;
    all_average(1,m) = average;
    all_errorbar(1,m) = SE;
    
    errorbar(choiceRatio,average,SE) %% choice ratio, average are 1*1
    title('average choice ratio');
    
end


%% divide into different sized bins to determine their median and error bar
%% at  25 % bound
binumber = 2; % use the number of bins to determine binsize.
binwidth = 1/binumber; % given 1 is the largest ratio possible staring at 0
binmatrix = zeros(2,length(ratio_matrix_x)); % binmatrix = 2*lengthofdata
binmatrix(1,:) = ratio_matrix_x;
binmatrix(2,:) = ratio_matrix_y; %% initials all binmatrix to be nan
sorted_matrix = sortrows(binmatrix'); %% sort the data by X bin
sorted_matrix = sorted_matrix'; % transverse sorted matrix 2*1792(all number of data)
% declare matrix to hold data in each bin; 3Dmatrix with
% 1*length(ratio_matrix_x)*binnumber
binned_matrix = zeros(2,length(ratio_matrix_x),binumber);
binned_matrix(:,:,:) = nan;

startbin = 1;
endbin = 1+binwidth;


for k = 1:(binumber)
    %% sort matrix into different bins
    startbin = startbin - binwidth;
    endbin = endbin-binwidth; % decrement endbin
    
    A = sorted_matrix(:,intersect(find(sorted_matrix(1,:)>= startbin),find(sorted_matrix(1,:)<= endbin)));
    binned_matrix(:,1:length(A),k) = A;
    % binumber = binumber - 1;
end

subplot(2,2,4)
CFBD_matrix = zeros(2,3,binumber); %% matrix to contain median and 2 confidence bounds
% for xaxis and y axis with each numbered bins
%% confidence interval
q = binned_matrix(~isnan(binned_matrix)); %% excluding all the nans from binned matrix
CFBD_matrix = prctile(binned_matrix,[25 50 75],2); %% matrix dimension = 2*3*binumber

% plot(CFBD_matrix)
ERU_y = squeeze(CFBD_matrix(2,3,:)) - squeeze(CFBD_matrix(2,2,:)); %% upper error limit
ERL_y = squeeze(CFBD_matrix(2,2,:)) - squeeze(CFBD_matrix(2,1,:)); %% lower error limit

ERU_x = squeeze(CFBD_matrix(1,3,:)) - squeeze(CFBD_matrix(1,2,:)); %% upper error limit
ERL_x = squeeze(CFBD_matrix(1,2,:)) - squeeze(CFBD_matrix(1,1,:)); %% upper error limit

plot(squeeze(CFBD_matrix(1,2,:)),squeeze(CFBD_matrix(2,2,:)),'*');
hold on;
errorbar(CFBD_matrix(1,2,:),CFBD_matrix(2,2,:),ERU_y,ERL_y,'*');
hold on
% herrorbar(CFBD_matrix(1,2,:),CFBD_matrix(2,2,:),ERU_x, ERL_x, '*')
% hold on
% plot(squeeze(CFBD_matrix(1,1,:)),squeeze(CFBD_matrix(1,2,:)),'r*')
%
% hold on
% plot(squeeze(CFBD_matrix(1,3,:)),squeeze(CFBD_matrix(1,2,:)),'r*')

% % plot summary of previous data points, with mean and error bar in one line
% % plot
% subplot(2,2,4)
% errorbar(all_choice_ratio,all_average,all_errorbar) %% choice ratio, average are 1*1
% title('average choice ratio');
% % plot(choiceRatio, ratio(1,:) ,'*');
%
% %% plot the scattered plot according to different size of bin histogram
%
% % 1. collect all the points from all rounds of simulation
%
% plot all data distribution
subplot(2,2,3);
plot(ratio_matrix_x, ratio_matrix_y ,'*');
title('<|CP - 0.5|>/<|CP(0.5) - 0.5|> ')
hold on
end