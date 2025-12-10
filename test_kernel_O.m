clear all
clc
close all
%% parameter setting
G.nT = 2;
G.number_orientations = 2;
G.number_locations = 1;
G.dimension_G = 64;

G.kappa_O = [0.9 0.1]; % attended and unattended
G.kappa_G = 3;
G.pO = ones(1,G.number_orientations)/G.number_orientations; % orientations prior (uniform)
G.pL = ones(1,G.number_locations)/G.number_locations; % locations prior (uniform)
G.phi_O = [(0:2:2*(G.number_orientations-1));...
             (1:2:2*G.number_orientations-1)] * pi/2/G.number_orientations;
dimG = G.dimension_G;
G.phi_g = (0:dimG-1) / dimG * pi;

%% single T (original code)
% dimension 1 is [attended, unattended], the rest are self-explanatory
kernel_O = zeros([2 G.nT, G.number_orientations dimG]);
% normalization such that sum over g of p(g | O) is 1 (on average one g is
% 'on' at a time)
norm_kernel_O = 1 ./ besseli(0, G.kappa_O) / dimG;
for T = 1:G.nT
    for O = 1:G.number_orientations
       
            % same for 'discrimination' and 'detection' tasks here
            delta_cos = cos(2 * (G.phi_O(T,O) - G.phi_g));
       
        
        % the difference between attended and unattended is in the width of
        % the von Mises, controlled by kappa_O (see S_Exp_Para)
        kernel_O(1,T,O,:) = norm_kernel_O(1) * exp(G.kappa_O(1)*delta_cos);
        kernel_O(2,T,O,:) = norm_kernel_O(2) * exp(G.kappa_O(2)*delta_cos);
    end
end
kernel_O(kernel_O > 1) = 1;

figure; hold on

plot(squeeze(kernel_O(1,1,1,:))); 
plot(squeeze(kernel_O(1,1,2,:)));
plot(squeeze(kernel_O(1,2,1,:))); 
plot(squeeze(kernel_O(1,2,2,:)));
ylim([0,0.05])

legend('T = cardinal, O = 1',...
       'T = cardinal, O = 2',...
       'T = oblique, O = 1',...
       'T = oblique, O = 2')

title('Original kernel_O','Interpreter','none')
%% two T, define two kernels
kernel_O_cardinal = zeros([2, 1, G.number_orientations, dimG]);
kernel_O_oblique  = zeros([2, 1, G.number_orientations, dimG]);
 % normalization such that sum over g of p(g | O) is 1 (on average one g is
% 'on' at a time)
norm_kernel_O = 1 ./ besseli(0, G.kappa_O) / dimG;
for O = 1:G.number_orientations
    delta_cos = cos(2 * (G.phi_O(1,O) - G.phi_g)); % the first row in phi_O corresponds to the cardinal task
        
    kernel_O_cardinal(1,1,O,:) = norm_kernel_O(1) * exp(G.kappa_O(1)*delta_cos);
    kernel_O_cardinal(2,1,O,:) = norm_kernel_O(2) * exp(G.kappa_O(2)*delta_cos);
end

for O = 1:G.number_orientations
    delta_cos = cos(2 * (G.phi_O(2,O) - G.phi_g)); % the second row in phi_O corresponds to the oblique task
        
    kernel_O_oblique(1,1,O,:) = norm_kernel_O(1) * exp(G.kappa_O(1)*delta_cos);
    kernel_O_oblique(2,1,O,:) = norm_kernel_O(2) * exp(G.kappa_O(2)*delta_cos);
end
kernel_O_uniform  = ones([2, 1, G.number_orientations, dimG]);
%%% define a uniform kernel (no task is on)
kernel_O_uniform             = kernel_O_uniform ./ besseli(0, 0) / dimG/2;

figure
hold on
plot(squeeze(kernel_O_cardinal(1,1,1,:))); 
plot(squeeze(kernel_O_cardinal(1,1,2,:)));
plot(squeeze(kernel_O_oblique(1,1,1,:))); 
plot(squeeze(kernel_O_oblique(1,1,2,:)));
ylim([0,0.05])
%%
ori = G.phi_g / pi * 180;
figure
T_cardinal = 1; T_oblique = 0; O = 1;
option = 'use-T';
posterior_cardinal = 0.5; posterior_oblique = 0.5;
kernel_O_use = get_kernel_O_use(kernel_O_cardinal,kernel_O_oblique, T_cardinal, T_oblique, O, posterior_cardinal, posterior_oblique,option);
plot(ori, squeeze(kernel_O_use),'Linewidth',2); hold on


T_cardinal = 0; T_oblique = 1; O = 1;
option = 'use-T';
posterior_cardinal = 0.5; posterior_oblique = 0.5;
kernel_O_use = get_kernel_O_use(kernel_O_cardinal,kernel_O_oblique, T_cardinal, T_oblique, O, posterior_cardinal, posterior_oblique,option);
plot(ori, squeeze(kernel_O_use),'Linewidth',2); hold on

T_cardinal = 0; T_oblique = 0; O = 1;
option = 'use-T';
posterior_cardinal = 0.5; posterior_oblique = 0.5;
kernel_O_use = get_kernel_O_use(kernel_O_cardinal,kernel_O_oblique, T_cardinal, T_oblique, O, posterior_cardinal, posterior_oblique,option);
plot(ori, squeeze(kernel_O_use),'Linewidth',2); hold on

T_cardinal = 1; T_oblique = 1; O = 1;
option = 'use-T';
posterior_cardinal = 0.5; posterior_oblique = 0.5;
kernel_O_use = get_kernel_O_use(kernel_O_cardinal,kernel_O_oblique, T_cardinal, T_oblique, O, posterior_cardinal, posterior_oblique,option);
plot(ori, squeeze(kernel_O_use),'Linewidth',2); hold on


set(gca,'fontsize',18);
xlabel('Orientation'); ylabel('Probability');
title('Kernel O, O = 1')
legend('T_{cardinal} = 1, T_{oblique} = 0',...
        'T_{cardinal} = 0, T_{oblique} = 1',...
         'T_{cardinal} = 0, T_{oblique} = 0',...
         'T_{cardinal} = 1, T_{oblique} = 1')
%%
ori = G.phi_g / pi * 180;
figure
T_cardinal = 1; T_oblique = 0; 
kernel_O_use = get_kernel_O_use_dual(kernel_O_cardinal, kernel_O_oblique,kernel_O_uniform, T_cardinal, T_oblique);
plot(ori, squeeze(kernel_O_use(1,1,1,:)),'Linewidth',2); hold on
sum(kernel_O_use(1,1,1,:))

T_cardinal = 0; T_oblique = 1; 
kernel_O_use = get_kernel_O_use_dual(kernel_O_cardinal, kernel_O_oblique,kernel_O_uniform, T_cardinal, T_oblique);
plot(ori, squeeze(kernel_O_use(1,1,1,:)),'Linewidth',2); hold on
sum(kernel_O_use(1,1,1,:))

T_cardinal = 0; T_oblique = 0; 
kernel_O_use = get_kernel_O_use_dual(kernel_O_cardinal, kernel_O_oblique,kernel_O_uniform, T_cardinal, T_oblique);
plot(ori, squeeze(kernel_O_use(1,1,1,:)),'Linewidth',2); hold on
sum(kernel_O_use(1,1,1,:))

T_cardinal = 1; T_oblique = 1; 
kernel_O_use = get_kernel_O_use_dual(kernel_O_cardinal, kernel_O_oblique,kernel_O_uniform, T_cardinal, T_oblique);
plot(ori, squeeze(kernel_O_use(1,1,1,:)),'Linewidth',2); hold on
sum(kernel_O_use(1,1,1,:))

set(gca,'fontsize',18);
xlabel('Orientation'); ylabel('Probability');
title('Kernel O, O = 1')
legend('T_{cardinal} = 1, T_{oblique} = 0',...
        'T_{cardinal} = 0, T_{oblique} = 1',...
         'T_{cardinal} = 0, T_{oblique} = 0',...
         'T_{cardinal} = 1, T_{oblique} = 1')
%% winner take all


function kernel_O_use = get_kernel_O_use(kernel_O_cardinal,kernel_O_oblique, T_cardinal, T_oblique, O, posterior_cardinal, posterior_oblique,option)
attn = 1;
dimG = size(kernel_O_cardinal,4);
switch option
    case 'weighted-avg'
        %%%% 1. just use weighted average based on posteriors, regardless T_cardinal and T_oblique
        kernel_O_use = (kernel_O_cardinal(attn,1,O,:) * posterior_cardinal + kernel_O_oblique(attn,1,O,:) * posterior_oblique) / (posterior_cardinal + posterior_oblique);
    case 'winner-take-all'
        %%%% 2. Use posterior in a winner-take-all manner
        if posterior_cardinal > posterior_oblique
            kernel_O_use = kernel_O_cardinal(attn,1,O,:);
        else
            kernel_O_use = kernel_O_oblique(attn,1,O,:);
        end
    case 'use-T'
        % %%%% 3. Use the T variables, do not use the posterior
        % if T_cardinal == 1 & T_oblique == 0
        %     kernel_O_use = kernel_O_cardinal(attn,1,O,:);
        % elseif T_cardinal == 0 & T_oblique == 1
        %     kernel_O_use = kernel_O_oblique(attn,1,O,:);
        % elseif T_cardinal == 0 & T_oblique == 0
        %     % uniform profile?
        %     norm_kernel_O = 1 ./ besseli(0, 0) / dimG;
        %     kernel_O_use = norm_kernel_O * exp(0 .* kernel_O_cardinal(attn,1,O,:));
        % elseif T_cardinal == 1 & T_oblique == 1
        %     % just average
        %     kernel_O_use = kernel_O_cardinal(attn,1,O,:) + kernel_O_oblique(attn,1,O,:);
        % end

        kernel_O_use = kernel_O_cardinal(attn,1,O,:) * T_cardinal + kernel_O_oblique(attn,1,O,:) * T_oblique + ...
                    (1 - T_cardinal) * (1 - T_oblique) * exp(0 .* kernel_O_cardinal(attn,1,O,:)) * 1 ./ besseli(0, 0) / dimG;
    case 'use-T-posterior'
        %%% 4. Use posterior in a winner-take-all manner only when two Ts
        %%% are the same
        if T_cardinal == 1 & T_oblique == 0
            kernel_O_use = kernel_O_cardinal(attn,1,O,:);
        elseif T_cardinal == 0 & T_oblique == 1
            kernel_O_use = kernel_O_oblique(attn,1,O,:);
        else
            if posterior_cardinal > posterior_oblique
                kernel_O_use = kernel_O_cardinal(attn,1,O,:);
            else
                kernel_O_use = kernel_O_oblique(attn,1,O,:);
            end
        end
end
end


function  kernel_O_use = get_kernel_O_use_dual(kernel_O_cardinal, kernel_O_oblique,kernel_O_uniform, T_cardinal, T_oblique)

%%%%% By Shizhao Liu 11/20/2025. Get the kernel_O to use under
%%%%% dual-task-switching mode, based on pre-defined kernel_O_cardinal,
%%%%% kernel_O_oblique, and the samples of two task variable T.

%%%%% This line of code is equivalent to:
% % if T_cardinal == 1 & T_oblique == 0
% %     % Just like cardinal
% %     kernel_O_use = kernel_O_cardinal + kernel_O_uniform;
% % elseif T_cardinal == 0 & T_oblique == 1
% %     % Just like oblique
% %     kernel_O_use = kernel_O_oblique + kernel_O_uniform;
% % elseif T_cardinal == 0 & T_oblique == 0
% %     % Uniform profile
% %     kernel_O = kernel_O_uniform + kernel_O_uniform;
% %     
% % elseif T_cardinal == 1 & T_oblique == 1
% %     % Sum of two kernels
% %     kernel_O_use = kernel_O_cardinal + kernel_O_oblique;
% % end
kernel_O_use = kernel_O_cardinal * T_cardinal + kernel_O_oblique * T_oblique + ...
            (1 - T_cardinal) * kernel_O_uniform +  (1 - T_oblique) * kernel_O_uniform;
end