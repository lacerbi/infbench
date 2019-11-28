function [] = demo_loglik_gp_modelling()
% Plots *Figure 1* in the paper. 
%
% Plots a 1x2 subplot figure. Left side shows the gp model for the loglik and right side
% how its uncertainty propagates to the likelihood. 

subplot = @(m,n,p) subtightplot (m, n, p, [0.07 0.09], [0.2 0.1], [0.1 0.1]); % better subplots

close all;
rng(1234);

%% set up model etc.
load_from_file = 1;
root = '../results/demo_loglik/';
N = 100;
[grid_th,sim_model] = get_test_model('ricker_1',[],N); 
t = 30;

gp_opt.noise_model = 1; % if 1, then use bootstrapped noise variance estimates in GP
gp_opt.meanf = 1; % 0 is zero mean GP, 1 enables const/lin/quadratic terms
gp_opt.display_type = 'off';

lik_opt.method = 'sl'; 
lik_opt.sl.estimator = 'sl'; % 'sl', 'ubsl', 'ublogsl'
lik_opt.sl.N = N;

%% generate some logSL values (here simply from the uniform distribution)
th_tr = acquire_init_batch(t,grid_th);
fn_simul = [root,'loglik_cached_simulations'];
if load_from_file
    % load cached simulations
    load(fn_simul);
else
    % simulate now and save the results
    loglik_tr = NaN(t,1);
    bootvar = NaN(t,1);
    for j = 1:t
        [loglik_tr(j),bootvar(j)] = noisy_loglik_estim(sim_model,lik_opt.sl,...
            th_tr(j,:),lik_opt.method,gp_opt.noise_model);
    end
    save(fn_simul,'loglik_tr','bootvar');
end
sigma_tr = sqrt(bootvar);

%% fit GP and get the loglik estimates from it etc. 
gp = []; gp_optim_opt = [];
[gp,gp_optim_opt] = fit_gp_model(gp, gp_optim_opt, gp_opt,grid_th, loglik_tr,th_tr,sigma_tr);
P = precompute_gp_pred(gp, th_tr, loglik_tr, gp_opt);

%% compute lik variance and IQR
estim_post = post_from_gp_surrogate(grid_th,sim_model,gp,gp_opt,loglik_tr,th_tr,P,[]);

%% plot everything
figure(1);
set(gcf,'Position',[40 40 900 350]);

% 1/2: gp on loglik
subplot(1,2,1);
thx = grid_th.theta;
my_shadedplot(thx, estim_post.loglik_lb, estim_post.loglik_ub, 0.8*[1,1,1], [1,1,1]); % loglik uncertainty
hold on;
plot(thx,estim_post.eloglik,'-r'); 
plot(th_tr,loglik_tr,'k.','MarkerSize',10); 
if 1
    % plot 'errorbars' of the loglik observations
    for i = 1:t
        plot(th_tr(i)*[1 1],loglik_tr(i)+[1 -1]*1.96*sigma_tr(i),'-k');
    end
end
hold off;
box on;
xlim([thx(1),thx(end)]);
ylim([-600,30]);
xlabel('\theta');
ylabel('log-synthetic likelihood');
title('(a)');


% 2/2: uncertainty of lik
subplot(1,2,2);
my_shadedplot(thx, estim_post.post_lb, estim_post.post_ub, 0.8*[1,1,1], [1,1,1]); % post uncertainty
hold on; 
plot(thx,estim_post.epost,'-r'); % median 
if 1
    % plot variance as an additional curve
    plot(thx,sqrt(estim_post.varpost),'b--');
    % plot IQR as an additional curve?
    % ...
end
if 0
    % plot true value
    plot(sim_model.true_theta,0,'dr','MarkerSize',6,'MarkerFaceColor','r');
end
hold off;
box on;
xlim([thx(1),thx(end)]);
ylim([0,3e-6]);
xlabel('\theta');
ylabel('synthetic likelihood (rescaled)');
title('(b)');

% save to file
if 1
    fn = [root,'fig_loglik_gp'];
    my_export_fig(fn,'-transparent','-pdf');
end
end


