function [results] = run_algorithm(nr_init, nr_iter, batch_size, graphics_on, th_grid,...
    sim_model, gp_opt, acq_opt, mcmc_opt, lik_opt, other_opt, return_lvl)
% This is the main function for inference with 'GP-ABC' where the log-likelihood is 
% modelled with GP surrogate. 

if nargin < 12
    return_lvl = 1;
end

%% initialize
gp = []; gp_optim_opt = []; P = [];
th_tr = []; % acquired parameters
loglik_tr = []; % corresponding log lik values
sigma_tr = []; % noise variances (in case they are obtained)
results.iter = cell(nr_iter,1); % results are saved here

%% run the main algorithm
for iter = 1:nr_iter
    if strcmp(other_opt.display_type,'iter')
        if iter == 1
            disp(' '); disp('initial batch:');
        else
            disp(' '); disp(['batch: ', num2str(iter), '/', num2str(nr_iter)]);
        end
        disp('--------------');
    end
    
    %% select next evaluation location(s)
    if iter == 1
        % select initial evaluation locations
        next_th = acquire_init_batch(nr_init,th_grid);
    else
        % select evaluation locations by using an acq function (aka design criterion)
        next_th = acquire_next_batch(th_grid,batch_size,iter,gp,loglik_tr,th_tr,...
            sigma_tr,P,sim_model,gp_opt,acq_opt,other_opt,graphics_on);
    end
    
    
    %% evaluate log-post at new locations (this is the comput. costly part and should be 
    %% in reality done in parallel)
    next_lik = zeros(size(next_th,1),1);
    next_bootvar = zeros(size(next_th,1),1);
    for j = 1:size(next_th,1) % different points in the current batch
        if strcmp(lik_opt.method,'lfire')
            if iter == 1
                % in lfire, we also need to generate data from the marginal lik
                %...
            end
            error('lfire not implemented');
        else % SL or exact
            [next_lik(j),next_bootvar(j)] = noisy_loglik_estim(sim_model,lik_opt.sl,...
                next_th(j,:),lik_opt.method,gp_opt.noise_model);
        end
    end

    %% update training data
    th_tr = [th_tr; next_th]; % each row == one acquired location
    loglik_tr = [loglik_tr; next_lik]; 
    sigma_tr = [sigma_tr; sqrt(next_bootvar)];
    if strcmp(other_opt.display_type,'iter')
        % print the above values
        next_th = next_th', next_lik = next_lik(:)'
        next_bootsd = sqrt(next_bootvar(:)')
    end
    
    %% estimate GP model hyperparameters
    if iter == 1 || rem(iter,gp_opt.hyp_upd_freq) == 0
        [gp,gp_optim_opt] = fit_gp_model(gp, gp_optim_opt, gp_opt,th_grid, loglik_tr,th_tr,sigma_tr);
    end
    
    %% precompute some GP quantities
    P = precompute_gp_pred(gp, th_tr, loglik_tr, gp_opt);
    
    %% compute current posterior using the GP surrogate and compare to the true baseline 
    %% (if available) by computing TV/KL
    if graphics_on >= 2 || (graphics_on == 1 && iter == nr_iter) || (graphics_on == -1 && any(other_opt.res_ind==iter))
    	[estim_post,mcmc_diag] = post_from_gp_surrogate(th_grid,sim_model,gp,gp_opt,loglik_tr,th_tr,P,mcmc_opt); 
        
        res_i = compare_to_true_baseline(th_grid, estim_post, sim_model);
        res_i.mcmc_diag = mcmc_diag;
        if return_lvl >= 2 && other_opt.res_ind(end)==iter
            % save also the final model-based posterior estimate!
            res_i.estim_post = estim_post.epost;
        end
        results.iter{iter} = res_i;
    end

    %% compare to the true baseline (if available) visually
    if (graphics_on >= 2 && any(other_opt.viz_ind==iter)) || (graphics_on == 1 && iter == nr_iter)
        figure(1);
        clf;
        plot_comparison_fig(th_grid,loglik_tr,th_tr,estim_post,sim_model,nr_init,batch_size,res_i,other_opt);
        drawnow;
        if graphics_on >= 3
            pause;
        end
        if other_opt.viz_save
            fn = [other_opt.output_folder,'/',other_opt.output_filename,'_iter',num2str(iter),'res'];
            my_export_fig(fn,'-transparent','-pdf');
        end
    end
end

% gather acquired training data and settings
if return_lvl >= 1
    results.theta_tr = th_tr;
    results.log_lik_tr = loglik_tr;
    results.sigma_tr = sigma_tr;
    results.opt_all = gather_settings(nr_init, nr_iter, batch_size, graphics_on, th_grid,...
        sim_model, gp_opt, acq_opt, mcmc_opt, lik_opt, other_opt);
end
end


function opt_all = gather_settings(nr_init, nr_iter, batch_size, graphics_on, th_grid,...
    sim_model, gp_opt, acq_opt, mcmc_opt, lik_opt, other_opt)
% Gather all settings to a one data structure.

opt_all.nr_init = nr_init;
opt_all.nr_iter = nr_iter;
opt_all.batch_size = batch_size;
opt_all.th_grid = th_grid;
opt_all.sim_model = sim_model;
opt_all.gp_opt = gp_opt;
opt_all.acq_opt = acq_opt;
opt_all.mcmc_opt = mcmc_opt;
opt_all.lik_opt = lik_opt;
opt_all.other_opt = other_opt;
end





