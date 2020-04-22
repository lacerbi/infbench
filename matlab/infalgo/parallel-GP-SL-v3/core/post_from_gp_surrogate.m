function [post,mcmc_diag] = post_from_gp_surrogate(th_grid,sim_model,gp,gp_opt,log_lik_tr,th_tr,P,mcmc_opt)
% Compute the estimate of (unnormalised) posterior from the fitted GP surrogate model of
% the log posterior. If dim <= 2, we compute it in a grid, otherwise we sample from it.

mcmc_diag = [];
d = th_grid.dim;
if ~isfield(mcmc_opt,'always_mcmc'); mcmc_opt.always_mcmc = 0; end

if d <= 2 && ~mcmc_opt.always_mcmc
    if d == 1
        th_eval_grid = th_grid.theta(:);
    else
        th_eval_grid = th_grid.theta2d';
    end
    [eft,varft] = gp_pred_fast(gp,th_tr,log_lik_tr,th_eval_grid,P);
    
    % log lik etc.
    post.eloglik = eft; % log lik
    post.varloglik = varft; % var of log lik
    post.varnloglik = varft + gp_noise_model_var(gp);
    
    % log lik CI (computed from  Gaussian densities)
    post.loglik_lb = eft - 1.96*sqrt(varft); % loglik (latent) uncertainty
    post.loglik_ub = eft + 1.96*sqrt(varft);
    post.nloglik_lb = eft - 1.96*sqrt(post.varnloglik); % new loglik measurement uncertainty
    post.nloglik_ub = eft + 1.96*sqrt(post.varnloglik);
    
    % posterior
    % NOTE: WE RESCALE THESE TO AVOID NUMERICAL OVER/UNDERFLOW!
    pr_val = sim_model.prior_eval(th_eval_grid);
    log_pr_val = log(pr_val);
    log_meanpost = log_pr_val + (eft + 0.5*varft);
    log_medpost = log_pr_val + (eft);
    % this is the same formula as the EV acq:
    log_varpost = real(2*log_pr_val + 2*(eft + varft) + (log1p(-exp(-varft)))); 
    
    %[epost,c] = exponentiate_log_array(log_meanpost); % this corresponds scaling by exp(-c)
    [epost,c] = exponentiate_log_array(log_medpost); % USE MEDIAN AS A POINT ESTIMATE INSTEAD OF MEAN
    post.epost = epost;
    post.meanpost = exp(log_meanpost - c); % could still overflow
    post.varpost = exp(log_varpost - 2*c); 
    
    % posterior CI (computed from lognormal density)
    post.post_lb = pr_val.*exp(post.loglik_lb - c);
    post.post_ub = pr_val.*exp(post.loglik_ub - c);
    post.samples = [];
    
else % dim > 2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% compute 'slices' at the true parameter
    if isfield(sim_model,'true_theta')
        gl = size(th_grid.theta,2);
        post.slice.eloglik = NaN(d,gl);
        post.slice.varloglik = NaN(d,gl);
        post.slice.epost = NaN(d,gl);
        for i = 1:d
            th_gri = th_grid.theta;
            ii = setdiff(1:d,i);
            th_gri(ii,:) = repmat(sim_model.true_theta(ii),gl,1)';
            [efti,varfti] = gp_pred_fast(gp,th_tr,log_lik_tr,th_gri',P);
            
            post.slice.eloglik(i,:) = efti;
            post.slice.varloglik(i,:) = varfti; 
            
            pr_vali = sim_model.prior_eval(th_gri');
            log_pr_vali = log(pr_vali);
            %log_meanposti = log_pr_vali + (efti + 0.5*varfti);
            log_eposti = log_pr_vali + efti; % USE MEDIAN AS A POINT ESTIMATE INSTEAD OF MEAN
            post.slice.epost(i,:) = exponentiate_log_array(log_eposti);
        end
    end
    
    %% sampling from the point estimate of posterior
    % If dim > 2, we sample from the estimated posterior density using MCMC and compute 
    % the marginals from the resulting samples using KDE.
    % Here we use the adaptive MCMC (DRAM) by Haario et al. 2006 but many other MCMC
    % algorithms could be alternatively used.
    
    display_type = mcmc_opt.display_type;
    
    % Set up variables and options for MCMC (some of these could be varied according to each particular problem)
    npar = sim_model.dim;
    
    model.ssfun = @(th,data) -2*log_post_point_estim(th,th_tr,log_lik_tr,gp,P,sim_model);
    data = [];
    model.N = 1;
    
    % Take the maximum posterior value evaluated at training data points as the initial
    % point
    f_all = model.ssfun(th_tr,[]);
    [f_opt,opt_ind] = min(f_all); % min because -2*log(.)
    init = th_tr(opt_ind,:);
    
    if ~strcmp(display_type,'off')
        disp('Initial point for MCMC:');
        init
    end
    
    params = cell(npar,1);
    for i = 1:npar
        params{i} = {sprintf('\\theta_{%d}',i),init(i),th_grid.theta(i,1),th_grid.theta(i,end)};
    end
    
    % Additional MCMC settings
    options.nsimu = mcmc_opt.nsimu;
    options.qcov = 1/10^2*diag((th_grid.theta(:,1)-th_grid.theta(:,end)).^2);
    options.method = 'am';
    options.updatesigma = 0;
    options.verbosity = ~strcmp(display_type,'off'); % no printing from mcmc
    options.waitbar = 0;
    %options.burnintime = 1000; % what is this number actually here?
    
    % Initialize results
    samples_all = NaN(mcmc_opt.nsimu,npar,mcmc_opt.nchains);
    results_all = cell(mcmc_opt.nchains,1);
    
    % Run MCMC chains!
    for i = 1:mcmc_opt.nchains
        if ~strcmp(display_type,'off')
            if i == 1
                disp('Running MCMC...');
            end
            if mcmc_opt.nchains > 1
                disp(['Chain ', num2str(i), '/', num2str(mcmc_opt.nchains)]);
            end
        end
        [results,samples] = mcmcrun(model,data,params,options);
        results_all{i} = results;
        samples_all(:,:,i) = samples;
        if i == mcmc_opt.nchains && ~strcmp(display_type,'off')
            disp('Done.');
        end
    end
    
    % Leave out burn-in (e.g. the first half of each chain)
    cl = size(samples_all,1);
    samples_all = samples_all(ceil(cl/2:cl),:,:);
    
    % Assess the convergence
    % psrf is taken from GPstuff/diag
    [R,neff,Vh,W,B,tau,thin] = psrf(samples_all);
    mcmc_diag.R = R;
    mcmc_diag.neff = neff;
    mcmc_diag.is_converged = (max(abs(R-1)) < 0.1); % one number summary of convergence assessment
    if ~strcmp(display_type,'off')
        disp(['nr chains = ', num2str(mcmc_opt.nchains)]);
        disp(['R = ',num2str(R)]);
        disp(['neff = ', num2str(neff)]);
    end
    if mcmc_diag.is_converged ~= 1
        warning('Convergence not reached when sampling from the model-based posterior estimate.');
    end
    
    % Add the chains together
    cl = size(samples_all,1);
    samples = NaN(mcmc_opt.nchains*cl,npar);
    for i = 1:mcmc_opt.nchains
        samples(1 + (i-1)*cl:i*cl,:) = samples_all(:,:,i);
    end
    
    % Thin to final length (to reduce the save size of the samples-matrix)
    final_cl = min(mcmc_opt.nfinal,size(samples,1));
    samples = samples(floor(linspace(1,size(samples,1),final_cl)),:);
    
    % Print and plot results to visually examine whether convergence is reached
    if ~strcmp(display_type,'off') 
        %results{1:end}, mean_samples = mean(samples)
        figure(50);
        clf;
        mcmcplot(samples,1:npar,results_all{1}.names,'chainpanel');
        suptitle('model-based posterior');
        set(gcf,'Position',[60 600 600 400]);
    end
    
    % Compute marginals in the grid using KDE
    post.epost = kde_for_abc(th_grid,samples,1);
    post.epost = post.epost';
    post.samples = samples;
end
end


function lp = log_post_point_estim(th_prop,th_tr,log_lik_tr,gp,P,sim_model)
% Wrapper for computing the log posterior value from the GP surrogate for the MCMC code.

[eft,varft] = gp_pred_fast(gp, th_tr, log_lik_tr, th_prop, P);
%ll = eft + 0.5*varft; % marginal mean type of estimator
ll = eft; % median i.e. plug-in type estimator, USE MEDIAN AS A POINT ESTIMATE INSTEAD OF MEAN
lp = log(sim_model.prior_eval(th_prop)) + ll;
end





