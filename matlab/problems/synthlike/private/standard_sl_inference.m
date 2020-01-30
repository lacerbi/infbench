function [estim_density,init_samples,mcmc_diag] = standard_sl_inference(sim_model,...
    grid_th,opt,mcmc_opt,use_grid)
% Performs grid or (pseudo-marginal) MCMC based inference with synthetic likelihood.
%
% if dim == 1 and 'use_grid' option is on, 'estim_density' shows SL values in a grid, 
% otherwise this variable contains posterior samples generated using MCMC-SL. 
% Variable 'init_samples' contains the first 1000 samples (whenever applicable). 

if use_grid && sim_model.dim == 1
    % evaluate SL in a 1d grid - for demonstration only
    
    grid_len = length(grid_th.theta);
    sl_grid = NaN(grid_len,1);
    for i = 1:grid_len
        sl_grid(i) = noisy_loglik_estim(sim_model,opt,grid_th.theta(i)) ...
            + log(sim_model.prior_eval(theta(i)));
    end
    estim_density = sl_grid;
    init_samples = [];
    mcmc_diag = [];
    
elseif use_grid && sim_model.dim == 2
    % evaluate SL in a 2d grid
    
    error('Not implemented.');
    
else % grid not used and/or dim > 2
    
    % Here we use the adaptive MCMC (DRAM) by Haario et al. 2006 but many other MCMC
    % algorithms could be alternatively used
    
    display_type = mcmc_opt.display_type;
    
    model.ssfun = @(x,data) -2*(noisy_loglik_estim(sim_model,opt,x) + log(sim_model.prior_eval(x)));
    data = [];
    model.N = 1;
    
    npar = sim_model.dim;
    params = cell(npar,1);
    for i = 1:npar
        params{i} = {sprintf('\\theta_{%d}',i), mcmc_opt.init(i), grid_th.theta(i,1), grid_th.theta(i,end)};
    end
    
    % Additional MCMC settings
    options.nsimu = mcmc_opt.nsimu;
    options.qcov = 1/10^2*diag((grid_th.theta(:,1)-grid_th.theta(:,end)).^2);
    options.method = 'am';
    options.updatesigma = 0;
    options.verbosity = 0; % no printing from mcmc
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
        warning('Convergence possibly not reached when sampling in SL.');
    end
    
    % Add the chains together
    cl = size(samples_all,1);
    samples = NaN(mcmc_opt.nchains*cl,npar);
    for i = 1:mcmc_opt.nchains
        samples(1 + (i-1)*cl:i*cl,:) = samples_all(:,:,i);
    end
    
    % Extract the first 1000 samples 
    % (If multiple chains, only 1000 samples from the first chain are returned)
    if nargout > 1
        init_n = min(1000,size(samples,1));
        init_samples = samples(1:init_n,:);
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
        set(gcf,'Position',[60 600 600 400]);
    end
    estim_density = samples;
end
end




