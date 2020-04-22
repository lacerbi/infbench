function [] = sl_inference_baseline(cluster_run, id, plot_res, m)
% Computes the baseline posterior density for the real test problems whose exact 
% likelihood is unavailable. The results are saved to file to be used for comparisons in
% the GP-ABC code. This same code can also be used for plotting the resulting posterior
% once the MCMC has first been run and results saved. 
%
% Note 1: prior bounds are not set here because the MCMC never really seem to run
% issues with the boundaries anyway. This may not be the case with some other models not
% considered here though.
% Note 2: index m tells which model results to plot when 'plot_res' is 1

if nargin < 3
    plot_res = 0; m = [];
end

if cluster_run
	warning off; 
    set_paths_for_triton();
    plot_res = 0; m = [];
else
    disp('Local run...');
end

TEST_MODE = 0;
SAVE_FIGS = 1;

root = '../results/sl_baseline/';
if TEST_MODE
    root = '../results/sl_baseline_test/'; % For testing...
end

nruns = 5; % 5 runs per model
models = {'lorenz','ricker','ricker_12','gk_model','ma2'};
Ns = [100,100,100,100,100]; % N for SL
%Ns = [50,80,80,20,20]; % N for SL, old values
if TEST_MODE
    % For testing...
    ind = 5;
    models = models(ind);
    Ns = Ns(ind);
    %nruns = 1;
end
nmodels = length(models);
use_grid = 0;

if ~plot_res
    % map given id to an index to the model and repeated run
    % [run_id,model_id] = run_inds(nruns,nmodels,id);
    run_id = 1;
    model_id = 2        
    if ~cluster_run
        % some debug printing:
        model_id, run_id, cur_model=models{model_id}, cur_N=Ns(model_id)
    end
    
    % get test model settings
    [grid_th,sim_model] = get_test_model(models{model_id},[],Ns(model_id));
    
    %% SL settings
    % mcmc related settings
    sl_mcmc_opt.init = sim_model.true_theta; % USE THE 'TRUE VALUE' AS INITIAL POINT FOR MCMC
    sl_mcmc_opt.nsimu = 100000;
    sl_mcmc_opt.nchains = 1;
    sl_mcmc_opt.nfinal = 10000;
    sl_mcmc_opt.display_type = 'on';
    %...
    
    sl_opt.N = Ns(model_id);
    sl_opt.estimator = 'sl'; % 'sl','ubsl','ublogsl'
    
    
    %% run SL-MCMC inference
    sl_seed = 42 + id; % different seed for each run
    rng(sl_seed);
    if ~cluster_run
        t0=tic;
    end
    try
        [sl_samples,init_sl_samples,mcmc_diag] = standard_sl_inference(sim_model,...
            grid_th,sl_opt,sl_mcmc_opt,use_grid);
    catch err
        % save which settings etc. caused the error
        fid = fopen('errorfile.txt','a+');
        fprintf(fid, '--------------------\n');
        fprintf(fid, 'run_id: %i, model_id: %i', run_id, model_id);
        fprintf(fid, '\n');
        fprintf(fid, '%s', err.getReport('extended','hyperlinks','off'));
        fprintf(fid, '\n--------------------\n');
        fprintf(fid, '\n');
        fclose(fid);
    end
    if ~cluster_run
        toc(t0)
    end
    
    % compute marginals from the SL-MCMC samples
    post_kde = kde_for_abc(grid_th,sl_samples);
    if sim_model.dim == 2
        post_kde = vec_to_grid_matrix(post_kde, grid_th);
    end
    
    %% save thinned set of samples, initial 1000 samples and the final estimated density
    fn = [root,'/',models{model_id},'_run',num2str(run_id)];
    save(fn,'sl_samples','init_sl_samples','post_kde','mcmc_diag','grid_th',...
        'sim_model','sl_mcmc_opt','sl_opt');
end


%% Plot results
if plot_res && ~cluster_run
    % load already computed results from files...
    post = cell(nruns,1);
    samples = cell(nruns,1);
    for i = 1:nruns
        fn = [root,'/',models{m},'_run',num2str(i)];
        load(fn,'sl_samples','post_kde','grid_th','sim_model','mcmc_diag'); 
        post{i} = post_kde;
        samples{i} = sl_samples;
        mcmc_diag.R
        mcmc_diag.is_converged
    end
    d = sim_model.dim;
    
    % plot SL-posterior
    close all;
    figure(1);
    cols = {'b','r','k','g','m','--b','--r','--k','--g','--m'};
    if d ~= 2
        % 1) plot all 1d marginals to one figure:
        set(gcf,'Position',[50 1000 400*d 450]);
        for i = 1:d
            subplot(1,d,i); hold on;
            for j = 1:nruns
                plot(grid_th.theta(i,:),post{j}(:,i),[cols{j}]);
            end
            hold off;
            if isfield(sim_model,'theta_names')
                xlabel(sim_model.theta_names{i});
            end
            box on;
        end
        drawnow;
        if SAVE_FIGS
            fn = [root,'/sl_marg_post_',models{m}];
            my_export_fig(fn,'-png');
        end
        
        % 2) plot also all 2d marginals:
        if d > 2
            for j = 1:nruns
                figure(1+j);
                set(gcf,'Position',[50 50 1000 1000]);
                bivariate_marginals_plot(grid_th,sim_model,samples{j},post{j});
                drawnow;
                if SAVE_FIGS && j == 1
                    % only first figure saved to file
                    fn = [root,'/sl_post_',models{m}];
                    my_export_fig(fn,'-png');
                end
            end
        end
        
    else % d == 2
        % plot joint 2d posterior:
        nr_contour_lines = 25;
        thx = grid_th.theta(1,:);
        thy = grid_th.theta(2,:);
        set(gcf,'Position',[50 1000 400*nruns 350]);
        for j = 1:nruns
            subplot(1,nruns,j);
            contour(thx, thy, post{j}, nr_contour_lines);
            if isfield(sim_model,'theta_names')
                xlabel(sim_model.theta_names{1}); ylabel(sim_model.theta_names{2});
            end
        end
        drawnow;
        if SAVE_FIGS
            fn = [root,'/sl_post_',models{m}];
            my_export_fig(fn,'-png');
        end
    end
    
    % compute all pairwise TV/KL/L2 distances to get a (rough) estimate of variability
    if 1
        % KL is not symmetric so in KL case we could in fact compute all combinations
        % but this is not done here.
        ds = NaN(nruns,nruns);
        for i = 1:nruns
            for j = (i+1):nruns
                [kl,tv,l2] = compute_dist(grid_th, post{i}', post{j}');
                ds(i,j) = mean(tv);
            end
        end
        
        % print results
        ds
        meand = nanmean(ds(:))
        mediand = nanmedian(ds(:))
    end
end

if cluster_run
    exit;
end
end




