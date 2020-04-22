function [] = test_main_code()
% An example script for inference with various test models.

clear; close('all'); format compact;

use_profiler = 0; % whether to use matlab profiler to analyze the code
seed = 12345;
rng(seed);

%% main settings
nr_init = 20; % number of initial locations (generated uniformly)
nr_iter = 41; % number of iterations of the algorithm
batch_size = 1; % batch size
graphics_on = 2; % interactive plotting: plot only final result=1 / plot after each batch=2 / plot all and pause=3 / compute estimates but no plotting=-1 (for cluster)


%% select test problem
%[grid_th,sim_model] = get_test_model('gaussian1d',[],100); 
%[grid_th,sim_model] = get_test_model('gaussian2d',[],100); 
%[grid_th,sim_model] = get_test_model('gaussian_3',[],100);
%[grid_th,sim_model] = get_test_model('ricker_1',[],500); 
%[grid_th,sim_model] = get_test_model('ricker_12',[],500); 
%[grid_th,sim_model] = get_test_model('ricker',[],100);
%[grid_th,sim_model] = get_test_model('simple2d',2,[]);
[grid_th,sim_model] = get_test_model('banana2d',2,[]);
%[grid_th,sim_model] = get_test_model('bimodal2d',2,[]);
%[grid_th,sim_model,samples] = get_test_model('simple6d',2,[]);
%[grid_th,sim_model,samples] = get_test_model('banana6d',2,[]);
%[grid_th,sim_model,samples] = get_test_model('bimodal6d',2,[]);
%[grid_th,sim_model] = get_test_model('gk_model',[],100);
%[grid_th,sim_model] = get_test_model('lorenz',[],100);
%[grid_th,sim_model] = get_test_model('ma2',[],100);


%% GP settings
gp_opt.noise_model = 1; % 0=constant GP noise term, 1=uses bootstrapped noise variance estimates in GP
gp_opt.meanf = 1; % 0 is zero mean GP prior, 1 enables const/lin/quadratic terms
gp_opt.hyp_upd_freq = 1; % how often GP hypers are updated using MAP-estimation
gp_opt.display_type = 'on';


%% acq settings
%acq_opt.method = 'MAXV'; % which acquisition function (aka design criterion)
%acq_opt.method = 'EIV';
%acq_opt.method = 'MAXIQR';
acq_opt.method = 'IMIQR';
%acq_opt.method = 'unif'; % this is the same as 'RAND' in the paper
acq_opt.batch_method = 'greedy'; % which method used to construct batches for parallel computation
%acq_opt.batch_method = 'joint';
acq_opt.optim_alg = 'fmincon'; % which optimisation algoritm for obtaining the design points ('fmincon' recommended, others not carefully tested)
%acq_opt.optim_alg = 'cmaes';
%acq_opt.optim_alg = 'grid';
%acq_opt.optim_alg = 'rs'; % random search
%acq_opt.optim_alg = 'direct';
acq_opt.rs.nr_init = 1000; % number of evals in random search
acq_opt.fmincon.nr_inits = 2000; % number of initial points for the multistart optimization
acq_opt.fmincon.nr_inits_grad = 10; % number of the best initial points that are actually used for multistart optimization
acq_opt.fmincon.tolf = 1e-5; % error tolerances for fmincon optimizer
acq_opt.fmincon.tolx = 1e-5;
acq_opt.direct.maxevals = 1000; % max. number of function evals  (default is 20)
acq_opt.direct.maxits = 100; % max. number of iterations  (default is 10)
acq_opt.direct.maxdeep = 100; % max. number of rect. divisions (default is 100)
acq_opt.exp.is_samples = 500; % how many samples from the importance distribution 
acq_opt.exp.nr_grid.dim1 = 100; % number of grid points for grid integration
acq_opt.exp.nr_grid.dim2 = 70;
acq_opt.display_type = 'off';


%% MCMC settings (for sampling from GP-based posterior when dim > 2, grid-based computations always used when dim <= 2)
mcmc_opt.nchains = 5; % how many chains
mcmc_opt.nsimu = 20000; % how many samples for each chain
mcmc_opt.nfinal = 10000; % final amount of samples after concatenating the chains and thinning
mcmc_opt.display_type = 'on';


%% SL estimator settings
lik_opt.method = 'sl'; % which method to obtain log-lik estimates: 'sl', 'lfire' ('lfire' not implemented')
if isfield(sim_model,'loglik_eval')
    lik_opt.method = 'exact'; % if exact log-likelihood evaluations are obtained (synthetic 2D examples only)
    gp_opt.noise_model = 0;
end
lik_opt.sl.estimator = 'sl'; % which SL estimator: 'sl', 'ubsl', 'ublogsl'
lik_opt.sl.N = sim_model.N; % number of repeated samples computing SL at each evaluation location
lik_opt.sl.robust_bootvar = 1; % if 1 uses robust variance estimator in bootstrap 

%% other misc settings
other_opt.res_ind = 1:10:nr_iter; % iterations when to compute TV/KL if cluster computation
other_opt.viz_ind = 1:nr_iter; % iterations when to plot the figures
other_opt.viz_save = 0; % whether the save the plotted figures
other_opt.viz_demomode2d = 0; % 'demomode' for plotting 2d illustrations
other_opt.display_type = 'iter';
other_opt.output_folder = ['../results/',sim_model.name];
other_opt.output_filename = [acq_opt.method,acq_opt.batch_method];

% Note that some settings related to the GP prior and MCMC are currently hard-coded and not included here.
% However, it may be useful to adjust also these for particular inference problem at hand. 

if use_profiler
    profile on;
end

if exist('samples','var')
    sim_model.samples = samples;
end

%% run the main algorithm
results = run_algorithm(nr_init, nr_iter, batch_size, graphics_on, grid_th, sim_model, ...
    gp_opt, acq_opt, mcmc_opt, lik_opt, other_opt);

if use_profiler
    profile viewer;
end
end




