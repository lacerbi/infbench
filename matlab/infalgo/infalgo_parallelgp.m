function [history,post,algoptions] = infalgo_parallelgp(algo,algoset,probstruct)

% Add algorithm to MATLAB path
BaseFolder = fileparts(mfilename('fullpath'));
AlgoFolder = 'parallel-GP-SL';
addpath(genpath([BaseFolder filesep() AlgoFolder]));

%% Base settings
algoptions.Algorithm = 'parallelgp';
algoptions.MaxFunEvals = probstruct.MaxFunEvals;
algoptions.MaxIter = Inf;
algoptions.Ninit = 10;      % # of initial locations (generated uniformly)
algoptions.BatchSize = 1;   % Batch size
algoptions.Niter = 1 + ceil((algoptions.MaxFunEvals - algoptions.Ninit) / algoptions.BatchSize);    % # of iterations of the algorithm
algoptions.Plot = 0;        % Interactive plotting: plot only final result=1 / plot after each batch=2 / plot all and pause=3 / compute estimates but no plotting=-1 (for cluster)
algoptions.AcqMethod = 'IMIQR';    % Acquisition method

if probstruct.Debug
    algoptions.TrueMean = probstruct.Post.Mean;
    algoptions.TrueCov = probstruct.Post.Cov;
end

% Options from current problem
switch algoset
    case {0,'debug'}; algoset = 'debug'; algoptions.Plot = 1;
    case {1,'base'}; algoset = 'base';           % Use defaults
    case {2,'maxiqr'}; algoset = 'maxiqr'; algoptions.AcqMethod = 'MAXIQR';
        
    otherwise
        error(['Unknown algorithm setting ''' algoset ''' for algorithm ''' algo '''.']);
end


%% GP settings
gp_opt.noise_model = 1; % 0=constant GP noise term, 1=uses bootstrapped noise variance estimates in GP, 2=known constant GP
gp_opt.meanf = 1; % 0 is zero mean GP prior, 1 enables const/lin/quadratic terms
gp_opt.hyp_upd_freq = 1; % how often GP hypers are updated using MAP-estimation
gp_opt.display_type = 'off';

%% acq settings
acq_opt.method = algoptions.AcqMethod;  % 'IMIQR', 'MAXIQR', 'EIV', 'MAXV', 'unif'
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
acq_opt.exp.nr_grid.dim2 = 50;
acq_opt.display_type = 'off';

%% MCMC settings (for sampling from GP-based posterior when dim > 2, grid-based computations always used when dim <= 2)
mcmc_opt.nchains = 5; % how many chains
mcmc_opt.nsimu = 10000; % how many samples for each chain
mcmc_opt.nfinal = 5000; % final amount of samples after concatenating the chains and thinning
mcmc_opt.display_type = 'on';

[grid_th,sim_model] = get_model(probstruct,algoptions);

%% SL estimator settings
lik_opt.method = 'sl'; % which method to obtain log-lik estimates: 'sl', 'lfire' ('lfire' not implemented')
if isfield(sim_model,'loglik_eval')
    lik_opt.method = 'exact'; % if exact log-likelihood evaluations are obtained (synthetic 2D examples only)
    gp_opt.noise_model = 2;
    gp_opt.sigma_n_const = probstruct.NoiseSigma;
end
lik_opt.sl.estimator = 'sl'; % which SL estimator: 'sl', 'ubsl', 'ublogsl'
lik_opt.sl.N = sim_model.N; % number of repeated samples computing SL at each evaluation location


%% other misc settings
other_opt.res_ind = 1:10:algoptions.Niter; % iterations when to compute TV/KL if cluster computation
other_opt.viz_ind = 1:algoptions.Niter; % iterations when to plot the figures
other_opt.viz_save = 0; % whether the save the plotted figures
other_opt.viz_demomode2d = 0; % 'demomode' for plotting 2d illustrations
other_opt.display_type = 'iter';
other_opt.output_folder = ['../results/',sim_model.name];
other_opt.output_filename = [algoptions.AcqMethod,acq_opt.batch_method];

x0 = probstruct.InitPoint;
D = size(x0,2);

grid_th.theta0 = x0;

% Do not add log prior to fcn evaluation, already passed to the algorithm
probstruct.AddLogPrior = false;

%% run the main algorithm
algo_timer = tic;
results = run_algorithm(algoptions.Ninit, algoptions.Niter, ...
    algoptions.BatchSize, algoptions.Plot, grid_th, sim_model, ...
    gp_opt, acq_opt, mcmc_opt, lik_opt, other_opt, 1);
TotalTime = toc(algo_timer);

X = results.theta_tr;
y = results.log_lik_tr;
if isfield(gp_opt,'sigma_n_const')
    s2 = gp_opt.sigma_n_const^2 * ones(size(y));
else
    s2 = results.sigma_tr.^2;    
end
s2 = max(s2,1e-12);

Niter = numel(probstruct.SaveTicks);
[history,post] = StoreAlgoResults(...
    probstruct,[],Niter,X,y,[],[],[],[],TotalTime,[],s2);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [grid_th,sim_model] = get_model(probstruct,algoptions)
%GET_MODEL Get simulation model information for parallel-GP-SL

PLB = probstruct.PLB;
PUB = probstruct.PUB;
LB = probstruct.LB;
UB = probstruct.UB;
x0 = probstruct.InitPoint;
D = size(x0,2);

bounds_range = PUB - PLB;
LB(~isfinite(LB)) = PLB(~isfinite(LB)) - 3*bounds_range(~isfinite(LB));
UB(~isfinite(UB)) = PUB(~isfinite(UB)) + 3*bounds_range(~isfinite(UB));

grid_th = theta_grid([LB',UB'],100);

sim_model.N = NaN; % for SL
sim_model.name = probstruct.Prob;
for i = 1:D; sim_model.theta_names{i} = ['x_' num2str(i)]; end

sim_model.true_theta = probstruct.Post.Mode;
sim_model.dim = D;

sim_model.loglik_eval = @(theta) infbench_func(theta,probstruct);
sim_model.prior_eval = @(theta) exp(infbench_lnprior(theta,probstruct));
if D == 2
    Nx = size(grid_th.theta2d,2);
    sim_model.true_post_pdf2d = zeros(Nx,1);
    probstruct.NoiseSigma = 0;  % Remove noise for ground truth
%     for i = 1:Nx
%         th = grid_th.theta2d(:,i)';
%         sim_model.true_post_pdf2d(i) = exp(infbench_func(th,probstruct) + infbench_lnprior(th,probstruct));
%     end
end

end
