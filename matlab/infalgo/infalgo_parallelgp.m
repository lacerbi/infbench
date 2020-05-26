function [history,post,algoptions] = infalgo_parallelgp(algo,algoset,probstruct)

%% Base settings
algoptions.Algorithm = 'parallelgp';
algoptions.MaxFunEvals = probstruct.MaxFunEvals;
algoptions.MaxIter = Inf;
algoptions.Ninit = 10;      % # of initial locations (generated uniformly)
algoptions.BatchSize = 1;   % Batch size
algoptions.Niter = 1 + ceil((algoptions.MaxFunEvals - algoptions.Ninit) / algoptions.BatchSize);    % # of iterations of the algorithm
algoptions.Plot = -1;       % Interactive plotting: plot only final result=1 / plot after each batch=2 / plot all and pause=3 / compute estimates but no plotting=-1 (for cluster)
algoptions.AcqMethod = 'IMIQR';    % Acquisition method
algoptions.AcqOptInit = 2000;  % # of initial points for the multistart optimization
algoptions.AcqOptNstarts = 10; % # of the best initial points that are actually used for multistart optimization
algoptions.AcqOptMCMCnchains = 5; % how many MCMC chains for importance sampling
algoptions.AcqOptMCMCnsamples = 1e4; % how many MCMC samples for each importance sampling chain
algoptions.MCMCnchains = 5;     % how many MCMC chains for importance sampling
algoptions.MCMCnsamples = 1e4;  % how many MCMC samples for each chain
algoptions.NsamplesIS = 500;    % how many samples from the importance distribution
algoptions.GridSizeD2 = 2500;   % # grid points for 2D integral
algoptions.GPUpdateFreq = 1;    % How often GP hyperparameters are updated via MAP estimation
algoptions.Version = 'v2';      % Original version tested
algoptions.MCMCmethod = 'dram'; % Default MCMC sampler
algoptions.AcqOptMCMCmethod = 'dram'; % Default MCMC sampler for importance sampling

if probstruct.Debug
    algoptions.TrueMean = probstruct.Post.Mean;
    algoptions.TrueCov = probstruct.Post.Cov;
end

% Options from current problem
switch algoset
    case {0,'debug'}; algoset = 'debug'; algoptions.Plot = 1;
    case {1,'base'}; algoset = 'base';           % Use defaults
    case {2,'maxiqr'}; algoset = 'maxiqr'; algoptions.AcqMethod = 'MAXIQR';
    case {3,'fast'}; algoset = 'fast'; algoptions.AcqOptInit = 500; algoptions.AcqOptNstarts = 3; algoptions.AcqOptMCMCnchains = 3; algoptions.AcqOptMCMCnsamples = 1e3; algoptions.NsamplesIS = 200; algoptions.GridSizeD2 = 500;
    case {4,'fast2'}; algoset = 'fast2'; algoptions.AcqOptInit = 500; algoptions.AcqOptNstarts = 3; algoptions.AcqOptMCMCnchains = 3; algoptions.AcqOptMCMCnsamples = 1e3; algoptions.NsamplesIS = 200; algoptions.GridSizeD2 = 500; algoptions.GPUpdateFreq = 5;
    case {101,'v3'}; algoset = 'v3'; algoptions.Version = 'v3';
    case {102,'v3ss'}; algoset = 'v3ss'; algoptions.Version = 'v3'; algoptions.MCMCmethod = 'slicesample'; algoptions.MCMCnchains = 3; algoptions.AcqOptMCMCmethod = 'slicesample'; algoptions.AcqOptMCMCnchains = 2;
        
    otherwise
        error(['Unknown algorithm setting ''' algoset ''' for algorithm ''' algo '''.']);
end

%% Load the desired version

% Repo: https://github.com/mjarvenpaa/parallel-GP-SL

BaseFolder = fileparts(mfilename('fullpath'));
switch algoptions.Version
    case 'v3' % Latest version as of Apr 2020
        AlgoFolder = 'parallel-GP-SL-v3';
    otherwise % v2
        AlgoFolder = 'parallel-GP-SL';
end
addpath(genpath([BaseFolder filesep() AlgoFolder]));


%% GP settings
if ~isempty(probstruct.Noise) || probstruct.IntrinsicNoisy
    gp_opt.noise_model = 1; % 0=constant GP noise term (to be inferred), 1=heteroskedastic noise (provided), 2 known constant GP noise term
else
    gp_opt.noise_model = 2;
    gp_opt.sigma_n_const = sqrt(1e-5);
end

gp_opt.meanf = 1; % 0 is zero mean GP prior, 1 enables const/lin/quadratic terms
gp_opt.hyp_upd_freq = algoptions.GPUpdateFreq; % how often GP hypers are updated using MAP-estimation
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
acq_opt.fmincon.nr_inits = algoptions.AcqOptInit; % number of initial points for the multistart optimization
acq_opt.fmincon.nr_inits_grad = algoptions.AcqOptNstarts; % number of the best initial points that are actually used for multistart optimization
acq_opt.fmincon.tolf = 1e-5; % error tolerances for fmincon optimizer
acq_opt.fmincon.tolx = 1e-5;
acq_opt.direct.maxevals = 1000; % max. number of function evals  (default is 20)
acq_opt.direct.maxits = 100; % max. number of iterations  (default is 10)
acq_opt.direct.maxdeep = 100; % max. number of rect. divisions (default is 100)
acq_opt.exp.is_samples = algoptions.NsamplesIS; % how many samples from the importance distribution 
acq_opt.exp.nr_grid.dim1 = 100; % number of grid points for grid integration
acq_opt.exp.nr_grid.dim2 = ceil(sqrt(algoptions.GridSizeD2));
acq_opt.exp.nchains = algoptions.AcqOptMCMCnchains;
acq_opt.exp.nsimu = algoptions.AcqOptMCMCnsamples;
acq_opt.exp.method = algoptions.AcqOptMCMCmethod;
acq_opt.display_type = 'off';

%% MCMC settings (for sampling from GP-based posterior when dim > 2, grid-based computations always used when dim <= 2)
mcmc_opt.nchains = algoptions.MCMCnchains; % how many chains
mcmc_opt.nsimu = algoptions.MCMCnsamples; % how many samples for each chain
mcmc_opt.nfinal = Inf; % final amount of samples after concatenating the chains and thinning
mcmc_opt.display_type = 'off';
mcmc_opt.always_mcmc = 1; % Always do MCMC (instead of gridding)
mcmc_opt.method = algoptions.MCMCmethod; % MCMC sampler

[grid_th,sim_model] = get_model(probstruct,algoptions);

%% SL estimator settings
lik_opt.method = 'noisy_exact'; % which method to obtain log-lik estimates
%if isfield(sim_model,'loglik_eval')
%    lik_opt.method = 'exact'; % if exact log-likelihood evaluations are obtained (synthetic 2D examples only)
%    gp_opt.noise_model = 2;
%    gp_opt.sigma_n_const = probstruct.NoiseSigma;
%end
lik_opt.sl.estimator = 'sl'; % which SL estimator: 'sl', 'ubsl', 'ublogsl'
lik_opt.sl.N = sim_model.N; % number of repeated samples computing SL at each evaluation location


%% other misc settings
other_opt.res_ind = 1:5:algoptions.Niter; % iterations when to compute TV/KL if cluster computation
other_opt.viz_ind = 1:algoptions.Niter; % iterations when to plot the figures
other_opt.viz_save = 0; % whether the save the plotted figures
other_opt.viz_demomode2d = 0; % 'demomode' for plotting 2d illustrations
other_opt.display_type = 'iter';
other_opt.output_folder = ['../results/',sim_model.name];
other_opt.output_filename = [algoptions.AcqMethod,acq_opt.batch_method];

PLB = probstruct.PLB;
PUB = probstruct.PUB;
x0 = probstruct.InitPoint;
D = size(x0,2);

% Initial point reservoir drawn within the plausible box
N0 = algoptions.Ninit - size(x0,1);
x0 = [x0; bsxfun(@plus,bsxfun(@times,rand(N0,D),PUB-PLB),PLB)];
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

% Store samples at each iteration
Nmax = results.funccount(end);
idx = probstruct.SaveTicks(probstruct.SaveTicks <= Nmax);
time_mcmc = zeros(size(idx));
time_algo = zeros(size(idx));
for i = 1:numel(idx)
    ii = find(idx(i) == results.funccount,1);
    if isempty(ii); continue; end
    samples_iter{i} = results.samples{ii};
    time_mcmc(i) = results.time_mcmc(ii);
    time_algo(i) = results.runtime(ii);
end

% Adjust total runtime (only consider the last sampling time)
TotalTime = TotalTime - cumsum(time_mcmc) + time_mcmc(end);

Niter = numel(probstruct.SaveTicks);
[history,post] = StoreAlgoResults(...
    probstruct,[],Niter,X,y,[],[],[],[],TotalTime,samples_iter,s2);

% Correct for sampling time (only consider the one for the current iteration)
history.ElapsedTime = time_algo - cumsum(time_mcmc) + time_mcmc;


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
