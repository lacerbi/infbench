function [] = demo_extra_viz2d()
% Some extra visualisations to see how the acq surfaces, acquired points etc. look like.
% Go through the 2 toy problems, different levels of noise, different acq functions with
% different batch strategies.
% Running this takes ~3h

close all;

gp_abc_seed = 42;

% different acq techniques
if 1
    %toy_models = {'banana2d'};
    toy_models = {'bimodal2d'};
    %toy_models = {'banana2d','bimodal2d'};
    noise_stdevs = [1,2,5];
    
    acqs = {'unif','MAXV','MAXV','EIV','EIV','EIV','MAXIQR','MAXIQR','IMIQR','IMIQR','IMIQR'};
    batch_methods = {'seq','seq','greedy','seq','greedy','joint','seq','greedy','seq','greedy','joint'};
    b = 4;
    m = 81;
    batch_sizes = [1, 1,b,   1,b,b,     1,b,   1,b,b    ];
    ib = (m-1)/b+1;
    n_iters =     [m, m,ib, m,ib,ib, m,ib, m,ib,ib];
    ba = ceil([2,ib/2,ib]); nba = ceil([2,m/2,m]);
    viz_ind =     {nba, nba,ba, nba,ba,ba, nba,ba, nba,ba,ba};
    
else
    % for testing this script...
    toy_models = {'banana2d'};
    noise_stdevs = [1,5];

    %acqs = {'unif','MIQR'};
    acqs = {'unif','MAXIQR'};
    batch_methods = {'greedy','greedy'};
    b = 4;
    m = 81;
    batch_sizes = [1, b ];
    ib = (m-1)/b+1;
    n_iters =     [m, ib];
    ba = ceil([2,ib/2,ib]); nba = ceil([2,m/2,m]);
    viz_ind =     {nba, ba};
end

%% main settings
nr_init = 10;
%%%nr_iter = 50;
%%%batch_size = 1;
graphics_on = 2; % plot only final result=1 / plot after each batch=2 / plot all and pause=3 / compute estimates but no plotting=-1 (for cluster)


%% gp settings
%%%gp_opt.noise_model = 1; % if 1, then use bootstrapped noise variance estimates in GP
gp_opt.meanf = 1; % 0 is zero mean GP, 1 enables const/lin/quadratic terms
gp_opt.hyp_upd_freq = 1;
gp_opt.display_type = 'off';


%% acq settings
%%%acq_opt.method = 'unif';
%%%acq_opt.batch_method = 'greedy';
acq_opt.optim_alg = 'fmincon';
acq_opt.rs.nr_init = 1000; % number of evals in random search
acq_opt.fmincon.nr_inits = 1000; % number of initial points for the multistart optimization
acq_opt.fmincon.nr_inits_grad = 10; % number of the best initial points that are actually used for multistart optimization
acq_opt.fmincon.tolf = 1e-5; % error tolerances for fmincon optimizer
acq_opt.fmincon.tolx = 1e-5;
acq_opt.direct.maxevals = 1000; % max. number of function evals  (default is 20)
acq_opt.direct.maxits = 100; % max. number of iterations  (default is 10)
acq_opt.direct.maxdeep = 100; % max. number of rect. divisions (default is 100)
acq_opt.exp.is_samples = 200; % how many samples from the importance distribution 
acq_opt.exp.nr_grid.dim1 = 100; % number of grid points for grid integration
acq_opt.exp.nr_grid.dim2 = 50;
acq_opt.display_type = 'off';


%% 2d so mcmc settings not needed here
mcmc_opt = [];


%% abc/sl lik estimator settings
lik_opt.method = 'exact'; % HERE WE HAVE EXACT EVALS
gp_opt.noise_model = 0;
lik_opt.sl.estimator = 'sl'; % 'sl', 'ubsl', 'ublogsl'
lik_opt.sl.N = [];


%% other misc settings
other_opt.res_ind = []; % iterations when to compute TV/KL if cluster computation
%%%other_opt.viz_ind = 1:nr_iter; % iterations when to plot the figures
other_opt.viz_save = 1; % whether the save the plotted figures
other_opt.viz_demomode2d = 1; % for plotting 2d illustrations
other_opt.display_type = 'off';
%%%other_opt.output_folder = ['../results/',sim_model.name];
%%%other_opt.output_filename = [acq_opt.method,acq_opt.batch_method];


%% run the computations + visualisation
nacqs = length(acqs); ntoy = length(toy_models); nnoise = length(noise_stdevs);
tot = nacqs*ntoy*nnoise;
tot
cur = 0;
for i = 1:ntoy
    for j = 1:nnoise
        % get toy model settings
        [grid_th,sim_model] = get_test_model(toy_models{i},noise_stdevs(j));
        
        t0=tic;
        for k = 1:nacqs
            cur = cur + 1
            
            % update settings with current values
            nr_iter = n_iters(k);
            batch_size = batch_sizes(k);
            
            acq_opt.method = acqs{k};
            acq_opt.batch_method = batch_methods{k};
            other_opt.viz_ind = viz_ind{k};
            
            other_opt.output_folder = ['../results/viz/',toy_models{i}];
            other_opt.output_filename = [num2str(noise_stdevs(j)),'_',acq_opt.method,'_',acq_opt.batch_method];
            
            % run GP-ABC in current test scenario
            close all;
            rng(gp_abc_seed);
            run_algorithm(nr_init, nr_iter, batch_size, graphics_on, grid_th, sim_model, ...
                gp_opt, acq_opt, mcmc_opt, lik_opt, other_opt);
        end
        toc(t0)
    end
end
end




