function [grid_th,sim_model,samples] = get_test_model(name, noise_stdev, N, data_seed)
% Returns settings for different test scenarios with simulator models. 
%
% INPUT:
% name: which simulation model/testing scenario
% noise_stdev: sets the amount of noise in log-lik evaluations (only applicable for synthetic 2D models)
% N: number of repeated simulations at each proposed point for SL
% data_seed: seed for generating data (not applicable for all cases)
%
% OUTPUT:
% grid_th: structure for grid-based computations
% sim_model: structure for various settings etc. related to the simulation model
% samples: SL posterior samples (returned separately so that they are not redundantly
% saved to multiple places (if not available, empty matrix is returned))

if nargin < 2
    noise_stdev = [];
end
if nargin < 3
    N = [];
end

% Fix seed so that the same 'true' data is always generated from the model (in some cases
% external data is used and this will not have any effect)
if nargin < 4
    data_seed = 1234567;
end
cur_seed = rng;
rng(data_seed);

samples = [];
sim_model.N = N;
if strcmp(name,'gaussian1d')
    % Simple 1d Gaussian test case, similar example as in LFIRE paper
    
    grid_th = theta_grid([-20,20],200);
    
    % simulation model specification
    sim_model.name = 'Gaussian 1d';
    sim_model.true_theta = 2.3;
    sim_model.dim = 1;
    sim_model.fixed_var = 3; 
    sim_model.gen = @(theta,n) theta + sqrt(sim_model.fixed_var)*randn(n,1); % generate one data set
    sim_model.n_data = 1;
    sim_model.summary_dim = 1;
    sim_model.comp_summaries = @(data,obs_data) mean(data);
    
    % data from the model and true posterior
    sim_model.data = sim_model.gen(sim_model.true_theta, sim_model.n_data);
    sim_model.data_avg = mean(sim_model.data);
    sim_model.summary_true = sim_model.data_avg;
    sim_model.prior_eval = @(theta)ones(size(theta));
    sim_model.true_post_pdf = normpdf(grid_th.theta, sim_model.data_avg, sqrt(sim_model.fixed_var/sim_model.n_data));
    
elseif strcmp(name,'gaussian2d')
    % Simple 2d Gaussian model
    
    grid_th = theta_grid([-20,20;-20,20],75);
    
    % simulation model specification
    sim_model.name = 'Gaussian 2d';
    sim_model.true_theta = [0,0];
    sim_model.dim = 2;
    cc = 0.25; % correlation between elements (same for all)
    cvar = 1^2; % variance of each element
    sim_model.fixed_var = cvar*(cc*ones(2) + (1-cc)*eye(2));
    sim_model.gen = @(theta,n) mvnrnd(theta(:)',sim_model.fixed_var);
    sim_model.n_data = 1;
    sim_model.summary_dim = 2;
    sim_model.comp_summaries = @(data,obs_data) data;
    
    % data from the model and true posterior
    sim_model.data = sim_model.gen(sim_model.true_theta, sim_model.n_data);
    sim_model.data_avg = sim_model.data;
    sim_model.summary_true = sim_model.data_avg;
    sim_model.prior_eval = @(theta) ones(size(theta,1),1);
    
    sim_model.post_mean = sim_model.data_avg;
    sim_model.post_var = sim_model.fixed_var/sim_model.n_data;
    sim_model.true_post_pdf = [normpdf(grid_th.theta(1,:), sim_model.post_mean(1), sqrt(sim_model.post_var(1,1))),...
                               normpdf(grid_th.theta(2,:), sim_model.post_mean(2), sqrt(sim_model.post_var(2,2)))];
    sim_model.true_post_pdf2d = mvnpdf(grid_th.theta2d',sim_model.post_mean,sim_model.post_var);
    
elseif strcmp(name(1:min(end,9)),'gaussian_')
    % Simple Gaussian test problem, 1 data point, name must be 'gaussian_X', where
    % X=dimension
    
    d = str2double(name(10:end));
    if(d<3); error('Wrong dimension.'); end;
    grid_th = theta_grid(repmat([-20,20],d,1),100);
    
    % simulation model specification
    sim_model.name = ['Gaussian ',num2str(d),'d'];
    sim_model.true_theta = ones(1,d);
    sim_model.dim = d;
    cc = 0.25; % correlation between elements (same for all)
    cvar = 1; % variance of each element
    sim_model.fixed_var = cvar*(cc*ones(d) + (1-cc)*eye(d));
    sim_model.gen = @(theta,n) mvnrnd(theta(:)',sim_model.fixed_var); % generate one data set
    sim_model.n_data = 1;
    sim_model.summary_dim = d;
    sim_model.comp_summaries = @(data,obs_data) data;
    
    % data from the model and true posterior
    sim_model.data = sim_model.gen(sim_model.true_theta, sim_model.n_data);
    sim_model.data_avg = sim_model.data;
    sim_model.summary_true = sim_model.data_avg;
    sim_model.prior_eval = @(theta) ones(size(theta,1),1);
    
    sim_model.post_mean = sim_model.data_avg;
    sim_model.post_var = diag(sim_model.fixed_var/sim_model.n_data);
    sim_model.true_post_pdf = NaN(d,100);
    for i = 1:d
        sim_model.true_post_pdf(i,:)=normpdf(grid_th.theta(i,:), sim_model.data_avg(i), sqrt(sim_model.post_var(i)));
    end
    
elseif strcmp(name,'ricker')
    % Ricker model benchmark problem from literature, 3 parameters
    
    warning off;
    
    grid_th = theta_grid([3,5; 4,20; 0,0.8],100);
    %grid_th = theta_grid([3.5,4.5; 7,13; 0,0.4],100); % TESTING: smaller param space
    
    % simulation model specification
    sim_model.name = 'Ricker';
    sim_model.theta_names = {'log(r)','\phi','\sigma_e'};
    sim_model.true_theta = [3.8,10,0.3];
    sim_model.dim = 3;
    T = 50;
    sim_model.gen = @(theta,n) simulate_ricker(theta,1,T); % generate one data set
    sim_model.n_data = T;
    sim_model.summary_dim = 13;
    sim_model.comp_summaries = @(data,obs_data) ricker_summstats(data,obs_data);
    
    % data from the model and true posterior
    sim_model.data = sim_model.gen(sim_model.true_theta, sim_model.n_data);
    sim_model.summary_true = ricker_summstats(sim_model.data,sim_model.data);
    sim_model.prior_eval = @(theta) ones(size(theta,1),1);
    [sim_model.true_post_pdf,samples] = sl_baseline_get(name);
    sim_model.true_post_pdf = sim_model.true_post_pdf';
    
elseif strcmp(name(1:min(7,end)),'ricker_')
    % Ricker model where one or two of the parameters are fixed to their 'true' values and
    % the third one is to be estimated.
    % ricker_X, where X=which parameter is estimated.
    
    warning off;
    
    i = str2double(name(8:end)); if(i>10); i=sort([floor(i/10),rem(i,10)]); end; 
    bds = [3,5; 4,20; 0,0.8];
    %bds = [3.5,5; 6,15; 0.1,0.6]; % TESTING: smaller param space
    true_th = [3.8,10,0.3];
    %true_th = [5,10,0]; % TESTING ONLY, DIFFICULT CASE!!!
    nms = {'log(r)','\phi','\sigma_e'};
    grid_th = theta_grid(bds(i,:),100); 
    
    % simulation model specification
    sim_model.name = 'Ricker';
    sim_model.theta_names = nms(i);
    sim_model.true_theta = true_th(i);
    sim_model.dim = length(i);
    T = 50;
    sim_model.gen = @(theta,n) ricker_wrapper(theta,n,i,true_th,T); % generate one data set
    sim_model.n_data = T;
    sim_model.summary_dim = 13;
    sim_model.comp_summaries = @(data,obs_data) ricker_summstats(data,obs_data);
    
    % data from the model and true posterior
    sim_model.data = sim_model.gen(sim_model.true_theta, sim_model.n_data);
    sim_model.summary_true = ricker_summstats(sim_model.data,sim_model.data);
    sim_model.prior_eval = @(theta)ones(size(theta,1),1);
    if strcmp(name,'ricker_12')
        sim_model.true_post_pdf2d = sl_baseline_get(name);
    else
        sim_model.true_post_pdf = [];
    end
    
elseif strcmp(name,'simple2d')
    % Simple 2d toy model where an exact but noisy log lik estimate can be computed
    
    grid_th = theta_grid(1*[-16,16;-16,16],75);
    
    sim_model.N = NaN; % for SL
    sim_model.name = 'Simple';
    sim_model.theta_names = {'x','y'};
    sigma_n = noise_stdev;
    if 0
        sigmaf = @(theta)sigma_n*(abs(theta(:,2)).^(1)); % noise depends on input (in a known way!)
    else
        sigmaf = @(theta)sigma_n*ones(size(theta,1),1); % constant noise
    end

    sim_model.true_theta = [1, 1];
    sim_model.dim = 2;
    
    sim_model.loglik_eval = @(theta) log_simple_pdf(theta) + sigmaf(theta).*randn(size(theta,1),1);
    sim_model.prior_eval = @(theta) ones(size(theta,1),1);
    sim_model.true_post_pdf2d = exp(log_simple_pdf(grid_th.theta2d'));
    
elseif strcmp(name,'banana2d')
    % Simple 2d toy model where an exact but noisy log lik estimate can be computed
    
    grid_th = theta_grid(2*[-3,3;-10,1],75);
    
    sim_model.N = NaN; % for SL
    sim_model.name = 'Banana';
    sim_model.theta_names = {'x','y'};
    ab = [1,1];
    sigma_n = noise_stdev;
    if 0
        sigmaf = @(theta)sigma_n*(abs(theta(:,2)).^(1)); % noise depends on input (in a known way!)
    else
        sigmaf = @(theta)sigma_n*ones(size(theta,1),1); % constant noise
    end

    sim_model.true_theta = [0, -ab(1)^2*ab(2)];
    sim_model.dim = 2;
    
    sim_model.loglik_eval = @(theta) log_banana_pdf(theta,ab) + sigmaf(theta).*randn(size(theta,1),1);
    sim_model.prior_eval = @(theta) ones(size(theta,1),1);
    sim_model.true_post_pdf2d = exp(log_banana_pdf(grid_th.theta2d',ab));

elseif strcmp(name,'bimodal2d')
    % Simple 2d toy model where an exact but noisy log lik estimate can be computed
    
    grid_th = theta_grid(2*[-3,3;-3,3],75);
    
    sim_model.N = NaN; % for SL
    sim_model.name = 'Bimodal';
    sim_model.theta_names = {'x','y'};
    sigma_n = noise_stdev;
    if 0
        sigmaf = @(theta)sigma_n*(abs(theta(:,2)).^(1)); % noise depends on input (in a known way!)
    else
        sigmaf = @(theta)sigma_n*ones(size(theta,1),1); % constant noise
    end

    sim_model.true_theta = [0, sqrt(2)];
    sim_model.dim = 2;
    
    sim_model.loglik_eval = @(theta) log_bimodal_pdf(theta) + sigmaf(theta).*randn(size(theta,1),1);
    sim_model.prior_eval = @(theta) ones(size(theta,1),1);
    sim_model.true_post_pdf2d = exp(log_bimodal_pdf(grid_th.theta2d'));
    
    %#####################################################################################
elseif strcmp(name,'simple6d')
    % Simple 2d toy model where an exact but noisy log lik estimate can be computed
    
    d = 6;
    grid_th = theta_grid(repmat([-16,16],d,1),200);
    
    sim_model.N = NaN; % for SL
    sim_model.name = ['Simple (',num2str(d),'D)'];
    sigma_n = noise_stdev;
    sigmaf = @(theta)sigma_n*ones(size(theta,1),1); % constant noise
    
    sim_model.true_theta = ones(1,d);
    sim_model.dim = d;
    
    sim_model.loglik_2d = @log_simple_pdf;
    sim_model.loglik_eval = @(theta) log_pdf_nd(theta,@log_simple_pdf,d) + sigmaf(theta).*randn(size(theta,1),1);
    sim_model.prior_eval = @(theta) ones(size(theta,1),1);
    [sim_model.true_post_pdf,samples] = high_dim_toy_model_inference(grid_th, sim_model);
    
elseif strcmp(name,'banana6d')
    % Simple 2d toy model where an exact but noisy log lik estimate can be computed
    
    d = 6;
    grid_th = theta_grid(repmat(2*[-3,3;-10,1],d/2,1),200);
    
    sim_model.N = NaN; % for SL
    sim_model.name = ['Banana (',num2str(d),'D)'];
    ab = [1,1];
    sigma_n = noise_stdev;
    sigmaf = @(theta)sigma_n*ones(size(theta,1),1); % constant noise

    sim_model.true_theta = repmat([0, -ab(1)^2*ab(2)],1,d/2);
    sim_model.dim = d;
    
    log_pdf = @(theta) log_banana_pdf(theta,ab);
    sim_model.loglik_2d = log_pdf;
    sim_model.loglik_eval = @(theta) log_pdf_nd(theta,log_pdf,d) + sigmaf(theta).*randn(size(theta,1),1);
    sim_model.prior_eval = @(theta) ones(size(theta,1),1);
    [sim_model.true_post_pdf,samples] = high_dim_toy_model_inference(grid_th, sim_model);

elseif strcmp(name,'bimodal6d')
    % Simple 2d toy model where an exact but noisy log lik estimate can be computed
    
    d = 6;
    grid_th = theta_grid(repmat(2*[-3,3],d,1),200);
    
    sim_model.N = NaN; % for SL
    sim_model.name = ['Bimodal (',num2str(d),'D)'];
    sigma_n = noise_stdev;
    sigmaf = @(theta)sigma_n*ones(size(theta,1),1); % constant noise

    sim_model.true_theta = repmat([0, sqrt(2)],1,d/2);
    sim_model.dim = d;
    
    sim_model.loglik_2d = @log_bimodal_pdf;
    sim_model.loglik_eval = @(theta) log_pdf_nd(theta,@log_bimodal_pdf,d) + sigmaf(theta).*randn(size(theta,1),1);
    sim_model.prior_eval = @(theta) ones(size(theta,1),1);
    [sim_model.true_post_pdf,samples] = high_dim_toy_model_inference(grid_th, sim_model);
    
    %#####################################################################################
elseif strcmp(name,'cell_model')
    % Cell biology model from Bayesian synthetic likelihood paper
    % NOTE: We use the data from the paper.
    % NOTE2: The posterior is very peaked here and the log-likelihood is essentially -Inf 
    % in large region of the parameter space -> this example is not very suitable for
    % current GP-ABC algorithms where loglik is modelled with GP.
    
    oned = 0; % if second param fixed
    
    if oned
        grid_th = theta_grid([0,1],100);
        %grid_th = theta_grid([0.2,0.6],100); % TESTING: smaller param space
    else
        grid_th = theta_grid([0,1; 0,1],100);
        %grid_th = theta_grid([0.3,0.38; 0,2.5e-3],100); % TESTING: smaller param space
    end
    load('all_locations_simulated.mat'); % Pm, Pp, S, Yinit
    
    % simulation model specification
    sim_model.name = 'Cell';
    sim_model.theta_names = {'P_m','P_p'};
    if oned
        sim_model.true_theta = 0.35;
        sim_model.dim = 1;
        sim_model.gen = @(theta,n) run_cell_model(S,Yinit,[theta,0.001],n); % generate one data set
    else
        sim_model.true_theta = [0.35, 0.001];
        sim_model.dim = 2;
        sim_model.gen = @(theta,n) run_cell_model(S,Yinit,theta,n); % generate one data set
    end
    sim_model.summary_dim = 145;
    sim_model.n_data = 1;
    sim_model.comp_summaries = @(data,obs_data) cell_model_summary(Yinit,data);
    
    % data from the model and true posterior
    sim_model.data = []; % NOTE: true data is always used and loaded from file
    sim_model.summary_true = cell_model_true_summary();
    sim_model.prior_eval = @(theta) ones(size(theta,1),1);
    sim_model.true_post_pdf = [];
    
elseif strcmp(name,'gk_model')
    % g-and-k model as in Bayesian synthetic likelihood paper
    % NOTE: We use the data from the paper.
    
    oned = 0; % if only first param estimated
    
    if oned
        grid_th = theta_grid([0 6],100);
        %grid_th = theta_grid([2.5 3.5],100); % small param space
    else
        grid_th = theta_grid([2 4; 0 3; 1 4; 0 2],100);
        grid_th = theta_grid([2.5 3.5; 0.5 1.5; 1.5 2.5; 0.3 0.7],100); % medium-sized param space
        %grid_th = theta_grid([2.5 3.5; 0.8 1.2; 1.8 2.3; 0.4 0.6],100); % small param space
    end
    
    % load true data
    load('data_gandk.mat'); % y
    load('Tskew_to_GandK_MLEs.mat'); % a,b,mut,sigmas, %skew t MLEs based on observed data
    
    % simulation model specification
    sim_model.name = 'g-and-k';
    sim_model.theta_names = {'a','b','g','k'};
    sim_model.n_data = 1;
    if oned
        sim_model.true_theta = [3];
        sim_model.dim = 1;
        sim_model.gen = @(theta,n) simulate_gandk(length(y),[theta,1,2,0.5]); % generate one data set
        sim_model.summary_dim = 4;
    else
        sim_model.true_theta = [3,1,2,0.5]; % was initially [2.99,1,2.08,0.5] which was rough estimate based on the figure
        sim_model.dim = 4;
        sim_model.gen = @(theta,n) simulate_gandk(length(y),theta); % generate one data set
        sim_model.summary_dim = 4;
    end
    sim_model.comp_summaries = @(data,obs_data) Scores(a,b,mut,sigmas,data);
    
    % data from the model and true posterior
    sim_model.data = []; % NOTE: true data is always used and loaded from file
    sim_model.summary_true = [0,0,0,0];
    sim_model.prior_eval = @(theta) ones(size(theta,1),1);
    [sim_model.true_post_pdf,samples] = sl_baseline_get(name);
    sim_model.true_post_pdf = sim_model.true_post_pdf';
    
elseif strcmp(name,'lorenz')
    % Lorenz model benchmark problem, 2 parameters, matlab code from LFIRE paper
    % NOTE: We generate data from the model and use it to infer the parameters.
    
    %grid_th = theta_grid([0.5,3.5; 0,0.3],75); % same as LFIRE paper
    %grid_th = theta_grid([1.3,3; 0,0.2],75); % TESTING: smaller param space
    grid_th = theta_grid([0,5; 0,0.5],75); % TESTING: larger param space
    
    % simulation model specification
    sim_model.name = 'Lorenz';
    sim_model.theta_names = {'\theta_1','\theta_2'};
    sim_model.n_data = 1;
    sim_model.true_theta = [2,0.1];
    sim_model.dim = 2;
    
    %l95init % generates l95s0.mat
    load l95s0.mat % loads s0
    T = 160;
    sim_model.gen = @(theta,n) l95_par_run(T,theta,s0); % generate one data set
    sim_model.summary_dim = 6;
    sim_model.comp_summaries = @(data,obs_data) myStats(data);
    
    % data from the model and true posterior
    sim_model.data = sim_model.gen(sim_model.true_theta, []);
    sim_model.summary_true = sim_model.comp_summaries(sim_model.data,[]);
    sim_model.prior_eval = @(theta) ones(size(theta,1),1);
    sim_model.true_post_pdf2d = sl_baseline_get(name);
    
elseif strcmp(name,'ma2')
    % Moving average model, 2 parameters, 2 summaries, as in Marin et al. 2012 ABC review paper.
    % NOTE: We generate data from the model and use it to infer the parameters.
    
    grid_th = theta_grid([-2,2; -1,1],75); % could be the triangle
    
    % simulation model specification
    sim_model.name = 'MA2';
    sim_model.theta_names = {'\theta_1','\theta_2'};
    sim_model.n_data = 100;
    sim_model.true_theta = [0.6,0.2];
    sim_model.dim = 2;
    
    sim_model.gen = @(theta,n) simul_ma2(theta,n); 
    sim_model.summary_dim = 2;
    sim_model.comp_summaries = @(data,obs_data) summaries_ma2(data,[]);
    
    % data from the model and true posterior
    sim_model.data = sim_model.gen(sim_model.true_theta, sim_model.n_data);
    sim_model.summary_true = sim_model.comp_summaries(sim_model.data,[]);
    sim_model.prior_eval = @(theta) ones(size(theta,1),1); % could limit to the triangle
    sim_model.true_post_pdf2d = sl_baseline_get(name);
    
else
    error('Incorrect test model.');
end

rng(cur_seed);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper functions for the test simulation models:

function [post,samples] = sl_baseline_get(model_name)
% Loads the baseline MCMC-SL results computed elsewhere using extensive simulations. 
fn = ['../results/sl_baseline/',model_name,'_run1.mat'];
try
    load(fn,'post_kde','sl_samples');
    post = post_kde; samples = sl_samples;
catch
    %error('Baseline posterior could not be loaded.');
    warning('Baseline posterior could not be loaded.');
    post = []; samples = []; % if the precomputed baseline posterior is not found
end
end


function s = ricker_wrapper(theta,n,i,true_th,T)
theta_full = true_th;
theta_full(i) = theta;
s = simulate_ricker(theta_full,1,T);
end


function lp = log_simple_pdf(theta)
% Evaluates log pdf of a simple, Gaussian-like distribution which small correlation.

if size(theta,2) ~= 2
    error('Incorrect parameter input to the simple density.');
end
rho = 0.25;
S = [1 rho; rho 1];
invS = inv(S);
y = theta;
lp = -0.5*sum(y*invS.*y,2);
end

function lp = log_banana_pdf(theta,ab)
% Evaluates log pdf of a banana distribution as in Example 1 of http://helios.fmi.fi/~lainema/dram/
% Note: The example specification of the banana-function is not the same as was actually
% used in the code.

if size(theta,2) ~= 2
    error('Incorrect parameter input to the banana density.');
end
a = ab(1);
b = ab(2);
rho = 0.9;
S = [1 rho; rho 1];
invS = inv(S);
y = [theta(:,1)/a, theta(:,2)*a + a*b*(theta(:,1).^2 + a^2)];
lp = -0.5*sum(y*invS.*y,2);
end


function lp = log_bimodal_pdf(theta)
% Evaluates log pdf of a bimodal distribution.

if size(theta,2) ~= 2
    error('Incorrect parameter input to the bimodal density.');
end
rho = 0.5;
S = [1 rho; rho 1];
invS = inv(S);
y = [theta(:,1), theta(:,2).^2 - 2];
lp = -0.5*sum(y*invS.*y,2);
end


function lp = log_pdf_nd(theta,log_pdf,d)
% Evaluates log pdf of the multidim test problem which are constructed from the 2D test
% problems. 

lp = 0;
for i = 1:(d/2)
    ind = (i-1)*2+1;
    lp = lp + log_pdf(theta(:,ind:(ind+1)));
end
end

function [true_post_pdf,samples] = high_dim_toy_model_inference(grid_th, sim_model)
% Inference for obtaining baseline i.e. ground truth posterior for the multidim toy
% problems. In principle we could just use pseudo-marginal MCMC but as the toy problems
% are constructed from dim/2 independent 2D densities, we here take advantage of that.
%
% NOTE: The seed is always fixed so the Monte Carlo error is always the same.

% we sample exactly (under a grid assumption) from each 2d block
n = 10000;
d = sim_model.dim;
grid2dm = theta_grid(grid_th.range(1:2,:),100); % same for each 2d block
loglik2dm = sim_model.loglik_2d(grid2dm.theta2d');
lik2dm = exp(loglik2dm-max(loglik2dm));
samples = NaN(n,d);
for i = 1:(d/2)
    ind = (i-1)*2+1;
    sampled_ind = randsample(1:numel(lik2dm),n,true,lik2dm);
    samples(:,ind:(ind+1)) = grid2dm.theta2d(:,sampled_ind)'; 
    % TODO: could also add some jitter because of the grid approximation
end

% compute marginals in the grid
true_post_pdf = kde_for_abc(grid_th,samples)';
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions for MA2-model

function data = simul_ma2(theta,n)
% Simulates data from MA2-model.

data = zeros(n,1);
u = randn(n+2,1);
for i = 1:n
    data(i) = u(i+2) + theta(1)*u(i+1) + theta(2)*u(i);
end
end


function tau_j = compute_autocov(y,j)
% Compute the autocovariance tau_j = sum_(k=j+1)^n y_k y_(k-1) for the MA model.

n = length(y);
tau_j = 0;
for k = j + 1:n
    tau_j = tau_j + y(k)*y(k-j);
end
end


function s = summaries_ma2(data,obs_data)
% Compute the summary statistics for MA2-model.

tau1 = compute_autocov(data,1);
tau2 = compute_autocov(data,2);
s = [tau1;tau2];
end




