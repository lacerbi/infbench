function [y,y_std] = infbench_norton2019(x,infprob,mcmc_params)
%INFBENCH_NORTON2019 Inference benchmark log pdf -- changepoint detection model.

if nargin < 3; mcmc_params = []; end

problem_name = 'norton2019';
infbench_fun = str2func(['infbench_' problem_name]);

if isempty(x)
    if isempty(infprob) % Generate this document        
        fprintf('\n');

        % Add sampler directory to MATLAB path
        pathstr = fileparts(mfilename('fullpath'));
        addpath([pathstr,filesep(),'..',filesep(),'infalgo',filesep(),'parallel-GP-SL',filesep(),'utils',filesep(),'mcmcstat-master']);
                
        for n = 1     
            name = ['S' num2str(n)];
            
            infprob = infbench_fun([],n);
            infprob.DeterministicFlag = true;
            if isempty(mcmc_params); id = 0; else; id = mcmc_params(1); end

            D = infprob.D;
            trinfo = infprob.Data.trinfo;
            
            % Prior used for sampling
            infprob.PriorMean = infprob.Prior.Mean;
            infprob.PriorVar = diag(infprob.Prior.Cov)';
            infprob.PriorVolume = prod(infprob.UB - infprob.LB);
            infprob.PriorType = 'uniform';            

            LB = infprob.LB;
            UB = infprob.UB;
            PLB = infprob.PLB;
            PUB = infprob.PUB;
            
            if id == 0
            
                % Compute optimum
                Nopts = 5;
                
                opts = struct('Display','iter','MaxFunEvals',2e3);

                for iOpt = 1:Nopts                
                    x0 = rand(1,D).*(PUB - PLB) + PLB;
                    fun = @(x) -infbench_fun(x,infprob);
                    [xnew(iOpt,:),fvalnew(iOpt)] = bads(fun,x0,LB,UB,PLB,PUB,[],opts);
%                    [xnew(iOpt,:),fvalnew(iOpt)] = fmincon(fun,x0,[],[],[],[],LB,UB,[],opts);
                end
                
                [fvalnew,idx_best] = min(fvalnew);
                xnew = xnew(idx_best,:);
                
                fvalnew = -fvalnew;
                xmin = warpvars(xnew,'inv',trinfo);
                fval = fvalnew + warpvars(xnew,'logp',trinfo);

                x0 = xmin;
                x0 = warpvars(x0,'d',trinfo);   % Convert to unconstrained coordinates            
                fun = @(x) -logpost(x,infprob);
                [xnew,fvalnew] = bads(fun,x0,LB,UB,PLB,PUB,[],opts);
%                [xnew,fvalnew] = fmincon(fun,x0,[],[],[],[],LB,UB,[],opts);

                fvalnew = -fvalnew;
                xmin_post = xmin;
                xmin_post = warpvars(xnew,'inv',trinfo);
                fval_post = fvalnew + warpvars(xnew,'logp',trinfo);

                fprintf('\t\t\tcase %d\n',n);
                fprintf('\t\t\t\tname = ''%s'';\n\t\t\t\txmin = %s;\n\t\t\t\tfval = %s;\n',name,mat2str(xmin),mat2str(fval));
                fprintf('\t\t\t\txmin_post = %s;\n\t\t\t\tfval_post = %s;\n',mat2str(xmin_post),mat2str(fval_post));
                
            elseif id > 0 && n == mcmc_params(2)

                rng(id);
                widths = 0.5*(PUB - PLB);
                logpfun = @(x) logpost(x,infprob);
                
                % Number of samples
                if numel(mcmc_params) > 2
                    Ns = mcmc_params(3);
                else
                    Ns = 1e3;
                end
                
                W = 2*(infprob.D+1);    % Number of walkers
                
                sampleopts.Thin = 17;
                sampleopts.Burnin = Ns*sampleopts.Thin;
                sampleopts.Display = 'notify';
                sampleopts.Diagnostics = false;
                sampleopts.VarTransform = false;
                sampleopts.InversionSample = false;
                sampleopts.FitGMM = false;
                sampleopts.TolX = 1e-5;
                % sampleopts.TransitionOperators = {'transSliceSampleRD'};

                x0 = infprob.Post.Mode;
                x0 = warpvars(x0,'d',trinfo);
                
                [Xs,lls,exitflag,output] = eissample_lite(logpfun,x0,Ns,W,widths,LB,UB,sampleopts);
                
%                 % Re-evaluate final lls with higher precision
%                 fprintf('Re-evaluate log posteriors...\n');
%                 infprob.Ng = 501;
%                 lls = zeros(size(Xs,1),1);
%                 for i = 1:size(Xs,1)
%                     lls(i) = logpost(Xs(i,:),infprob);
%                     if mod(i,ceil(size(Xs,1)/10))==0; fprintf('%d/%d..',i,size(Xs,1)); end
%                 end
%                 fprintf('\nDone.\n');
                
                
                filename = [problem_name '_mcmc_n' num2str(n) '_id' num2str(id) '.mat'];
                save(filename,'Xs','lls','exitflag','output');                
            end
            
        end
        
    else
        % Initialization call -- define problem and set up data
        n = infprob(1);
        
        % Are we transforming the entire problem to unconstrained space?
        transform_to_unconstrained_coordinates = false;
        
        % Add problem directory to MATLAB path
        pathstr = fileparts(mfilename('fullpath'));
        addpath([pathstr,filesep(),problem_name]);
        
        switch n
			case {1,101}
                nid = 5;
                modelIndex = 12;
				xmin = [1.46277127696673 -2.1681896767972 0.0644547353054086 -2.63903192475833 2.98550740964576 0.527213339997871 -0.486009929908597 1.44165159656731 -1.09762755985263];
				fval = -6119.75845967798;
				xmin_post = [1.46277115993661 -2.16818934149361 0.0644541883964579 -2.63903178693316 2.98550739097048 0.527212363432901 -0.486009373135257 1.441651287652 -1.0976269939063];
				fval_post = -6135.88429490272;       
                % R_max = 1.001. Ntot = 100000. Neff_min = 98320.1. Total funccount = 9600356.
                Mean_mcmc = [1.46435616822208 -2.17136309719058 0.0662371968199702 -0.334478257778285 0.680593825602192 0.527985463809655 -0.487041253601341 0.361789560357904 -1.10430955104543];
                Cov_mcmc = [0.00264442949861135 -0.00227324311869048 2.71619884653618e-05 -0.000115242744900933 0.00021988909418795 0.00012133224919512 -0.000135268063605132 0.000113173308505702 -0.00111293895938135;-0.00227324311869048 0.00385308530641299 -0.000114417511142836 -0.000836605351275131 0.000659813611845422 -0.00035512382197266 0.000331776148388622 0.000232540708263256 0.00191824949333633;2.71619884653618e-05 -0.000114417511142836 0.000562234083187425 0.00141007313249507 -0.00140509847545958 7.47780222149517e-05 -8.45600609635136e-05 -0.000628820759690247 -0.000580051167672241;-0.000115242744900933 -0.000836605351275131 0.00141007313249507 2.29055616793922 -2.28984729970018 0.000173538096446579 -0.00030630980873009 -1.07337857780648 -0.00502169644942798;0.00021988909418795 0.000659813611845422 -0.00140509847545958 -2.28984729970018 2.29094259326485 -0.000102674868697585 0.000232214833236278 1.07341045776547 0.0047958485563488;0.00012133224919512 -0.00035512382197266 7.47780222149517e-05 0.000173538096446579 -0.000102674868697585 0.00300124182724716 -0.00225836300475242 1.04698866459393e-05 -0.00161661488501393;-0.000135268063605132 0.000331776148388622 -8.45600609635136e-05 -0.00030630980873009 0.000232214833236278 -0.00225836300475242 0.00331194970939271 5.97198826346781e-05 0.0014688542168854;0.000113173308505702 0.000232540708263256 -0.000628820759690247 -1.07337857780648 1.07341045776547 1.04698866459393e-05 5.97198826346781e-05 0.503656767723943 0.00161252474960342;-0.00111293895938135 0.00191824949333633 -0.000580051167672241 -0.00502169644942798 0.0047958485563488 -0.00161661488501393 0.0014688542168854 0.00161252474960342 0.00661121593996541];
                lnZ_mcmc = -6152.2651430423;                                
        end
        
        % Get model and preprocess data        
        fitinfo = norton2019_getinfo(nid,modelIndex);
        
        % Define parameter upper/lower bounds
        lb = fitinfo.LB;
        ub = fitinfo.UB;
        plb = fitinfo.LB + 0.1*(fitinfo.UB - fitinfo.LB);
        pub = fitinfo.LB + 0.9*(fitinfo.UB - fitinfo.LB);
        noise = [];
        
        D = numel(lb);
        xmin = NaN(1,D);       fval = Inf;
        xmin_post = NaN(1,D);  fval_post = Inf;
        Mean_mcmc = NaN(1,D);  Cov_mcmc = NaN(D,D);    lnZ_mcmc = NaN;
        
        if transform_to_unconstrained_coordinates
            trinfo = warpvars(D,lb,ub,plb,pub);     % Transform to unconstrained space
            trinfo.mu = zeros(1,D);     % Necessary for retro-compatibility
            trinfo.delta = ones(1,D);
        else
            trinfo = [];
        end
                
        data.Ntrials = size(x,1);
        data.lapse = 0.02;  % Fix lapse rate
        
        xmin = warpvars(xmin,'d',trinfo);
        fval = fval + warpvars(xmin,'logp',trinfo);
        
        Mean = zeros(1,D);
        Cov = eye(D);
        Mode = xmin;
                
        y.D = D;
        y.LB = warpvars(lb,'d',trinfo);
        y.UB = warpvars(ub,'d',trinfo);
        y.PLB = warpvars(plb,'d',trinfo);
        y.PUB = warpvars(pub,'d',trinfo);
        
        y.lnZ = 0;              % Log normalization factor
        y.Mean = Mean;          % Distribution moments
        y.Cov = Cov;
        y.Mode = Mode;          % Mode of the pdf
        y.ModeFval = fval;
        
        priorMean = 0.5*(y.PUB + y.PLB);
        priorSigma2 = (0.5*(y.PUB - y.PLB)).^2;
        priorCov = diag(priorSigma2);
        y.Prior.Mean = priorMean;
        y.Prior.Cov = priorCov;
        
        % Compute each coordinate separately
        xmin_post = warpvars(xmin_post,'d',trinfo);
        fval_post = fval_post + warpvars(xmin_post,'logp',trinfo);
                
        y.Post.Mean = Mean_mcmc;
        y.Post.Mode = xmin_post;          % Mode of the posterior
        y.Post.ModeFval = fval_post;        
        y.Post.lnZ = lnZ_mcmc;
        y.Post.Cov = Cov_mcmc;
                
        if n > 100
            fitinfo.IBSNreps = 0; % Deterministic problems            
        else
            fitinfo.IBSNreps = 500; % Reps used for IBS estimator
        end
        
        % Read marginals from file
%        marginals = load([problem_name '_marginals.mat']);
%        y.Post.MarginalBounds = marginals.MarginalBounds{nid};
%        y.Post.MarginalPdf = marginals.MarginalPdf{nid};
        
        % Save data and coordinate transformation struct
        fitinfo.trinfo = trinfo;
        y.Data = fitinfo;
        y.DeterministicFlag = (fitinfo.IBSNreps == 0);
                
    end
    
else
    
    % Iteration call -- evaluate objective function
    
    % Transform unconstrained variables to original space
    x_orig = warpvars(x,'i',infprob.Data.trinfo);
    dy = warpvars(x,'logpdf',infprob.Data.trinfo);   % Jacobian correction
    
    % Compute log likelihood of data and possibly std of log likelihood
    if infprob.DeterministicFlag        
        LL = changeprob_nll(x_orig,infprob.Data);
        y_std = 0;
    else
        
        
        
        Nibs = infprob.Data.IBSNreps/50;
        IBSNreps = 50; % Split to avoid memory errors (should fix ibslike)
        ibs_opts = struct('Nreps',IBSNreps,...
            'ReturnPositive',true,'ReturnStd',false);        
        for iter = 1:Nibs
            [LL(iter),y_var(iter)] = ibslike(@akrami2018_gendata,x_orig,infprob.Data.y,[],ibs_opts,infprob.Data);
        end
        LL = mean(LL);
        y_std = sqrt(mean(y_var)/Nibs);
    end
    y = LL + dy;
    
end

end

%--------------------------------------------------------------------------
function y = logpost(x,infprob)    
    y = infbench_norton2019(x,infprob);
    lnp = infbench_lnprior(x,infprob);
    y = y + lnp;
end

function ll = akrami2018_llfun(theta,data)

% Model parameters
w = theta(1:end);
lambda = data.lapse;

p_left = 1./(1+exp(sum(bsxfun(@times,w,data.X),2)));
p_lresp = (data.y == 1).*p_left + (data.y == 2).*(1-p_left);
p_lresp = (1-lambda)*p_lresp + lambda/2;

ll = sum(log(p_lresp));

end

function R = akrami2018_gendata(theta,idx,data)

% Model parameters
w = theta(1:end);
lambda = data.lapse;

p_left = 1./(1+exp(sum(bsxfun(@times,w,data.X),2)));
p_lresp = (1-lambda)*p_left + lambda/2;

pl_vec = p_lresp(idx);

R = double(rand(size(pl_vec)) < pl_vec);
R(R == 0) = 2;

end

function akrami2018_test(infprob)

theta = randn(1,9);
LL = akrami2018_llfun(theta,infprob.Data);
ibs_opts = struct('Nreps',infprob.Data.IBSNreps,...
    'ReturnPositive',true,'ReturnStd',true);
[LL_ibs,y_std] = ibslike(@akrami2018_gendata,theta,infprob.Data.y,[],ibs_opts,infprob.Data);

[LL,LL_ibs,y_std]
[LL - LL_ibs, y_std]

end

%--------------------------------------------------------------------------
function fitinfo = norton2019_getinfo(subIndex,modelIndex)
%RUNFIT Runs model comparison for changing probability experiment

subID = {'CWG','EGC','EHN','ERK','GK','HHL','JKT','JYZ','RND','SML','SQC'};

% Models to fit (only Bayesian models available)
models = {[],'idealBayesian',[],[],[],[],[], ...
    'subBayesian_rlprior','subBayesian_conservative','subBayesian_pVec',...
    'subBayesian_betahyp','subBayesian_3param','subBayesian_4param','subBayesian_5param','subBayesian_flex',...
    [],[],[]};

runSubject = subID{subIndex};
runModel = models{modelIndex};
task = 2;       % Covert task

% Load data
temp = load(['ChangingProbabilities_', runSubject]);
data = temp.data;

MaxParams = 15;
switch(runModel)
    case 'subBayesian_rlprior'
        parameters = [1 0 0 0 0 0 1 0 0 0 0 0 0];
    case 'subBayesian_conservative'
        parameters = [1 0 0 0 0 1 0 0 0 0 0 0 0];
    case 'subBayesian_pVec'
        parameters = [1 0 0 0 0 0 0 1 0 0 0 0 0];
    case 'subBayesian_betahyp'
        parameters = [1 0 0 0 0 0 0 0 1 0 0 0 0];
    case 'subBayesian_3param'
        parameters = [1 0 0 0 0 0 1 0 1 0 0 0 0];
    case 'subBayesian_flex'
        parameters = [1 0 0 0 0 0 1 1 1 0 0 0 0 1 1];
    case 'subBayesian_4param'
        parameters = [1 0 0 0 0 0 1 1 1 0 0 0 0 0 0];
    case 'subBayesian_5param'
        parameters = [1 0 0 0 0 0 1 1 1 0 0 0 0 0 1];
    otherwise
        parameters = zeros(1,MaxParams);
end
parameters(1) = 1;
parameters = [parameters,zeros(1,MaxParams-numel(parameters))];

runModel = 'idealBayesian';

fitinfo.task = task;
fitinfo.runModel = runModel;
fitinfo.parameters = parameters;


% Parameters to be fit
paramNames = {'sigma_ellipse','sigma_criterion','lambda','gamma','alpha','w','Tmax', ...
    'pVec','beta','delta1','delta2','hRate','nu_p','delta_Tmax','delta_pVec'};
I_params = find(parameters);

% Lower and upper parameter bounds
paramBounds_def = [1,30; 1,30; 0,0.1; -Inf,Inf; 0,1; 0,1; 2,200; ...
    0.01,0.49; 0,10; 1.01,5; 1.01,14; 0,1; 0,5; 1,200; 0.01,0.49];    
paramBounds = paramBounds_def(I_params,:);

% Integer parameters (used only by VBMC)
%paramInteger = false(1,numel(paramNames));
%paramInteger([7,14]) = true;
%paramInteger = paramInteger(I_params);

%% Get session parameters

[NumTrials,sigma_ellipseData,mu,sigma_s,C,S,p_true,resp_obs,score] = ...
    changeprob_getSessionParameters(data,task);

fitinfo.NumTrials = NumTrials;
fitinfo.sigma_ellipseData = sigma_ellipseData;
fitinfo.mu = mu;
fitinfo.sigma_s = sigma_s;
fitinfo.C = C;
fitinfo.S = S;
fitinfo.p_true = p_true;
fitinfo.resp_obs = resp_obs;
fitinfo.score = score;
fitinfo.I_params = I_params;

X = [];

% Default parameter values for those not fitted
inputParams = zeros(1,numel(parameters));
I_notFit = find(parameters == 0);
sigmacriterion_def = 5; % Default sigma_criterion
lapse_def = 1e-4;       % Default lapse (i.e., tiny lapse)
gamma_def = Inf;        % Default gamma (i.e., BDT)
alpha_def = 0.2;        % Default alpha
w_def     = 1;          % Default w (i.e., no bias)
Tmax_def  = 0;          % Use default prior window
pVec_def = 0;           % Use default probability vector
beta_def = 0;           % Use default hyperprior, [0,0]
delta_def = 2;          % Use default node distance
hRate_def = .01;        % Use default hazard rate (average rate of change)
nu_p_def = log(2);      % Use default nu_p (Beta(1,1))
delta_Tmax_def = 0;     % Use default delta_Tmax
delta_pVec_def = 0;     % Use default delta_pVec

notFit_def = [sigma_ellipseData, sigmacriterion_def, lapse_def, gamma_def, ...
    alpha_def, w_def, Tmax_def, pVec_def, beta_def, delta_def, delta_def, ...
    hRate_def, nu_p_def, delta_Tmax_def, delta_pVec_def];

inputParams(I_notFit) = notFit_def(I_notFit);

fitinfo.inputParams = inputParams;
fitinfo.LB = paramBounds(:,1)';
fitinfo.UB = paramBounds(:,2)';
%inputParams
    
%% Recompute additional outputs from best fit params

%inputParams(I_params) = fitParams;
% sigma = sqrt(sigma_s^2 + inputParams(1)^2);
%[nLL,~,resp_model,p_estimate,~] = nLL_fun(fitParams);

end

function [nLL,resp_model,p_estimate] = changeprob_nll(params,fitinfo)
%CHANGEPROB_NLL Negative log likelihood for the specified model

idx_params = fitinfo.I_params;
inputParams = fitinfo.inputParams;
NumTrials = fitinfo.NumTrials;
mu = fitinfo.mu;
sigma_s = fitinfo.sigma_s;
C = fitinfo.C;
S = fitinfo.S;
p_true = fitinfo.p_true;
resp_obs = fitinfo.resp_obs;
task = fitinfo.task;
score = fitinfo.score;

inputParams(idx_params) = params;
sigma = sqrt(sigma_s^2 + inputParams(1)^2);

if inputParams(7) == 0
    prior_rl = [];
elseif (inputParams(7) ~= 0) && (inputParams(14) == 0)
    prior_rl = [max(1,inputParams(7)*2/3),inputParams(7)];
else
    prior_rl = [(inputParams(7)-1), inputParams(7)-1+inputParams(14)];
end
if inputParams(8) == 0
    p_vec = [];
elseif (inputParams(8) == 1) && (inputParams(15) == 0)
    p_vec = linspace(inputParams(8), 1-inputParams(8), 5);
else
    p_vec = linspace(inputParams(8), 0.5+inputParams(15), 5);
end
if inputParams(9) == 0
    beta_hyp = [];
else
    beta_hyp = inputParams(9)^2;
end
[nLL,rmse,p_estimate,resp_model] = ...
    ChangeProb_bocpd_nll_v2(inputParams(1:6), ...
    NumTrials,mu,sigma,C,S,p_true,resp_obs,score,task,prior_rl,p_vec,beta_hyp);


end

function [NumTrials,sigma_ellipseData,mu,sigma_s,C,S,p_true,resp,score] = changeprob_getSessionParameters(data,task)
%CHANGEPROB_GETSESSIONPARAMETERS Gets session parameters from an existing
%data struct or creates a fake dataset
        
sigma_ellipse = [];
sigma_criterion = [];

col = data.SessionOrder(task);     % Column of overt-criterion task
NumTrials = data.NumTrials;    

% Noise parameters
sigma_s = data.StdDev;
if isempty(sigma_ellipse) || ~isfinite(sigma_ellipse)
    sigma_ellipseData = data.EllipseNoise;
    sigma = sqrt(sigma_s^2 + sigma_ellipseData^2);
else
    sigma = sqrt(sigma_s^2 + sigma_ellipse^2);
end

% Category information 
% (in the data Category B/Red is coded as 1, Category A/Green is coded as 2)
C = (data.Category(:,col) == 2);        % Category A/Green
C = double(C);
p_true = data.pA(:,col);
mu = [data.GreenMean(col),data.RedMean(col)];

% Shift coordinate system to zero
mu_bar = mean(mu);
mu = bsxfun(@minus, mu, mu_bar);
S = bsxfun(@minus, data.StimulusAngle(:,col), mu_bar);

% Get task-relevant responses
switch task
    case 1  % Overt-criterion task    
        resp = bsxfun(@minus, data.Criterion(:,col), mu_bar);   % Reported criterion
    case 2  % Covert-criterion task
        resp = data.Response(:,col) == 2;
        resp = double(resp);
end
score = data.Score(:,col);
    
end


function [nLL,rmse,p_estimate,resp_model,post] = ChangeProb_bocpd_nll_v2(parameters, NumTrials, mu, sigma, C, S, p_true, resp_obs, score, task, prior_rl, p_vec, beta_hyp)
%CHANGEPROB_BOCPD_NLL Bayesian online changepoint detection observer.

% Parameter vector:
% #1 is SIGMA_ELLIPSE, #2 is SIGMA_CRITERION, #3 is LAPSE, #4 is GAMMA,
% #5 is ALPHA, and #6 is W
if nargin < 1
    parameters = [];
    prior_rl = [];
    p_vec = [];
    [NumTrials, sigma_ellipse, mu, sigma, C, S, p_true, resp_obs] = changeprob_getSessionParameters();
    task = 1; 
end

% Data struct or random seed for fake data generation
if nargin < 2; error('You must specify the session parameters.'); end

% Task (1 overt, 2 covert, 3 mixed)
if nargin < 10 || isempty(task); task = 1; end
if task ~= 1 && task ~= 2 && task ~=3; error('TASK can only be 1 (overt-criterion), 2 (covert-criterion), or 3 (mixed).'); end

if nargin < 11 || isempty(prior_rl)
    prior_rl = [80,120];    % Default runlengths
end

if nargin < 12 || isempty(p_vec)
    p_vec = linspace(0.2,0.8,5); % Default states
end

if nargin < 13 || isempty(beta_hyp)
    beta_hyp = [0,0];   % Default beta hyperprior after changepoint (maximum-likelihood)
end
if isscalar(beta_hyp); beta_hyp = beta_hyp*[1,1]; end

%% Experiment constants

% Run length ~ Uniform[prior_rl(1),prior_rl(2)]
% if isfield(options,'prior_rl') && ~isempty(options.prior_rl)
%     prior_rl = options.prior_rl;
% else
%end
runlength_min = prior_rl(1);
runlength_max = prior_rl(2);

% Probability states
% if isfield(options,'p_vec') && ~isempty(options.p_vec)
%     p_vec = options.p_vec;
% else
%     p_vec = linspace(0.2,0.8,5);    % Default states
%end
Nprobs = numel(p_vec);              % # states

p_bias = .5;

if task == 3
    idx_Overt = 5:5:NumTrials;
    idx_Covert = 1:NumTrials;
    idx_Covert(idx_Overt) = [];
end

%% Observer model parameters

% Default values
lambda = 0;
gamma = Inf;
w = 1;

switch numel(parameters)
    case 0  % Empty parameter vector
        sigma_criterion = sigma_ellipse;
    case 1
        sigma_ellipse = parameters(1);
        sigma_criterion = sigma_ellipse;
    case 2
        sigma_ellipse = parameters(1);
        sigma_criterion = parameters(2);        
    case 3
        sigma_ellipse = parameters(1);
        sigma_criterion = parameters(2);
        lambda = parameters(3);
    case 4
        sigma_ellipse = parameters(1);
        sigma_criterion = parameters(2);
        lambda = parameters(3);
        gamma = parameters(4);
    case 5
        sigma_ellipse = parameters(1);
        sigma_criterion = parameters(2);
        lambda = parameters(3);
        gamma = parameters(4);
    case 6
        sigma_ellipse = parameters(1);
        sigma_criterion = parameters(2);
        lambda = parameters(3);
        gamma = parameters(4);
        w = parameters(6);
    otherwise
        error('PARAMETERS should be a vector with up to six parameters.');
end

%% Initialize inference

% Hazard function (defined only where nonzero)
L = (ceil(runlength_max)-floor(runlength_min))+1;  % Number of elements
rlprior = ones(L,1)/L;      % Constant changepoint prior

% Trick to allow non-integer run length min and max
frac_min = 1 - (runlength_min - floor(runlength_min));
rlprior(1) = rlprior(1)*frac_min;
frac_max = 1 - (ceil(runlength_max) - runlength_max);
rlprior(end) = rlprior(end)*frac_max;
H = rlprior(:)./flipud(cumsum(flipud(rlprior(:))));

% Posterior over run lengths (from 0 to PRIOR_RL(2)-1)
post = zeros(ceil(runlength_max),Nprobs);
post(1,:) = 1;  % Change in the first trial
post = post ./ sum(post(:));

% Table of binomial count/probability
Psi = zeros(size(post,1),1,Nprobs); Psi(1,1,:) = 1/Nprobs;

% Transition matrix
Tmat(1,:,:) = (ones(Nprobs,Nprobs) - eye(Nprobs)) / (Nprobs-1);
% Tmat(1,:,:) = 1;
p_vec3(1,1,:) = p_vec;

% Auxiliary variables for covert-criterion task
if task == 2 || task == 3
    nx = 1001;  % Measurement grid size
    MAXSD = 5;  % When integrating a Gaussian go up to this distance
    X = bsxfun(@plus, S, sigma_ellipse*MAXSD*linspace(-1,1,nx));
    W = normpdf(linspace(-MAXSD,MAXSD,nx));
    W = W./qtrapz(W);
end

%% Begin loop over trials
P = zeros(NumTrials+1,Nprobs);
P(1,:) = ones(1,Nprobs)/Nprobs;
last = zeros(NumTrials,size(post,1));
PCx = [];

for t = 1:NumTrials
    %t
    if task == 2 || and(task == 3, mod(t,5) ~= 0); Xt = X(t,:); else Xt = []; end    
    [post,Psi,pi_post,PCxA] = bayesianOCPDupdate(Xt,C(t),post,Psi,Tmat,H,p_vec3,beta_hyp,mu,sigma,task,t);
    tt = nansum(post,2);

    % The predictive posterior is about the next trial
    P(t+1,:) = pi_post;
    last(t,:) = tt/sum(tt);
    
    % Record conditional posterior for covert-criterion task
    if task == 2 || and(task == 3, mod(t,5) ~=0)
        if isempty(PCx); PCx = zeros(NumTrials,size(PCxA,2)); end
        PCx(t,:) = PCxA;
    end
end

%% Compute log likelihood

MIN_P = 1e-4;   % Minimum lapse/error probability

% Compute predicted criterion
pA_t = sum(bsxfun(@times, P, p_vec(:)'),2);
pB_t = sum(bsxfun(@times, P, 1-p_vec(:)'), 2);
pA_t_bias = bsxfun(@plus, w*pA_t, (1-w)*p_bias);
pB_t_bias = bsxfun(@plus, w*pB_t, (1-w)*p_bias);
Gamma_t = pA_t_bias./pB_t_bias;
z_opt = sigma^2 * log(Gamma_t) / diff(mu);
z_opt = z_opt(1:NumTrials);
        
switch task        
    case 2  % Covert-criterion task
        
        Pcx = 1./(1 + ((1-PCx)./PCx).^gamma);       % Softmax        
        PChatA = qtrapz(bsxfun(@times,Pcx,W),2);    % Marginalize over noise
        PChatA_bias = w*PChatA + (1-w)*p_bias;      % Conservative bias when w < 1
        lambda = max(1e-4,lambda);                  % Minimum lapse to avoid numerical trouble
        PChatA_bias = lambda/2 + (1-lambda)*PChatA_bias;
        
        % Log probability of covert task responses
        log_PChat = log(PChatA_bias).*(resp_obs == 1) + log(1-PChatA_bias).*(resp_obs ~= 1);
        
        % Sum negative log likelihood
        log_PChat(~isfinite(log_PChat)) = log(MIN_P);
        nLL = -nansum(log_PChat);   
        resp_model = PChatA_bias(1:NumTrials);
end

% RMSE between predictive posterior probability and true category probability
meanP = sum(bsxfun(@times,p_vec,P(1:NumTrials,:)),2);
p_estimate = meanP;
rmse = sqrt(mean((meanP - p_true).^2));

end

%--------------------------------------------------------------------------
function [post,Psi,pi_post,PCxA] = bayesianOCPDupdate(X,C,post,Psi,Tmat,H,p_vec,beta_hyp,mu,sigma,task,trial)
%BAYESIANCPDUPDATE Bayesian online changepoint detection update

    %if mod(t,100) == 0
    %    t
    %end
    L = size(H,1);  % Width of prior over run lengths
    % Slice with nonzero hazard function (probability of change)
    idxrange = (size(post,1)-L+1:size(post,1));

    % Posterior over pi_i before observing C_t
    pi_post = bsxfun(@times, Psi, Tmat);
    predCatA(:,:) = sum(bsxfun(@times, pi_post, p_vec),3);
    predCatB(:,:) = sum(bsxfun(@times, pi_post, 1 - p_vec),3);
    
    %----------------------------------------------------------------------
    % 2. Compute probability of response for covert-criterion task
    if task == 2 || (task == 3 && mod(trial,5) ~= 0)
        pxCatA = exp(-0.5*((X-mu(1))./sigma).^2);
        pxCatB = exp(-0.5*((X-mu(2))./sigma).^2);

        PCxA = nansum(nansum(bsxfun(@times,predCatA,post),2),1).*pxCatA;
        PCxB = nansum(nansum(bsxfun(@times,predCatB,post),2),1).*pxCatB;    
        PCxA = PCxA./(PCxA + PCxB);
    else
        PCxA = [];
    end
    
    %----------------------------------------------------------------------
    % 3. Observe Ct (do nothing)
    
    %----------------------------------------------------------------------
    % 4a. Evaluate predictive probability
    if C == 1
        predCat = predCatA ./ (predCatA + predCatB);
    else
        predCat = predCatB ./ (predCatA + predCatB);         
    end
    
    % 4b. Multiply posterior by predictive probability
    post = bsxfun(@times, post, predCat);
        
    % 4c. Evaluate posterior probability over state (only for relevant range)
    if C == 1
        pi_postC = bsxfun(@times, pi_post(idxrange,:,:), p_vec);
    else
        pi_postC = bsxfun(@times, pi_post(idxrange,:,:), 1-p_vec);
    end
    pi_postC = bsxfun(@rdivide, pi_postC, sum(pi_postC,3));

    %----------------------------------------------------------------------    
    % 5a. Calculate unnormalized changepoint probabilities
    slice = post(idxrange,:);   % Slice out trials where change is possible
    currenttrial = nansum( ...
        bsxfun(@times, sum(bsxfun(@times, slice, pi_postC),2), H), ...
        1);
    
    % 5b. Calculate unnormalized growth probabilities    
    post(idxrange,:) = bsxfun(@times, slice, 1-H); 
    
    % Shift posterior by one step
    post = circshift(post,1,1);
    post(1,:) = currenttrial;
    post(isnan(post)) = 0;
    
    % 5c. Calculate normalization
    Z = sum(post(:));
    
    % 5d. Normalize posterior
    post = post ./ Z;
    
    %----------------------------------------------------------------------
    % 6a. Update sufficient statistics
    Psi = circshift(Psi,1,1);
    if C == 1
        Psi = bsxfun(@times, Psi, p_vec);
    else
        Psi = bsxfun(@times, Psi, 1 - p_vec);
    end
    
    % Hyperprior
    Psi(1,1,:) = exp(beta_hyp(1).*log(p_vec) + beta_hyp(2).*log(1-p_vec));
    Psi = bsxfun(@rdivide, Psi, sum(Psi,3));
    
    % 6b. Store predictive posterior over pi_t
    pi_post = nansum(sum(bsxfun(@times,bsxfun(@times, Psi, Tmat), post),2),1);
    pi_post = pi_post / sum(pi_post);

end

