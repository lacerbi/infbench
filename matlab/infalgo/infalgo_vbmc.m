function [history,post,algoptions] = infalgo_vbmc(algo,algoset,probstruct)

algoptions = vbmc('all');                   % Get default settings

ControlRunFlag = false;     % Do NOT run in control mode

% VBMC old defaults -- some of these may have changed
algoptions.FunEvalsPerIter = 5;
algoptions.AcqFcn = '@vbmc_acqskl';

algoptions.MinFunEvals = probstruct.MaxFunEvals;
algoptions.MaxFunEvals = probstruct.MaxFunEvals;
algoptions.MinIter = 0;     % No limits on iterations
algoptions.MaxIter = Inf;
algoptions.WarpNonlinear = 'off';   % No nonlinear warping for now
algoptions.BestSafeSD = 5;
algoptions.BestFracBack = 0.25;
algoptions.Diagnostics = 'on';
algoptions.InitDesign = 'plausible';    % Initial design uniform in plausible box
algoptions.EmpiricalGPPrior = 'yes';
algoptions.WarmupNoImproThreshold = Inf; 
algoptions.TolStableExceptions = 1;
algoptions.TolStableIters = 8;
algoptions.WarmupCheckMax = false;
algoptions.SGDStepSize = 0.01;
algoptions.RankCriterion = false;
algoptions.EmpiricalGPPrior = true;
algoptions.gpQuadraticMeanBound = false;
algoptions.WarmupOptions = [];

if probstruct.Debug
    algoptions.TrueMean = probstruct.Post.Mean;
    algoptions.TrueCov = probstruct.Post.Cov;
end

% Use prior as proposal function
algoptions.ProposalFcn = @(X_) exp(infbench_lnprior(X_,probstruct));

% Default options for variational active sampling
if numel(algoset) >= 3 && strcmpi(algoset(1:3),'vas')
    algoptions.VarActiveSample = true;
    algoptions.VariableMeans = false;
    algoptions.Plot = 'on';
    algoptions.NSent = 0;
    algoptions.NSentFast = 0;
    algoptions.NSentFine = '@(K) 2^15*round(sqrt(K))';
    algoptions.DetEntTolOpt = 1e-3;
    algoptions.EntropySwitch = true;
    algoptions.DetEntropyMinD = 0;
    algoptions.EntropyForceSwitch = Inf;
    algoptions.TolWeight = 0;
    algoptions.NSelbo = '@(K) 50*sqrt(K)';
    algoptions.SearchAcqFcn = @acqvasreg_vbmc;
    algoptions.NSsearch = 2^7;
    algoptions.Warmup = false;
    algoptions.FunEvalsPerIter = 1;
end

% Options from current problem
switch algoset
    case {0,'debug'}; algoset = 'debug'; algoptions.Debug = true; algoptions.Plot = 'on'; algoptions.FeatureTest = true;
    case {1,'base'}; algoset = 'base';                                                      % Use defaults
    case {2,'acqusreg'}; algoset = 'acqusreg'; algoptions.SearchAcqFcn = @acqusreg_vbmc;    % Vanilla uncertainty sampling
    case {3,'acqproreg'}; algoset = 'acqproreg'; algoptions.SearchAcqFcn = @acqfreg_vbmc;   % Prospective uncertainty sampling
    case {4,'control'}; algoset = 'control'; ControlRunFlag = true;                         % Control experiment
    case {5,'test'}; algoset = 'test'; algoptions.FeatureTest = true;                       % Test feature
    case {6,'narrow'}; algoset = 'narrow'; algoptions.InitDesign = 'narrow';                % Narrow initialization
    case {7,'control2'}; algoset = 'control2'; ControlRunFlag = true;                       % Control experiment, repeated
    case {8,'test2'}; algoset = 'test2'; algoptions.FeatureTest = true;                     % Test feature (second case)

    % Fixed number of mixture components
    case {11,'K1'}; algoset = 'K1'; algoptions.Kfun = 1; algoptions.KfunMax = 1; algoptions.Kwarmup = 1;
    case {12,'K2'}; algoset = 'K2'; algoptions.Kfun = 2; algoptions.KfunMax = 2; algoptions.Kwarmup = 2;
    case {13,'K5'}; algoset = 'K5'; algoptions.Kfun = 5; algoptions.KfunMax = 5; algoptions.Kwarmup = 5;
        
    % Ongoing research and testing
    case {21,'acqfregvlnn'}; algoset = 'acqfregvlnn'; algoptions.SearchAcqFcn = @vbmc_acqfregvlnn;
    case {22,'acqfregvsqrtn'}; algoset = 'acqfregvsqrtn'; algoptions.SearchAcqFcn = @vbmc_acqfregvsqrtn;
    case {23,'acqfregt'}; algoset = 'acqfregt'; algoptions.SearchAcqFcn = @vbmc_acqfregt;   % annealed prospective uncertainty sampling
    case {24,'se'}; algoset = 'se'; algoptions.gpMeanFun = 'se';           
    case {25,'const'}; algoset = 'const'; algoptions.gpMeanFun = 'const';        
    case {26,'acqf2reg'}; algoset = 'acqf2reg'; algoptions.SearchAcqFcn = @vbmc_acqf2reg;
    case {27,'acqpropnew'}; algoset = 'acqpropnew'; algoptions.SearchAcqFcn = @acqpropreg_vbmc; algoptions.gpQuadraticMeanBound = 1; algoptions.EmpiricalGPPrior = 0; algoptions.WarmupNoImproThreshold = 20 + 5*numel(probstruct.InitPoint); algoptions.TolStableExceptions = 2; algoptions.TolStableIters = 10; algoptions.WarmupCheckMax = true; algoptions.SGDStepSize = 0.005;
    case {28,'bowarmup'}; algoset = 'bowarmup';
        w.SearchCacheFrac = 0.1; w.HPDSearchFrac = 0.9; w.HeavyTailSearchFrac = 0; w.MVNSearchFrac = 0; w.SearchAcqFcn = @acqpropreg_vbmc; w.StopWarmupThresh = 0.1; w.SearchCMAESVPInit = false;
        algoptions.WarmupOptions = w; algoptions.Plot = 0;
        algoptions.TolStableWarmup = 55; algoptions.BOWarmup = true; algoptions.NSgpMaxWarmup = 8;
        algoptions.gpQuadraticMeanBound = 1; algoptions.EmpiricalGPPrior = 0; algoptions.WarmupNoImproThreshold = 220 + 5*numel(probstruct.InitPoint); algoptions.TolStableExceptions = 2; algoptions.TolStableIters = 10; algoptions.WarmupCheckMax = true; algoptions.SGDStepSize = 0.005;        
    case {29,'acqfheavy'}; algoset = 'acqfheavy'; algoptions.Plot = 1; algoptions.SearchAcqFcn = @acqfheavyreg_vbmc; algoptions.gpQuadraticMeanBound = 1; algoptions.EmpiricalGPPrior = 0; algoptions.WarmupNoImproThreshold = 20 + 5*numel(probstruct.InitPoint); algoptions.TolStableExceptions = 2; algoptions.TolStableIters = 10; algoptions.WarmupCheckMax = true; algoptions.SGDStepSize = 0.005;
    case {30,'gpswitch'}; algoset = 'gpswitch'; algoptions.Plot = 0; algoptions.StableGPSampling = '100+10*nvars'; algoptions.WarmupOptions.StableGPSampling = '200+10*nvars'; algoptions.gpQuadraticMeanBound = 1; algoptions.EmpiricalGPPrior = 0; algoptions.WarmupNoImproThreshold = 20 + 5*numel(probstruct.InitPoint); algoptions.TolStableExceptions = 2; algoptions.TolStableIters = 10; algoptions.WarmupCheckMax = true; algoptions.SGDStepSize = 0.005;        
        
    % New defaults
    case {100,'newdef'}; algoset = 'newdef'; algoptions.Plot = 0; algoptions.gpQuadraticMeanBound = 1; algoptions.EmpiricalGPPrior = 0; algoptions.WarmupNoImproThreshold = 20 + 5*numel(probstruct.InitPoint); algoptions.TolStableExceptions = 2; algoptions.TolStableIters = 10; algoptions.WarmupCheckMax = true; algoptions.SGDStepSize = 0.005;        
        
    % Variational active sampling
    case {1000,'vas'}; algoset = 'vas'; 
        
    otherwise
        error(['Unknown algorithm setting ''' algoset ''' for algorithm ''' algo '''.']);
end

% Increase base noise with noisy functions
if ~isempty(probstruct.Noise) || probstruct.IntrinsicNoisy
    algoptions.UncertaintyHandling = 'on';
    NoiseEstimate = probstruct.NoiseEstimate;
    if isempty(NoiseEstimate); NoiseEstimate = 1; end    
    algoptions.NoiseSize = NoiseEstimate(1);
else
    algoptions.UncertaintyHandling = 'off';
end

PLB = probstruct.PLB;
PUB = probstruct.PUB;
LB = probstruct.LB;
UB = probstruct.UB;
x0 = probstruct.InitPoint;
D = size(x0,2);

% Add log prior to function evaluation 
% (the current version of VBMC is agnostic of the prior)
probstruct.AddLogPrior = true;

algo_timer = tic;
[vp,elbo,elbo_sd,exitflag,~,output,stats] = ...
    vbmc(@(x) infbench_func(x,probstruct),x0,LB,UB,PLB,PUB,algoptions);
TotalTime = toc(algo_timer);

% Remove training data from GPs, too bulky (can be reconstructed)
for i = 1:numel(stats.gp)
    stats.gp(i).X = [];
    stats.gp(i).y = [];
end

if ~ControlRunFlag
    
    history = infbench_func(); % Retrieve history
    history.scratch.output = output;
    history.TotalTime = TotalTime;
    history.Output.stats = stats;
    
    % Store computation results (ignore points discarded after warmup)
    history.Output.X = output.X_orig(output.X_flag,:);
    history.Output.y = output.y_orig(output.X_flag);
    post.lnZ = elbo;
    post.lnZ_var = elbo_sd^2;
    [post.gsKL,post.Mean,post.Cov,post.Mode] = computeStats(vp,probstruct);

    % Return estimate, SD of the estimate, and gauss-sKL with true moments
    Nticks = numel(history.SaveTicks);
    for iIter = 1:Nticks
        idx = find(stats.N == history.SaveTicks(iIter),1);
        if isempty(idx); continue; end
        
        % Compute variational solution that VBMC would return at a given iteration
        [vp,elbo,elbo_sd,idx_best] = ...
            best_vbmc(stats,idx,algoptions.BestSafeSD,algoptions.BestFracBack,algoptions.RankCriterion);
                
        history.Output.N(iIter) = history.SaveTicks(iIter);
        history.Output.lnZs(iIter) = elbo;
        history.Output.lnZs_var(iIter) = elbo_sd^2;
        [gsKL,Mean,Cov,Mode] = computeStats(vp,probstruct);
        history.Output.Mean(iIter,:) = Mean;
        history.Output.Cov(iIter,:,:) = Cov;
        history.Output.gsKL(iIter) = gsKL;
        history.Output.Mode(iIter,:) = Mode;    
    end
else
    % Control condition -- run VBMC as normal but compute marginal likelihood
    % and posterior via other methods

    % Add WSABI algorithm to MATLAB path
    BaseFolder = fileparts(mfilename('fullpath'));
    AlgoFolder = 'wsabi';
    addpath(genpath([BaseFolder filesep() AlgoFolder]));
        
    history = infbench_func(); % Retrieve history
    
    % Start warmup
    
    % Store all points (no warmup pruning)
    X = output.X_orig(1:output.Xmax,:);
    y = output.y_orig(1:output.Xmax);
    Niter = find(size(X,1) == history.SaveTicks,1);
    N = history.SaveTicks(1:Niter);
    
    % Find when warm-up ends
    idx = find(stats.warmup == 0,1);
    if isempty(idx); endWarmupN = Inf; else; endWarmupN = stats.N(idx); end
    
    mu = zeros(1,Niter);
    ln_var = zeros(1,Niter);

    % Parameters for WSABI
    diam = probstruct.PUB - probstruct.PLB;
    kernelCov = diag(diam/10);     % Input length scales for GP likelihood model
    lambda = 1;                     % Ouput length scale for GP likelihood model
    alpha = 0.8;
    
    for iIter = 1:Niter        
        X_train = X(1:N(iIter),:);
        y_train = y(1:N(iIter));
        % Prune trials after warmup
        if N(iIter) >= endWarmupN
            idx_keep = output.X_flag(1:N(iIter));
            X_train = X_train(idx_keep,:);
            y_train = y_train(idx_keep);
        end
        X_iter{iIter} = X_train;
        y_iter{iIter} = y_train;
        lnp = infbench_lnprior(X_train,probstruct);
        y_train = y_train - lnp;  % Remove log prior for WSABI
        [mu(iIter),ln_var(iIter)] = ...
            wsabi_oneshot('L',probstruct.PriorMean,diag(probstruct.PriorVar),kernelCov,lambda,alpha,X_train,y_train);
    end
    vvar = max(real(exp(ln_var)),0);
    
    [history,post] = ...
        StoreAlgoResults(probstruct,[],Niter,X_iter{Niter},y_iter{Niter},mu,vvar,X_iter,y_iter,TotalTime);
    history.scratch.output = output;
    history.Output.stats = stats;    
    
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [gsKL,Mean,Cov,Mode] = computeStats(vp,probstruct)
%COMPUTE_STATS Compute additional statistics.
    
% Compute Gaussianized symmetric KL-divergence with ground truth
Ns_moments = 1e6;
xx = vbmc_rnd(vp,Ns_moments,1,1);
Mean = mean(xx,1);
Cov = cov(xx);
[kl1,kl2] = mvnkl(Mean,Cov,probstruct.Post.Mean,probstruct.Post.Cov);
gsKL = 0.5*(kl1 + kl2);

% Compute mode
Mode = vbmc_mode(vp,1);

end
