function [history,post,algoptions] = infalgo_wsabiplus(algo,algoset,probstruct)

% Add algorithm to MATLAB path
BaseFolder = fileparts(mfilename('fullpath'));
AlgoFolder = 'wsabiplus';
addpath(genpath([BaseFolder filesep() AlgoFolder]));

algoptions.Method = 'L';    % Default is WSABI-L
algoptions.AlphaFactor = 0.8;     % Fractional offset, as in paper.
algoptions.Nsearch = 2^10;  % Starting search points
algoptions.GPInputHypVar  = log(20)^2;   % Variance of GP input hyperparameters (default)
algoptions.GPOutputHypVar  = log(10)^2;   % Variance of GP output hyperparameter (default)
algoptions.Debug = false;
algoptions.MaxFunEvals = probstruct.MaxFunEvals;
algoptions.IgnoreNoise = false;
algoptions.LCBFactor = 2; % Lower Confidence Bound parameter
algoptions.PreciseSearch = true; % Assume zero noise during active search
algoptions.BoundedTransform = 'logit';
algoptions.NMCMCsamples = 2e4;
algoptions.SavePosterior = [10:10:probstruct.MaxFunEvals,Inf];
%algoptions.SavePosterior = [probstruct.SaveTicks,Inf];

% Options from current problem
switch algoset
    case {0,'debug'}; algoset = 'debug'; algoptions.Debug = 1;
    case {1,'base'}; algoset = 'base';           % Use defaults
    case {2,'mm'}; algoset = 'mm'; algoptions.Method = 'M';
    case {3,'ldet'}; algoset = 'ldet'; algoptions.Method = 'L'; algoptions.IgnoreNoise = true;        
    case {4,'probit'}; algoset = 'probit'; algoptions.BoundedTransform = 'probit';
        
    otherwise
        error(['Unknown algorithm setting ''' algoset ''' for algorithm ''' algo '''.']);
end

% Increase base noise with noisy functions
if (~isempty(probstruct.Noise) || probstruct.IntrinsicNoisy) && ~algoptions.IgnoreNoise
    algoptions.SpecifyTargetNoise = true;
end
%     algoptions.UncertaintyHandling = 'on';
%     NoiseEstimate = probstruct.NoiseEstimate;
%     if isempty(NoiseEstimate); NoiseEstimate = 1; end    
%     algoptions.NoiseSize = NoiseEstimate(1);
% else
%     algoptions.UncertaintyHandling = 'off';
% end

PLB = probstruct.PLB;
PUB = probstruct.PUB;
LB = probstruct.LB;
UB = probstruct.UB;
x0 = probstruct.InitPoint;
D = size(x0,2);

diam = probstruct.PUB - probstruct.PLB;

algoptions.GPInputLengths = diam/sqrt(10);   % Input length scales for GP likelihood model
algoptions.GPOutputLength = 1;                     % Ouput length scale for GP likelihood model

switch lower(probstruct.PriorType)
    case 'gaussian'
        algoptions.PriorMean = probstruct.PriorMean;
        algoptions.PriorCov = probstruct.PriorVar;
        % Do not add log prior to function evaluation, already passed to WSABI 
        probstruct.AddLogPrior = false;
    
    case 'uniform'
        algoptions.PriorMean = [];
        algoptions.PriorCov = [];
        probstruct.AddLogPrior = true;    
end

if algoptions.Debug; algoptions.Display = 2; else; algoptions.Display = 1; end

% Run WSABI algorithm
algo_timer = tic;
[mu,ln_var,output,mcmc_samples] = ...
    wsabiplus(@(x) infbench_func(x,probstruct),...
    PLB,PUB,LB,UB,x0,algoptions);
TotalTime = toc(algo_timer);
vvar = max(real(exp(ln_var)),0);
tt = output.clktime;

% Postprocessing
history = infbench_func(); % Retrieve history
history.scratch.output = output;
history.TotalTime = TotalTime;

%s2 = output.fval_sd.^2;

% Store computation results (ignore points discarded after warmup)
history.Output.X = output.X_transformed;
history.Output.y = output.fval_transformed;
post.lnZ = mu(end);
post.lnZ_var = vvar(end);
fprintf('Calculating WSABI output at iteration...\n');
fprintf('%d..',0);
[post.gsKL,post.Mean,post.Cov,post.Mode,post.MTV,post.samples,post.Test] = ...
    computeStats(mcmc_samples{end},probstruct);

% Return estimate, SD of the estimate, and gauss-sKL with true moments
Nticks = numel(probstruct.SaveTicks);
for iIter = 1:Nticks
    fprintf('%d..',iIter);

    idx = find(output.mcmc_iters == history.SaveTicks(iIter),1);
    if isempty(idx); continue; end

    history.ElapsedTime(iIter) = tt(history.SaveTicks(iIter));

    history.Output.N(iIter) = history.SaveTicks(iIter);
    history.Output.lnZs(iIter) = mu(history.SaveTicks(iIter));
    history.Output.lnZs_var(iIter) = vvar(history.SaveTicks(iIter));
    [gsKL,Mean,Cov,Mode,MTV,~,Test] = computeStats(mcmc_samples{idx},probstruct);
    history.Output.Mean(iIter,:) = Mean;
    history.Output.Cov(iIter,:,:) = Cov;
    history.Output.gsKL(iIter) = gsKL;
    history.Output.Mode(iIter,:) = Mode;
    history.Output.MTV(iIter,:) = MTV;
    history.Output.Test{iIter} = Test;
end
fprintf('\n');

%X_transformed = output.X_transformed;
%y_transformed = output.fval_transformed;
%s2 = output.fval_sd.^2;
%hyp = output.gp_hyp;
%inverse_transform = @(x) warpvars_wsabi(x,'inv',output.trinfo);

%Niter = numel(probstruct.SaveTicks);
%Nmax = numel(mu);
%idx = probstruct.SaveTicks(probstruct.SaveTicks <= Nmax);
%mu = mu(idx);
%vvar = vvar(idx);

%[history,post] = StoreAlgoResults(...
%    probstruct,[],Niter,X_transformed,y_transformed,mu,vvar,[],[],TotalTime,[],s2,[],inverse_transform);

history.Output.stats.tt = tt;
history.Output.stats.hyp = output.gp_hyp;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [gsKL,Mean,Cov,Mode,MTV,xx,Test] = computeStats(xx,probstruct)
%COMPUTE_STATS Compute additional statistics.
    
D = size(xx,2);

% Compute Gaussianized symmetric KL-divergence with ground truth
Mean = mean(xx,1);
Cov = cov(xx);
[kl1,kl2] = mvnkl(Mean,Cov,probstruct.Post.Mean,probstruct.Post.Cov);
gsKL = 0.5*(kl1 + kl2);

% Compute mode
Mode = NaN(1,D);

% Compute marginal total variation
try
    MTV = ComputeMarginalTotalVariation(xx,probstruct);
catch
    MTV = NaN(1,D);
end

% Compute test statistics (problem-specific)
if probstruct.ComputeTestStatistics
    fun = str2func(['infbench_' probstruct.Prob]);
    Test = fun('test',probstruct.ProbInfo,xx);
    Test
else
    Test = [];
end

gsKL
MTV

end