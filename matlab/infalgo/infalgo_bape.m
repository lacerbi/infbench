function [history,post,algoptions] = infalgo_bape(algo,algoset,probstruct)

% Add algorithm to MATLAB path
BaseFolder = fileparts(mfilename('fullpath'));
AlgoFolder = 'bapegp';
addpath(genpath([BaseFolder filesep() AlgoFolder]));

algoptions.Algorithm = 'bape';
algoptions.MaxFunEvals = probstruct.MaxFunEvals;
algoptions.MaxIter = Inf;
algoptions.Nsamples = 5e3;                 % Number of samples per iteration
algoptions.GPsamples = 80;
algoptions.AcqFun = @acqbapeEV;            % Exponentiated Variance acquisition function
algoptions.SamplingMethod = 'parallel';    % MCMC sampler for approximate posterior
algoptions.Plot = 0;                       % Make diagnostic plots at each iteration
algoptions.NcompMax = 30;                  % Maximum number of mixture components
algoptions.FracExpand = 0.1;               % Expand search box by this amount
algoptions.Meanfun = 'const';

if probstruct.Debug
    algoptions.TrueMean = probstruct.Post.Mean;
    algoptions.TrueCov = probstruct.Post.Cov;
end

% Use prior as proposal function
% algoptions.ProposalFcn = @(X_) exp(infbench_lnprior(X_,probstruct));

% Options from current problem
switch algoset
    case {0,'debug'}; algoset = 'debug'; algoptions.Plot = 1;
    case {1,'base'}; algoset = 'base';           % Use defaults
    case {2,'long'}; algoset = 'long'; algoptions.Nsamples = 2e4;
    case {3,'negquad'}; algoset = 'negquad'; algoptions.Meanfun = 'negquad';
    case {4,'nqdebug'}; algoset = 'nqdebug'; algoptions.Meanfun = 'negquad'; algoptions.Plot = 1;
        
    otherwise
        error(['Unknown algorithm setting ''' algoset ''' for algorithm ''' algo '''.']);
end

PLB = probstruct.PLB;
PUB = probstruct.PUB;
LB = probstruct.LB;
UB = probstruct.UB;
x0 = probstruct.InitPoint;
D = size(x0,2);

% BAPE tends to diverge on unconstrained problems, try with hard bounds
bounds_range = PUB - PLB;
idx = ~isfinite(bounds_range);
LB(idx) = PLB(idx) - 3*bounds_range(idx);
UB(idx) = PUB(idx) + 3*bounds_range(idx);

% Add log prior to function evaluation 
% (BAPE is agnostic of the prior)
probstruct.AddLogPrior = true;

algo_timer = tic;
[X,y,exitflag,output,vbmodel] = ...
    bapegp(@(x) infbench_func(x,probstruct),x0,LB,UB,PLB,PUB,algoptions);
TotalTime = toc(algo_timer);

Niter = numel(probstruct.SaveTicks);
%Nmax = numel(mu);
%idx = probstruct.SaveTicks(probstruct.SaveTicks <= Nmax);
%mu = mu(idx);

[history,post] = StoreAlgoResults(...
    probstruct,[],Niter,X,y,[],[],[],[],TotalTime);

stats = output.stats;

% Remove training data from GPs, too bulky (can be reconstructed)
for i = 1:numel(stats)
     stats(i).gp.X = [];
     stats(i).gp.y = [];
end

history.Output.stats = diagnostics;

% Remove vbmodel from stats
for i = 1:numel(history.Output.stats)
    history.Output.stats(i).vbmodel = [];
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [gsKL,Mean,Cov,Mode] = computeStats(stats,probstruct)
%COMPUTE_STATS Compute additional statistics.
    
vbmodel = stats.vbmodel;

% Compute Gaussianized symmetric KL-divergence with ground truth
Ns_moments = 1e6;   Nb = 20;
xx = [];
for i = 1:Nb; xx = [xx; vbgmmrnd(vbmodel,Ns_moments/Nb)']; end
Mean = mean(xx,1);
Cov = cov(xx);
[kl1,kl2] = mvnkl(Mean,Cov,probstruct.Post.Mean,probstruct.Post.Cov);
gsKL = 0.5*(kl1 + kl2);

% Compute mode
opts = optimoptions('fminunc','GradObj','off','Display','off');
Mode = fminunc(@(x) -vbgmmpdf(vbmodel,x')',stats.mode,opts);

end