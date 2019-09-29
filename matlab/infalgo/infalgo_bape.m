function [history,post,algoptions] = infalgo_bape(algo,algoset,probstruct)

% Add algorithm to MATLAB path
BaseFolder = fileparts(mfilename('fullpath'));
AlgoFolder = 'bapegp';
addpath(genpath([BaseFolder filesep() AlgoFolder]));

algoptions.Algorithm = 'bape';
algoptions.MaxFunEvals = probstruct.MaxFunEvals;
algoptions.MaxIter = Inf;
algoptions.NsMax_gp = 80;
algoptions.AcqFun = @acqbapeEV;            % Exponentiated Variance acquisition function
algoptions.SamplingMethod = 'parallel';    % MCMC sampler for approximate posterior
algoptions.Plot = 0;                       % Make diagnostic plots at each iteration
algoptions.NcompMax = 30;                  % Maximum number of mixture components
algoptions.FracExpand = 0.1;               % Expand search box by this amount

% These modifications are necessary for BAPE to converge
algoptions.Meanfun = 'negquad';             % Negative quadratic mean function
algoptions.TolGPVar = 1e-4;                 % GP regularization

if probstruct.Debug
    algoptions.TrueMean = probstruct.Post.Mean;
    algoptions.TrueCov = probstruct.Post.Cov;
end

% Use prior as proposal function
% algoptions.ProposalFcn = @(X_) exp(infbench_lnprior(X_,probstruct));

% Options from current problem
switch algoset
    case {0,'debug'}; algoset = 'debug'; algoptions.Plot = 1;
    case {1,'base'}; algoset = 'base';           % Use good defaults
    case {2,'original'}; algoset = 'original'; algoptions.Meanfun = 'const'; algoptions.TolGPVar = 0; algoptions.NsMax_gp = 0;        
    case {3,'step1'}; algoset = 'step1'; algoptions.Nstep = 1;
        
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
LB(idx) = PLB(idx) - 5*bounds_range(idx);
UB(idx) = PUB(idx) + 5*bounds_range(idx);

% Add log prior to function evaluation 
% (BAPE is agnostic of the prior)
probstruct.AddLogPrior = true;

algo_timer = tic;
[X,y,exitflag,output,vbmodel] = ...
    bapegp(@(x) infbench_func(x,probstruct),x0,LB,UB,PLB,PUB,algoptions);
TotalTime = toc(algo_timer);

Niter = numel(probstruct.SaveTicks);

N = [output.stats.N];
Nmax = max(N);
idx = probstruct.SaveTicks(probstruct.SaveTicks <= Nmax);
gpiter = [];
for ii = 1:numel(idx)
    pos = find(N == idx(ii),1);    
    if ~isempty(pos)
        gpiter{end+1} = output.stats(pos).gp;
    end
end

[history,post] = StoreAlgoResults(...
    probstruct,[],Niter,X,y,[],[],[],[],TotalTime,gpiter);

stats = output.stats;

% Remove training data from GPs, too bulky (can be reconstructed)
for i = 1:numel(stats)
     stats(i).gp.X = [];
     stats(i).gp.y = [];
end

history.Output.stats = stats;

% Remove vbmodel from stats
for i = 1:numel(history.Output.stats)
    history.Output.stats(i).vbmodel = [];
end

end