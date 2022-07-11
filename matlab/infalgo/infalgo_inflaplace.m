function [history,post,algoptions] = infalgo_inflaplace(algo,algoset,probstruct)

% Add algorithm to MATLAB path
BaseFolder = fileparts(mfilename('fullpath'));
AlgoFolder = 'inflaplace';
addpath(genpath([BaseFolder filesep() AlgoFolder]));

algoptions.MaxFunEvals = probstruct.MaxFunEvals;

% Options from current problem
switch algoset
    case {0,'debug'}; algoset = 'debug'; algoptions.Debug = 1; algoptions.Plot = 'scatter';
    case {1,'base'}; algoset = 'base';           % Use defaults    
    case {2,'probit'}; algoset = 'probit'; algoptions.BoundedTransform = 'probit';
    otherwise
        error(['Unknown algorithm setting ''' algoset ''' for algorithm ''' algo '''.']);
end

PLB = probstruct.PLB;
PUB = probstruct.PUB;
LB = probstruct.LB;
UB = probstruct.UB;
x0 = probstruct.InitPoint;
D = size(x0,2);

% Add log prior to function evaluation 
probstruct.AddLogPrior = true;

% Compute Laplace approximation (assume mode is already given)
%--------------------------------------------------------------------------
algo_timer = tic;

% SaveTicks = probstruct.SaveTicks(probstruct.SaveTicks <= algoptions.MaxFunEvals);
fun = @(x_) infbench_func(x_,probstruct);
x0 = probstruct.Post.Mode;

[mu,Sigma,lnZ,output] = inflaplace(fun,x0,PLB,PUB,LB,UB,algoptions);

TotalTime = toc(algo_timer);
%--------------------------------------------------------------------------

post = [];

history = infbench_func(); % Retrieve history
% history.scratch.output = output;
history.TotalTime = TotalTime;
history.Output.X = [];
history.Output.y = [];

post.lnZ = lnZ;
post.lnZ_var = 0;
post.trinfo = output.trinfo;

[post.gsKL,post.Mean,post.Cov,post.Mode,post.MTV,post.samples,post.Test] = ...
        computeStats(mu,Sigma,probstruct,post.trinfo);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [gsKL,Mean,Cov,Mode,MTV,xx,Test] = computeStats(mu,Sigma,probstruct,trinfo)
%COMPUTE_STATS Compute additional statistics.
    
% Compute Gaussianized symmetric KL-divergence with ground truth
Ns_moments = 1e6;
xx = mvnrnd(mu,Sigma,Ns_moments);
xx = warpvars_laplace(xx,'i',trinfo);

Mean = mean(xx,1);
Cov = cov(xx);
[kl1,kl2] = mvnkl(Mean,Cov,probstruct.Post.Mean,probstruct.Post.Cov);
gsKL = 0.5*(kl1 + kl2);

% Compute mode
Mode = warpvars_laplace(mu,'i',trinfo);

% Compute marginal total variation
try
    MTV = ComputeMarginalTotalVariation(xx,probstruct);
catch
    MTV = NaN(1,numel(mu));
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