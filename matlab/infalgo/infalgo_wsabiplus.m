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

algo_timer = tic;
[mu,ln_var,output] = ...
    wsabiplus(@(x) infbench_func(x,probstruct),...
    PLB,PUB,LB,UB,x0,algoptions);
TotalTime = toc(algo_timer);

tt = output.clktime;
X = output.X;
y = output.fval;
s2 = output.fval_sd.^2;
hyp = output.gp_hyp;

vvar = max(real(exp(ln_var)),0);

Niter = numel(probstruct.SaveTicks);
Nmax = numel(mu);
idx = probstruct.SaveTicks(probstruct.SaveTicks <= Nmax);
mu = mu(idx);
vvar = vvar(idx);

[history,post] = StoreAlgoResults(...
    probstruct,[],Niter,X,y,mu,vvar,[],[],TotalTime,[],s2);

history.Output.stats.tt = tt;
history.Output.stats.hyp = hyp;


end