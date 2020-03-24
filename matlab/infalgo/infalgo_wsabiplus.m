function [history,post,algoptions] = infalgo_wsabiplus(algo,algoset,probstruct)

% Add algorithm to MATLAB path
BaseFolder = fileparts(mfilename('fullpath'));
AlgoFolder = 'wsabiplus';
addpath(genpath([BaseFolder filesep() AlgoFolder]));

HypVar = [100^2,4^2*ones(1,probstruct.D)];

algoptions.Method = 'L';    % Default is WSABI-L
algoptions.Alpha = 0.8;     % Fractional offset, as in paper.
algoptions.Nsearch = 0;     % Extra search points
algoptions.HypVar  = 1;     % Variance of GP hyperparamters (WSABI default)
algoptions.Debug = false;
algoptions.MaxFunEvals = probstruct.MaxFunEvals;
algoptions.HypVar = HypVar;

% Options from current problem
switch algoset
    case {0,'debug'}; algoset = 'debug'; algoptions.Debug = 1;
    case {1,'base'}; algoset = 'base';           % Use defaults
    case {2,'mm'}; algoset = 'mm'; algoptions.Method = 'M';
        
    otherwise
        error(['Unknown algorithm setting ''' algoset ''' for algorithm ''' algo '''.']);
end

% Increase base noise with noisy functions
if ~isempty(probstruct.Noise) || probstruct.IntrinsicNoisy
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

range = [PLB - 3*diam; PUB + 3*diam];

% Do not add log prior to function evaluation, already passed to WSABI 
probstruct.AddLogPrior = false;

if algoptions.Debug; algoptions.Display = 2; else; algoptions.Display = 1; end

algo_timer = tic;
[mu, ln_var, tt, X, y, hyp] = ...
    wsabiplus(@(x) infbench_func(x,probstruct),...
    probstruct.PriorMean,diag(probstruct.PriorVar),...
    range(1,:),range(2,:),x0,...
    algoptions);
TotalTime = toc(algo_timer);

vvar = max(real(exp(ln_var)),0);

Niter = numel(probstruct.SaveTicks);
Nmax = numel(mu);
idx = probstruct.SaveTicks(probstruct.SaveTicks <= Nmax);
mu = mu(idx);
vvar = vvar(idx);

[history,post] = StoreAlgoResults(...
    probstruct,[],Niter,X(end:-1:1,:),y(end:-1:1),mu,vvar,[],[],TotalTime);

history.Output.stats.tt = tt;
history.Output.stats.hyp = hyp;


end