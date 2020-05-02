function changeprob_runfit()
%RUNFIT Runs model comparison for changing probability experiment
%   Input:
    % jobNumber: unique model, task, data type, and subject identifier
    % fitType: fits data via log marginal
    % likelihood ('logmarglike') or max likelihood ('maxlike')
    % fixNoise: fixes the noise parameter
    % gridSize: determines the size of the grid for computing the log
    % marginal likelihood ('logmarglike' only)
    % Nopts: # optimization runs for maximum-likelihood fits ('maxlike' only)
    
% Output:
   % logmarglikelihood: log marginal likelihood
   % modelPost: the model posterior approximated over a grid and used to
   % compute the log marginal likelihood
   % nLL: negative log likelihood for the best fit parameters
   % rmse: root mean squared error between the model and data
   % fitParams: vector of the best fit parameters
   % resp_model: model prediction on each trial
   % resp_obs: observers response on each trial
   % p_true: true probability of category A on each trial
   % p_estimate: model estimate of the probability of category A on each
   % trial
   % post: ideal Bayesian model posterior

subID = {'CWG','EGC','EHN','ERK','GK','HHL','JKT','JYZ','RND','SML','SQC'};

% Models to fit (only Bayesian models available)
models = {[],'idealBayesian',[],[],[],[],[], ...
    'subBayesian_rlprior','subBayesian_conservative','subBayesian_pVec',...
    'subBayesian_betahyp','subBayesian_3param',[],[],'subBayesian_flex',...
    [],[],[]};

subIndex = 1;   % Pick subject
runSubject = subID{subIndex};
task = 2;       % Covert task
modelIndex = 2; % Pick model
runModel = models{modelIndex};

% Add project directory and subdirs to path
matlabdir = fileparts(which('changeprob_logmarglike'));
basedir = matlabdir(1:find(matlabdir == filesep(), 1, 'last')-1);
addpath(genpath(basedir));

% Load data
load(['ChangingProbabilities_', runSubject]);

parameters = [];
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
end
runModel = 'idealBayesian';

nLL = changeprob_maxlike(runModel,data,task,parameters,[],1);

end

