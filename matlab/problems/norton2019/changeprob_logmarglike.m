function [lmL, modelPost, nLL, rmse, fitParams, resp_model,...
    resp_obs, p_true, p_estimate, post, v_estimate] = changeprob_logmarglike(model, data, task, parameters, gridSize, paramBounds, fixNoise, simulatedData)

%% CHANGEPROB_LOGMARGLIKE Log marginal likelihood for specified model

% Here we compute the log marginal likelihood of the specfied model. 
% To calculate the marginal likelihood, we multiply the probability of 
% the data given the model by the prior of the model, integrate over 
% the model parameters, and then take the logarithm.

% INPUT:
    % model: a character string indicating the model you want to fit, 
    % which is limited to the following:
        % 'idealBayesian'
        % 'fixed'
        % 'exponential'
        % 'RL_probability'
        % 'RL_criterion'
        % 'gold'
        % 'behrens'
    % data: data struct from changing probability experiment
    % task: lets you choose which task to fit
        % 1 - overt-criterion task
        % 2 - covert-criterion task
        % 3 - mixed design (overt-criterion task on every 5th trial)
    % parameters: Vector indicating the parameters to-be-fit (1 - fit, 0 - not fit)
        % [sigma_ellipse, sigma_criterion, lapse, gamma, alpha, w, Tmax, pVec, betahyp, delta1, delta2, hRate, nu_p]
    % gridSize: size of parameter grid (e.g., n or [n x m x p])
    % paramBounds: lower and upper parameter bounds (numel(gridSize) x 2)
    % fixNoise: fix noise from measurement session (leave empty for default,
    %           which is 0 for covert, 1 for overt, and 1 for mixed)

% OUTPUT:
    % lmL: a measure of model fit
    % modelPost: proportional to the p(parameters | data, model)
    % nLL: negative log likelihood that correponds to the best fitting
    % parameters
    % rmse: root mean squared error between the model probability and
    % the true probability (Note: does not apply to the RL_criterion model)
    % fitParams: best fitting model parameters computed by taking the
    % MAP of modelPost
    % resp_model: predicted response/criterion
    % resp_obs: observer's response/criterion
    % p_true: true probability values
    % p_estimate: model's estimate of probability

% Authors:  Elyse Norton, Luigi Acerbi
% Email:    {elyse.norton,luigi.acerbi}@gmail.com
% Date:     4/21/2017
    
if nargin < 2; data = []; end
if nargin < 3; task = []; end
if nargin < 4; parameters = []; end
if nargin < 5; gridSize = []; end
if nargin < 6; paramBounds = []; end
if nargin < 7; fixNoise = []; end
if nargin < 8; simulatedData = []; end

%tic
% Model to be fit
if nargin < 1; error('Please indicate the model you want to fit.'); end
potentialModels = {'idealBayesian', 'fixed', 'exponential', 'RL_probability', ...
    'RL_criterion', 'gold', 'behrens', 'behrens_jump'};
model = find(strcmp(model, potentialModels)); % recode model to numeric value

% Data struct or random seed for fake data generation
if isempty(data); data = 0; end
if isnumeric(data); rng(data); data = []; end

% Task (1 overt, 2 covert)
if isempty(task); task = 1; end % overt-criterion task is the default
if task ~= 1 && task ~= 2 && task ~= 3; error('TASK can only be 1 (overt), 2 (covert), or 3 (mixed).'); end
if task == 1; taskName = 'overt'; elseif task == 2; taskName = 'covert'; else taskName = 'mixed'; end

NumSamples = 5000;
MaxParams = 13;

% Parameters to be fit
if isempty(parameters)
    parameters = zeros(1,MaxParams);    
    switch task
        case 1
            if or(model < 3, model >= 7)
                parameters(2) = 1;
            elseif and(model > 2, model < 6)
                parameters([2,5]) = 1;
            else
                parameters([2,10,11]) = 1;
            end
        case 2
            if or(model < 3, model >= 7)
                parameters(1) = 1;
            elseif and(model > 2, model < 6)
                parameters([1,5]) = 1;
            else
                parameters([1,10,11]) = 1;
            end
        case 3
            if or(model < 3, model >= 7)
                parameters(2) = 1;
            elseif and(model > 2, model < 6)
                parameters([2,5]) = 1;
            else
                parameters([2,10,11]) = 1;
            end
    end
elseif sum(parameters) == 0
    error('You must specify at least one parameter to fit.')
end

if isempty(fixNoise)
    switch task
        case 1; parameters(1) = 0;  % Overt task: by default noise is fixed
        case 2; parameters(1) = 1;  % Covert task: by default noise is free
        case 3; parameters(1) = 0;  % Mixed task: by default noise is fixed
    end
else
    if fixNoise; parameters(1) = 0; else parameters(1) = 1; end
end

paramNames = {'sigma_ellipse', 'sigma_criterion', 'lambda', 'gamma', 'alpha', 'w', 'Tmax', 'pVec', 'beta', 'delta1', 'delta2', 'hRate', 'nu_p'};
NumParams = sum(parameters);
if NumParams > 0; fitParamNames = paramNames{logical(parameters)}; end 
I_params = find(parameters ~= 0);
% Sampling rate of parameter grid
if isempty(gridSize)
    gridSize = 100*ones(1, NumParams); % Default grid is 100 x 100
elseif isscalar(gridSize)
    gridSize = gridSize*ones(1, NumParams);
elseif numel(gridSize) ~= NumParams
    error('Matrix dimensions do not agree. numel(gridSize) must equal sum(parameters).');
end
fprintf('Grid for computation of the marginal likelihood has %d nodes.\n', prod(gridSize));

% Lower and upper parameter bounds
if isempty(paramBounds)    
    paramBounds_def = [1,30; 1,30; 0,0.1; -Inf,Inf; 0,1; 0,1; 2,200; 0,.5; 0,10; 1.01,5; 1.01,14; 0,1; 0,5];    
    paramBounds = paramBounds_def(I_params,:);
end

if NumParams ~= size(paramBounds,1)
    error('Please specify parameter bounds for each parameter to be fit.');
end

%% Print session
if ~isempty(fixNoise) || and(task ~= 2, parameters(1) == 0); noiseString = 'FIXED'; else noiseString = 'FREE'; end
fprintf('Fitting model %s, %s-criterion task; sensory noise is %s, %d free parameters.\n', potentialModels{model}, taskName, noiseString, NumParams);

%% Get session parameters
if isempty(simulatedData)
    [NumTrials, sigma_ellipseData, mu, sigma_s, C, S, p_true, resp_obs, score] = changeprob_getSessionParameters(data, task);
else
    NumTrials = data.NumTrials;
    sigma_ellipseData = data.sigmaEllipse;
    mu = data.mu;
    sigma_s = data.sigma_s;
    C = data.Category;
    S = data.Stimulus;
    p_true = data.pA;
    resp_obs = data.response;
    score = data.score;
    if size(resp_obs,1) == 1
        resp_obs = resp_obs';
        score = score';
    end
end
if and(task == 1, model == 5) || and(task == 3, model == 5)
    X = bsxfun(@plus, S, sigma_ellipseData*randn(numel(S), NumSamples));
else
    X = [];
end

%% Observer model parameters
if NumParams > 0
    params2fit = zeros(NumParams, gridSize(1));
    for iParam = 1:NumParams
        switch I_params(iParam)
            case {1, 2}
                params2fit(iParam,:) = log(linspace(paramBounds(iParam,1), paramBounds(iParam,2), gridSize(iParam))); % sigma_noise
                % Edit: should have been the following to place a uniform prior in log space:
                % params2fit(iParam,:) = linspace(log(paramBounds(iParam,1)), log(paramBounds(iParam,2)), gridSize(iParam))); % sigma_noise
            case 3
                params2fit(iParam,:) = linspace(paramBounds(iParam,1), paramBounds(iParam,2), gridSize(iParam)); % lambda
            case 4
                params2fit(iParam,:) = linspace(paramBounds(iParam,1), paramBounds(iParam,2), gridSize(iParam)); % gamma
            case 5
                params2fit(iParam,:) = linspace(paramBounds(iParam,1), paramBounds(iParam,2), gridSize(iParam)); % alpha
            case 6
                params2fit(iParam,:) = linspace(paramBounds(iParam,1), paramBounds(iParam,2), gridSize(iParam)); % w
            case 7
                params2fit(iParam,:) = round(linspace(paramBounds(iParam,1), paramBounds(iParam,2), gridSize(iParam))); % Tmax (discrete)
            case 8
                params2fit(iParam,:) = linspace(paramBounds(iParam,1), paramBounds(iParam,2), gridSize(iParam)); % Range for minimum probability (pVec(1))
            case 9
                params2fit(iParam,:) = linspace(paramBounds(iParam,1), paramBounds(iParam,2), gridSize(iParam)); % beta (hyperprior)
            case {10, 11}
                params2fit(iParam,:) = linspace(paramBounds(iParam,1), paramBounds(iParam,2), gridSize(iParam)); % distance between nodes (log of the true deltas)
            case 12
                params2fit(iParam,:) = linspace(paramBounds(iParam,1), paramBounds(iParam,2), gridSize(iParam)); % hazard rate
            case 13
                params2fit(iParam,:) = linspace(paramBounds(iParam,1), paramBounds(iParam,2), gridSize(iParam)); % nu_p
        end
    end

    % Create parameter list for fitting
    for iParam = 1:NumParams; params2fit_cell{iParam} = params2fit(iParam,:); end
    params2fit_list = combvec(params2fit_cell{:})';    
    exp_idx = I_params == 1 | I_params == 2;
    params2fit_list(:,exp_idx) = exp(params2fit_list(:,exp_idx));   % Check this
    clear params2fit_cell;
end

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
notFit_def = [sigma_ellipseData, sigmacriterion_def, lapse_def, gamma_def, alpha_def, w_def, Tmax_def, pVec_def, beta_def, delta_def, delta_def, hRate_def, nu_p_def];
inputParams(I_notFit) = notFit_def(I_notFit);

inputParams

%% Separate case if there are no free parameters
if NumParams == 0
    sigma = sqrt(sigma_s^2 + inputParams(1)^2);
    [nLL,rmse,resp_model,p_estimate,post] = ...
        changeprob_nll(inputParams, NumTrials, mu, sigma, C, S, p_true, resp_obs, task, score, model, X);
    lmL = -nLL;     % Log marginal likelihood is simply log likelihood
    modelPost = [];
    fitParams = [];    
    return;
end

%% Choose priors - start with uniformative priors for all parameters
    % These will be uniform for all priors since we took the log of the
    % sigma parameters (edit: actually in the end we are using uniform
    % priors for all parameters)
    
% Uniform prior for all parameters
prior = 1/prod(params2fit(:,end) - params2fit(:,1));
logprior = log(prior);

% Edit: Note that since we are using uniform priors also for SIGMA against
% our original intentions (but we took the log above), we need then to 
% apply a correction to the log marginal likelihood which turns out to be 
% log(log(30) - log(1)) - (log(30)-log(1)) = -2.1771
% for each SIGMA parameter activated in the model.

%% Compute the negative log likelihood for all parameter sets

nLL_mat = zeros([gridSize,1]);
maxii = prod(gridSize);
timestart = tic;
for ii = 1:maxii
    % Timing test
    if ii == 5; fprintf('Timing test: %.4g s.\n', toc(timestart)); end
    if rem(ii,500) == 0; fprintf('%.1f%%..', 100*ii/maxii); end
    % Parameters in params2fit_list are already transformed
%     inputParams(I_params) = params2fit_list(ii,:);
%     sigma = sqrt(sigma_s^2 + inputParams(1)^2);
    nLL_mat(ii) = changeprob_nll(params2fit_list(ii,:), I_params, inputParams, NumTrials, mu, sigma_s, C, S, p_true, resp_obs, task, score, model, X);
end
fprintf('\n');
clear params2fit_list;  % Free up some memory

%% Compute the marginal likelihood

% Add the log likelihood to the log prior, subtract the max, and exponentiate
logUpost = bsxfun(@plus, -nLL_mat, logprior);
maxlogUpost = max(logUpost(:));
modelUpost = exp(logUpost - maxlogUpost);   % Unnormalized posterior       

% Trapezoidal integration across all parameters
temp = modelUpost;
for iParam = 1:NumParams; temp = qtrapz(temp); end
dV = prod(params2fit(:,2) - params2fit(:,1));   % volume element
Z = temp*dV;    % Normalization constant

% Normalized posterior
modelPost = modelUpost / Z;

% Log marginal likelihood (add back max logUpost from before)
lmL = log(Z) + maxlogUpost;

% Find the MAP
[~,linidx] = max(logUpost(:));
[idx(1),idx(2),idx(3),idx(4),idx(5)] = ind2sub(size(logUpost),linidx);
for iParam = 1:NumParams; bestFit_param(iParam) = params2fit(iParam,idx(iParam)); end

% Transform parameters
fitParams = bestFit_param;
fitParams(exp_idx) = exp(fitParams(exp_idx));

% Recompute additional outputs from best fit params
% inputParams(I_params) = fitParams;
% sigma = sqrt(sigma_s^2 + inputParams(1)^2);
[nLL,rmse,resp_model,p_estimate,post, v_estimate] = ...
    changeprob_nll(fitParams, I_params, inputParams, NumTrials, mu, sigma_s, C, S, p_true, resp_obs, task, score, model, X);

%toc
end