function [nLL, fitParams, resp_model, resp_obs, p_true, p_estimate, lmL, vbmc_fit] = changeprob_maxlike(data,parameters)
%% CHANGEPROB_MAXLIKE Max log likelihood for specified model

% [Documentation to be written]

model = 1;
task = 2; % Covert task
MaxParams = 15;

% Parameters to be fit
if isempty(parameters); parameters = zeros(1,MaxParams); end
parameters(1) = 1;
parameters = [parameters,zeros(1,MaxParams-numel(parameters))];

paramNames = {'sigma_ellipse','sigma_criterion','lambda','gamma','alpha','w','Tmax', ...
    'pVec','beta','delta1','delta2','hRate','nu_p','delta_Tmax','delta_pVec'};
I_params = find(parameters);

% Lower and upper parameter bounds
paramBounds_def = [1,30; 1,30; 0,0.1; -Inf,Inf; 0,1; 0,1; 2,200; ...
    1e-6,0.5-1e-6; 0,10; 1.01,5; 1.01,14; 0,1; 0,5; 1,200; 1e-6,0.5-1e-6];    
paramBounds = paramBounds_def(I_params,:);

% Integer parameters (used only by VBMC)
%paramInteger = false(1,numel(paramNames));
%paramInteger([7,14]) = true;
%paramInteger = paramInteger(I_params);

%% Get session parameters

[NumTrials,sigma_ellipseData,mu,sigma_s,C,S,p_true,resp_obs,score] = ...
    changeprob_getSessionParameters(data,task);

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

%inputParams

%% Compute the min negative log likelihood using BADS

% DEFINE PARAMETER BOUNDS
LB = paramBounds(:,1);
UB = paramBounds(:,2);
PLB = LB + 0.1*(UB-LB);
PUB = LB + 0.9*(UB-LB);
  
% Define function to fit
nLL_fun = @(params2fit)changeprob_nll(params2fit(:)', I_params, inputParams, NumTrials, mu, ...
    sigma_s, C, S, p_true, resp_obs, task, score, model, X);
    
%% Recompute additional outputs from best fit params

%inputParams(I_params) = fitParams;
% sigma = sqrt(sigma_s^2 + inputParams(1)^2);
%[nLL,~,resp_model,p_estimate,~] = nLL_fun(fitParams);


end

