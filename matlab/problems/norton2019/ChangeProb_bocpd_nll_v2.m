function [nLL, rmse, p_estimate, resp_model, post] = ChangeProb_bocpd_nll_v2(parameters, NumTrials, mu, sigma, C, S, p_true, resp_obs, score, task, prior_rl, p_vec, beta_hyp)
%CHANGEPROB_BOCPD_NLL Bayesian online changepoint detection observer.
% (Documentation to be written.)
%
% Author:   Luigi Acerbi
% Email:    luigi.acerbi@gmail.com
% Date:     Oct/1/2016

% Modified by Elyse Norton on Oct/27/2016

% Parameter vector:
% #1 is SIGMA_ELLIPSE, #2 is SIGMA_CRITERION, #3 is LAPSE, #4 is GAMMA,
% #5 is ALPHA, and #6 is W
if nargin < 1
    parameters = [];
    prior_rl = [];
    p_vec = [];
    [NumTrials, sigma_ellipse, mu, sigma, C, S, p_true, resp_obs] = changeprob_getSessionParameters();
    task = 1; 
end

% Data struct or random seed for fake data generation
if nargin < 2; error('You must specify the session parameters.'); end

% Task (1 overt, 2 covert, 3 mixed)
if nargin < 10 || isempty(task); task = 1; end
if task ~= 1 && task ~= 2 && task ~=3; error('TASK can only be 1 (overt-criterion), 2 (covert-criterion), or 3 (mixed).'); end

if nargin < 11 || isempty(prior_rl)
    prior_rl = [80,120];    % Default runlengths
end

if nargin < 12 || isempty(p_vec)
    p_vec = linspace(0.2,0.8,5); % Default states
end

if nargin < 13 || isempty(beta_hyp)
    beta_hyp = [0,0];   % Default beta hyperprior after changepoint (maximum-likelihood)
end
if isscalar(beta_hyp); beta_hyp = beta_hyp*[1,1]; end

%% Experiment constants

% Run length ~ Uniform[prior_rl(1),prior_rl(2)]
% if isfield(options,'prior_rl') && ~isempty(options.prior_rl)
%     prior_rl = options.prior_rl;
% else
%end
runlength_min = prior_rl(1);
runlength_max = prior_rl(2);

% Probability states
% if isfield(options,'p_vec') && ~isempty(options.p_vec)
%     p_vec = options.p_vec;
% else
%     p_vec = linspace(0.2,0.8,5);    % Default states
%end
Nprobs = numel(p_vec);              % # states

p_bias = .5;

if task == 3
    idx_Overt = 5:5:NumTrials;
    idx_Covert = 1:NumTrials;
    idx_Covert(idx_Overt) = [];
end

%% Observer model parameters

% Default values
lambda = 0;
gamma = Inf;
w = 1;

switch numel(parameters)
    case 0  % Empty parameter vector
        sigma_criterion = sigma_ellipse;
    case 1
        sigma_ellipse = parameters(1);
        sigma_criterion = sigma_ellipse;
    case 2
        sigma_ellipse = parameters(1);
        sigma_criterion = parameters(2);        
    case 3
        sigma_ellipse = parameters(1);
        sigma_criterion = parameters(2);
        lambda = parameters(3);
    case 4
        sigma_ellipse = parameters(1);
        sigma_criterion = parameters(2);
        lambda = parameters(3);
        gamma = parameters(4);
    case 5
        sigma_ellipse = parameters(1);
        sigma_criterion = parameters(2);
        lambda = parameters(3);
        gamma = parameters(4);
    case 6
        sigma_ellipse = parameters(1);
        sigma_criterion = parameters(2);
        lambda = parameters(3);
        gamma = parameters(4);
        w = parameters(6);
    otherwise
        error('PARAMETERS should be a vector with up to six parameters.');
end

%% Initialize inference

% Hazard function (defined only where nonzero)
L = (ceil(runlength_max)-floor(runlength_min))+1;  % Number of elements
rlprior = ones(L,1)/L;      % Constant changepoint prior

% Trick to allow non-integer run length min and max
frac_min = 1 - (runlength_min - floor(runlength_min));
rlprior(1) = rlprior(1)*frac_min;
frac_max = 1 - (ceil(runlength_max) - runlength_max);
rlprior(end) = rlprior(end)*frac_max;
H = rlprior(:)./flipud(cumsum(flipud(rlprior(:))));

% Posterior over run lengths (from 0 to PRIOR_RL(2)-1)
post = zeros(ceil(runlength_max),Nprobs);
post(1,:) = 1;  % Change in the first trial
post = post ./ sum(post(:));

% Table of binomial count/probability
Psi = zeros(size(post,1),1,Nprobs); Psi(1,1,:) = 1/Nprobs;

% Transition matrix
Tmat(1,:,:) = (ones(Nprobs,Nprobs) - eye(Nprobs)) / (Nprobs-1);
% Tmat(1,:,:) = 1;
p_vec3(1,1,:) = p_vec;

% Auxiliary variables for covert-criterion task
if task == 2 || task == 3
    nx = 1001;  % Measurement grid size
    MAXSD = 5;  % When integrating a Gaussian go up to this distance
    X = bsxfun(@plus, S, sigma_ellipse*MAXSD*linspace(-1,1,nx));
    W = normpdf(linspace(-MAXSD,MAXSD,nx));
    W = W./qtrapz(W);
end

%% Begin loop over trials
P = zeros(NumTrials+1,Nprobs);
P(1,:) = ones(1,Nprobs)/Nprobs;
last = zeros(NumTrials,size(post,1));
PCx = [];

for t = 1:NumTrials
    %t
    if task == 2 || and(task == 3, mod(t,5) ~= 0); Xt = X(t,:); else Xt = []; end    
    [post,Psi,pi_post,PCxA] = bayesianOCPDupdate(Xt,C(t),post,Psi,Tmat,H,p_vec3,beta_hyp,mu,sigma,task,t);
    tt = nansum(post,2);

    % The predictive posterior is about the next trial
    P(t+1,:) = pi_post;
    last(t,:) = tt/sum(tt);
    
    % Record conditional posterior for covert-criterion task
    if task == 2 || and(task == 3, mod(t,5) ~=0)
        if isempty(PCx); PCx = zeros(NumTrials,size(PCxA,2)); end
        PCx(t,:) = PCxA;
    end
end

%% Compute log likelihood

MIN_P = 1e-4;   % Minimum lapse/error probability

% Compute predicted criterion
pA_t = sum(bsxfun(@times, P, p_vec(:)'),2);
pB_t = sum(bsxfun(@times, P, 1-p_vec(:)'), 2);
pA_t_bias = bsxfun(@plus, w*pA_t, (1-w)*p_bias);
pB_t_bias = bsxfun(@plus, w*pB_t, (1-w)*p_bias);
Gamma_t = pA_t_bias./pB_t_bias;
z_opt = sigma^2 * log(Gamma_t) / diff(mu);
z_opt = z_opt(1:NumTrials);
        
switch task
    case 1  % Overt-criterion task

        % Log probability of overt task responses
        log_Pz = -0.5*log(2*pi*sigma_criterion) - 0.5*((resp_obs-z_opt)./sigma_criterion).^2;
        if lambda > 0
            log_Pz = log(lambda/360 + (1-lambda)*exp(log_Pz));    
        end

        % Sum negative log likelihood
        log_Pz(~isfinite(log_Pz)) = log(MIN_P/360);
        nLL = -sum(log_Pz);
        resp_model = z_opt;
        
    case 2  % Covert-criterion task
        
        Pcx = 1./(1 + ((1-PCx)./PCx).^gamma);       % Softmax        
        PChatA = qtrapz(bsxfun(@times,Pcx,W),2);    % Marginalize over noise
        PChatA_bias = w*PChatA + (1-w)*p_bias;      % Conservative bias when w < 1
        lambda = max(1e-4,lambda);                  % Minimum lapse to avoid numerical trouble
        PChatA_bias = lambda/2 + (1-lambda)*PChatA_bias;
        
        % Log probability of covert task responses
        log_PChat = log(PChatA_bias).*(resp_obs == 1) + log(1-PChatA_bias).*(resp_obs ~= 1);
        
        % Sum negative log likelihood
        log_PChat(~isfinite(log_PChat)) = log(MIN_P);
        nLL = -nansum(log_PChat);   
        resp_model = PChatA_bias(1:NumTrials);
        
    case 3  % Mixed design task
        
        log_P = zeros(1, NumTrials);
        resp_model = zeros(1, NumTrials);
        % Log probability of overt task responses
        log_P(idx_Overt) = -0.5*log(2*pi*sigma_criterion) - 0.5*((resp_obs(idx_Overt)-z_opt(idx_Overt))./sigma_criterion).^2;
        if lambda > 0
            log_P(idx_Overt) = log(lambda/360 + (1-lambda)*exp(log_P(idx_Overt)));    
        end
        log_P(~isfinite(log_P)) = log(MIN_P/360);        
        resp_model(idx_Overt) = z_opt(idx_Overt);
        
        % Log probability of covert task responses
        PCx = PCx(idx_Covert,:);
        Pcx = 1./(1 + ((1-PCx)./PCx).^gamma);       % Softmax        
        PChatA = qtrapz(bsxfun(@times,Pcx,W),2);    % Marginalize over noise
        PChatA_bias = w*PChatA + (1-w)*p_bias;      % Conservative bias when w < 1
        lambda = max(1e-4,lambda);                  % Minimum lapse to avoid numerical trouble
        PChatA_bias = lambda/2 + (1-lambda)*PChatA_bias;
        
        log_P(idx_Covert) = log(PChatA_bias).*(resp_obs(idx_Covert) == 1) + log(1-PChatA_bias).*(resp_obs(idx_Covert) ~= 1);
        log_P(~isfinite(log_P)) = log(MIN_P);        
        resp_model(idx_Covert) = PChatA_bias;
        
        % Sum negative log likelihood
        nLL = -sum(log_P);
end

% RMSE between predictive posterior probability and true category probability
meanP = sum(bsxfun(@times,p_vec,P(1:NumTrials,:)),2);
p_estimate = meanP;
rmse = sqrt(mean((meanP - p_true).^2));

end

%--------------------------------------------------------------------------
function [post,Psi,pi_post,PCxA] = bayesianOCPDupdate(X,C,post,Psi,Tmat,H,p_vec,beta_hyp,mu,sigma,task,trial)
%BAYESIANCPDUPDATE Bayesian online changepoint detection update

    %if mod(t,100) == 0
    %    t
    %end
    L = size(H,1);  % Width of prior over run lengths
    % Slice with nonzero hazard function (probability of change)
    idxrange = (size(post,1)-L+1:size(post,1));

    % Posterior over pi_i before observing C_t
    pi_post = bsxfun(@times, Psi, Tmat);
    predCatA(:,:) = sum(bsxfun(@times, pi_post, p_vec),3);
    predCatB(:,:) = sum(bsxfun(@times, pi_post, 1 - p_vec),3);
    
    %----------------------------------------------------------------------
    % 2. Compute probability of response for covert-criterion task
    if task == 2 || (task == 3 && mod(trial,5) ~= 0)
        pxCatA = exp(-0.5*((X-mu(1))./sigma).^2);
        pxCatB = exp(-0.5*((X-mu(2))./sigma).^2);

        PCxA = nansum(nansum(bsxfun(@times,predCatA,post),2),1).*pxCatA;
        PCxB = nansum(nansum(bsxfun(@times,predCatB,post),2),1).*pxCatB;    
        PCxA = PCxA./(PCxA + PCxB);
    else
        PCxA = [];
    end
    
    %----------------------------------------------------------------------
    % 3. Observe Ct (do nothing)
    
    %----------------------------------------------------------------------
    % 4a. Evaluate predictive probability
    if C == 1
        predCat = predCatA ./ (predCatA + predCatB);
    else
        predCat = predCatB ./ (predCatA + predCatB);         
    end
    
    % 4b. Multiply posterior by predictive probability
    post = bsxfun(@times, post, predCat);
        
    % 4c. Evaluate posterior probability over state (only for relevant range)
    if C == 1
        pi_postC = bsxfun(@times, pi_post(idxrange,:,:), p_vec);
    else
        pi_postC = bsxfun(@times, pi_post(idxrange,:,:), 1-p_vec);
    end
    pi_postC = bsxfun(@rdivide, pi_postC, sum(pi_postC,3));

    %----------------------------------------------------------------------    
    % 5a. Calculate unnormalized changepoint probabilities
    slice = post(idxrange,:);   % Slice out trials where change is possible
    currenttrial = nansum( ...
        bsxfun(@times, sum(bsxfun(@times, slice, pi_postC),2), H), ...
        1);
    
    % 5b. Calculate unnormalized growth probabilities    
    post(idxrange,:) = bsxfun(@times, slice, 1-H); 
    
    % Shift posterior by one step
    post = circshift(post,1,1);
    post(1,:) = currenttrial;
    post(isnan(post)) = 0;
    
    % 5c. Calculate normalization
    Z = sum(post(:));
    
    % 5d. Normalize posterior
    post = post ./ Z;
    
    %----------------------------------------------------------------------
    % 6a. Update sufficient statistics
    Psi = circshift(Psi,1,1);
    if C == 1
        Psi = bsxfun(@times, Psi, p_vec);
    else
        Psi = bsxfun(@times, Psi, 1 - p_vec);
    end
    
    % Hyperprior
    Psi(1,1,:) = exp(beta_hyp(1).*log(p_vec) + beta_hyp(2).*log(1-p_vec));
    Psi = bsxfun(@rdivide, Psi, sum(Psi,3));
    
    % 6b. Store predictive posterior over pi_t
    pi_post = nansum(sum(bsxfun(@times,bsxfun(@times, Psi, Tmat), post),2),1);
    pi_post = pi_post / sum(pi_post);

end
