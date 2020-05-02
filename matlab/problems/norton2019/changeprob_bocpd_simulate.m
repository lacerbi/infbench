function dataSim = changeprob_bocpd_simulate(data, task, model, parameters)
%CHANGEPROB_BOCPD_SIMULATE Simulates a Bayesian online changepoint detection observer.
% (Documentation to be written.)
%
% Author:   Elyse Norton
% Email:    elyse.norton@gmail.com
% Date:     Jan/11/2017

% Check input arguments

% Check for data file
if nargin < 1; data = []; task = 1; parameters = []; model = {'idealBayesian'}; end

% Determine task to simulate
if nargin < 2; task = 1; parameters = []; model = {'idealBayesian'}; end
if isempty(task); task = 1; end
if task ~= 1 && task ~= 2 && task ~=3; error('TASK can only be 1 (overt-criterion), 2 (covert-criterion), or 3 (mixed design).'); end

% Choose the Bayesian model to simulate
if nargin < 3; model = {'idealBayesian'}; parameters = []; end

% Parameter vector: #1 is SIGMA_ELLIPSE or SIGMA_CRITERION; #2 is rl_prior
% (subBayesian_rlprior), w (subBayesian_conservative), or pvec
% (subBayesian_pVec)
if nargin < 4; parameters = []; end
if isempty(parameters) && strcmp(model, 'idealBayesian') 
    parameters = 5; 
elseif isempty(parameters) && strcmp(model, 'subBayesian_rlprior')
    parameters = [5, 120]; 
elseif isempty(parameters) && strcmp(model, 'subBayesian_conservative')
    parameters = [5, 1];
elseif isempty(parameters) && strcmp(model, 'subBayesian_pVec') 
    parameters = [5, .2]; 
end

% Get session parameters
if isempty(data) 
    [NumTrials, sigma_ellipseData, mu, sigma_s, C, S, p_true, resp_obs] = changeprob_getSessionParameters([], task, [parameters(1), parameters(1)]); % Generate fake data using specified parameters
else
    [NumTrials, sigma_ellipseData, mu, sigma_s, C, S, p_true, resp_obs] = changeprob_getSessionParameters(data, task); % Determine session parameters from provided data set
end

% Observer model parameters
if task ~= 2
    sigma_criterion = parameters(1);
    sigma_ellipse = sigma_ellipseData;
else
    sigma_ellipse = parameters(1);
end
sigma = sqrt(sigma_s^2 + sigma_ellipse^2);
lambda = 1e-4;      % Default lapse (i.e., tiny lapse)
gamma = Inf;        % Default gamma (i.e., BDT)

if strcmp(model, 'idealBayesian')
    prior_rl = [80, 120];
    w = 1;
    p_vec = linspace(.2, .8, 5);
    beta_hyp = [0,0];
elseif strcmp(model, 'subBayesian_rlprior')
    prior_rl = [max(1,floor(parameters(2)*2/3)),parameters(2)];
    w = 1;
    p_vec = linspace(.2, .8, 5);
    beta_hyp = [0,0];
elseif strcmp(model, 'subBayesian_conservative')
    prior_rl = [80, 120];
    w = parameters(2);
    p_vec = linspace(.2, .8, 5);
    beta_hyp = [0,0];
elseif strcmp(model, 'subBayesian_pVec')
    prior_rl = [80, 120];
    w = 1;
    p_vec = linspace(parameters(2), 1-parameters(2), 5);
    beta_hyp = [0,0];
elseif strcmp(model, 'subBayesian_betahyp')
    prior_rl = [80, 120];
    w = 1;
    p_vec = linspace(.2, .8, 5);
    beta_hyp = parameters(2)*[1,1];
elseif strcmp(model, 'subBayesian_3param')
    prior_rl = [max(1,floor(parameters(2)*2/3)),parameters(2)];
    w = 1;
    p_vec = linspace(.2, .8, 5);
    beta_hyp = parameters(3)*[1,1];
end

%% Experiment constants

% Run length ~ Uniform[prior_rl(1),prior_rl(2)]
% if isfield(options,'prior_rl') && ~isempty(options.prior_rl)
%     prior_rl = options.prior_rl;
% else
%end
L = diff(prior_rl)+1;           % Number of elements

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

%% Initialize inference

% Hazard function (defined only where nonzero)
H = ones(L,1)/L;
H = H./flipud(cumsum(flipud(H)));

% Posterior over run lengths (from 0 to PRIOR_RL(2)-1)
post = zeros(prior_rl(2),Nprobs);
post(1,:) = 1;  % Change in the first trial
post = post ./ sum(post(:));

% Table of binomial count/probability
Psi = zeros(size(post,1),1,Nprobs); Psi(1,1,:) = 1/Nprobs;

% Transition matrix
Tmat(1,:,:) = (ones(Nprobs,Nprobs) - eye(Nprobs)) / (Nprobs-1);
% Tmat(1,:,:) = 1;
p_vec3(1,1,:) = p_vec;

% Auxiliary variables for covert-criterion task
if task ~= 1
    X = S + randn(size(S))*sigma_ellipse;
end
% if task == 2 || task == 3
%     nx = 1001;  % Measurement grid size
%     MAXSD = 5;  % When integrating a Gaussian go up to this distance
%     X = bsxfun(@plus, S, sigma_ellipse*MAXSD*linspace(-1,1,nx));
%     W = normpdf(linspace(-MAXSD,MAXSD,nx));
%     W = W./qtrapz(W);
% end

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

% Compute predicted criterion
pA_t = sum(bsxfun(@times, P, p_vec(:)'),2);
pB_t = sum(bsxfun(@times, P, 1-p_vec(:)'), 2);
pA_t_bias = bsxfun(@plus, w*pA_t, (1-w)*p_bias);
pB_t_bias = bsxfun(@plus, w*pB_t, (1-w)*p_bias);
Gamma_t = pA_t_bias./pB_t_bias;
z_opt = sigma^2 * log(Gamma_t) / diff(mu);
z_opt = z_opt(1:NumTrials);
resp_obs = zeros(size(z_opt));
score = resp_obs;
        
switch task
    case 1  % Overt-criterion task
        
        resp_obs = z_opt + randn(size(z_opt))*sigma_criterion;
        score = double(or(C == 1 & S <= resp_obs, C == 0 & S > resp_obs));
        
    case 2  % Covert-criterion task
        
        Pcx = 1./(1 + ((1-PCx)./PCx).^gamma);        % Softmax        
        %PChatA = qtrapz(bsxfun(@times,Pcx,W),2);    % Marginalize over noise
        %PChatA_bias = w*PChatA + (1-w)*p_bias;      % Conservative bias when w < 1
        PChatA_bias = w*Pcx + (1-w)*p_bias;
        PChatA_bias = lambda/2 + (1-lambda)*PChatA_bias;
        
        idx_A = find(rand(size(PChatA_bias)) <= PChatA_bias);
        resp_obs(idx_A) = 1;
        score = double(C == resp_obs);
        
    case 3  % Mixed design task
        
        resp_obs(idx_Overt) = z_opt(idx_Overt) + randn(size(z_opt(idx_Overt)))*sigma_criterion;
        score(idx_Overt) = double(or(C(idx_Overt) == 1 & S(idx_Overt) <= resp_obs(idx_Overt), C(idx_Overt) == 0 & S(idx_Overt) > resp_obs(idx_Overt)));
        
        % Log probability of covert task responses
        PCx = PCx(idx_Covert,:);
        Pcx = 1./(1 + ((1-PCx)./PCx).^gamma);       % Softmax        
        %PChatA = qtrapz(bsxfun(@times,Pcx,W),2);    % Marginalize over noise
        %PChatA_bias = w*PChatA + (1-w)*p_bias;      % Conservative bias when w < 1
        PChatA_bias = w*Pcx + (1-w)*p_bias;
        PChatA_bias = lambda/2 + (1-lambda)*PChatA_bias;
        
        idx_A = find(rand(size(PChatA_bias)) <= PChatA_bias);
        
        resp_obs(idx_Covert(idx_A)) = 1;
        score(idx_Covert) = double(C(idx_Covert) == resp_obs(idx_Covert));
end

dataSim.NumTrials = NumTrials;
dataSim.mu = mu;
dataSim.sigma_s = sigma_s;
dataSim.sigma = sigma;
dataSim.sigmaEllipse = sigma_ellipse;
dataSim.Category = C;
dataSim.Stimulus = S;
dataSim.pA = p_true;
dataSim.response = resp_obs;
dataSim.score = score;

end

%--------------------------------------------------------------------------
function [post,Psi,pi_post,PCxA] = bayesianOCPDupdate(X,C,post,Psi,Tmat,H,p_vec,beta_hyp,mu,sigma,task, trial)
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
    Psi(1,1,:) = 1;
    Psi = bsxfun(@rdivide, Psi, sum(Psi,3));
    
    % Hyperprior
    Psi(1,1,:) = exp(beta_hyp(1).*log(p_vec) + beta_hyp(2).*log(1-p_vec));
    Psi = bsxfun(@rdivide, Psi, sum(Psi,3));
    
    % 6b. Store predictive posterior over pi_t
    pi_post = nansum(sum(bsxfun(@times,bsxfun(@times, Psi, Tmat), post),2),1);
    pi_post = pi_post / sum(pi_post);

end