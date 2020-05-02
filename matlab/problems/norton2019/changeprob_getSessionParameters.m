function [NumTrials,sigma_ellipseData,mu,sigma_s,C,S,p_true,resp,score] = changeprob_getSessionParameters(data,task,parameters)
%CHANGEPROB_GETSESSIONPARAMETERS Gets session parameters from an existing
%data struct or creates a fake dataset

%   INPUT: 
        % data
            % Experimental data struct to be decomposed
        % task: 1 - overt (default), 2 - covert, 3 - mixed
        % parameters used to generate fake data
            % parameters(1): sensory noise (sigma_ellipse)
            % parameters(2): adjustment noise (sigma_criterion)
        
%   OUTPUT:
        % NumTrials: total number of trials
        % sigma_ellipseData: sensory noise from calibration data
        % mu: vector containing the category means [muA, muB]
        % sigma_s: std dev of the category distributions
        % C: vector of category values (1 - A, 0 - B)
        % S: vector of true stimulus orientations
        % p_true: vector containing the probability of A
        % resp: vector containing the observer's criterion setting 
        % (overt task), categorizations (covert), or both (mixed)
        % score: 0 - wrong, 1 - correct
        
%   Author: Elyse Norton
%   Date: 1/11/17
%   email: elyse.norton@gmail.com
        
sigma_ellipse = parameters(1);
sigma_criterion = parameters(2);

col = data.SessionOrder(task);     % Column of overt-criterion task
NumTrials = data.NumTrials;    

% Noise parameters
sigma_s = data.StdDev;
if isempty(sigma_ellipse) || ~isfinite(sigma_ellipse)
    sigma_ellipseData = data.EllipseNoise;
    sigma = sqrt(sigma_s^2 + sigma_ellipseData^2);
else
    sigma = sqrt(sigma_s^2 + sigma_ellipse^2);
end

% Category information 
% (in the data Category B/Red is coded as 1, Category A/Green is coded as 2)
C = (data.Category(:,col) == 2);        % Category A/Green
C = double(C);
p_true = data.pA(:,col);
mu = [data.GreenMean(col),data.RedMean(col)];

% Shift coordinate system to zero
mu_bar = mean(mu);
mu = bsxfun(@minus, mu, mu_bar);
S = bsxfun(@minus, data.StimulusAngle(:,col), mu_bar);

% Get task-relevant responses
switch task
    case 1  % Overt-criterion task    
        resp = bsxfun(@minus, data.Criterion(:,col), mu_bar);   % Reported criterion
    case 2  % Covert-criterion task
        resp = data.Response(:,col) == 2;
        resp = double(resp);
end
score = data.Score(:,col);
    
end

