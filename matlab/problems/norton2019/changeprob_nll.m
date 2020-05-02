function [nLL, rmse, resp_model, p_estimate, post, v_estimate] = changeprob_nll(params,fitinfo)
%CHNAGEPROB_NLL Computes the negative log likelihood for the specified
%model
%(Documentation to be written)

idx_params = fitinfo.I_params;
inputParams = fitinfo.inputParams;
NumTrials = fitinfo.NumTrials;
mu = fitinfo.mu;
sigma_s = fitinfo.sigma_s;
C = fitinfo.C;
S = fitinfo.S;
p_true = fitinfo.p_true;
resp_obs = fitinfo.resp_obs;
task = fitinfo.task;
score = fitinfo.score;

inputParams(idx_params) = params;
sigma = sqrt(sigma_s^2 + inputParams(1)^2);
v_estimate = [];

if inputParams(7) == 0
    prior_rl = [];
elseif (inputParams(7) ~= 0) && (inputParams(14) == 0)
    prior_rl = [max(1,inputParams(7)*2/3),inputParams(7)];
else
    prior_rl = [(inputParams(7)-1), inputParams(7)-1+inputParams(14)];
end
if inputParams(8) == 0
    p_vec = [];
elseif (inputParams(8) == 1) && (inputParams(15) == 0)
    p_vec = linspace(inputParams(8), 1-inputParams(8), 5);
else
    p_vec = linspace(inputParams(8), 0.5+inputParams(15), 5);
end
if inputParams(9) == 0
    beta_hyp = [];
else
    beta_hyp = inputParams(9)^2;
end
[nLL, rmse, p_estimate, resp_model, post] = ...
    ChangeProb_bocpd_nll_v2(inputParams(1:6), NumTrials, mu, sigma, C, S, p_true, resp_obs, score, task, prior_rl, p_vec, beta_hyp);


end
