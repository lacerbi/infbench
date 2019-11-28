function sigma2 = gp_noise_model_var(gp, th_tr, sigma_tr, theta, tol)
% Returns the noise variance of the GP model at input parameter theta. If standard GP with
% constant noise parameter, then it is returned.

if nargin < 4
    theta = NaN;
end
if nargin < 5
    tol = 1e-9;
end

if strcmp(gp.lik.type, 'Gaussian-smt')
    % Special noise model
    % Estimates for noise at theta could be modelled and computed from a separate GP but
    % we now assume that the noise is actually negligible although this typically does not
    % quite hold in practice. 
    sigma2 = tol * ones(size(theta));
else
    % standard GP case: constant noise
    sigma2 = gp.lik.sigma2 * ones(size(theta));
end
end



