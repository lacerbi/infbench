function [gp,gp_optim_opt] = fit_gp_model(gp, gp_optim_opt, gp_opt, th_grid, ...
    log_post_tr, theta_tr, sigma_tr)
% Optimizes the GP hyperparameters using GPstuff 'gp_optim'.
%
% NOTE: GP_OPTIM_OPT are the settings for GP hyperparam optimisation, 
% GP_OPT are the general settings structure used for specifying GP model settings etc.

if nargin < 7
    sigma_tr = [];
end

%% Set up the GP model priors, initial settings etc.
init = isempty(gp);
if ~init && gp_opt.noise_model
    gp_old = gp;
end

if init || gp_opt.noise_model 
    % NOTE: if gaussian_smt, then variances need to be updated. This is done in a bit
    % clumsy way: create new gp and then copy only the updated 'lik'-part back to the 
    % original gp structure. This way e.g. the old hyperparam estimates are used as new 
    % initial points for hyperparameter optimisation.
    [gp,gp_optim_opt] = init_gp(gp_opt, th_grid, log_post_tr, theta_tr, sigma_tr);
end

if ~init && gp_opt.noise_model
    gp_old.lik = gp.lik;
    gp = gp_old;
end

%% fit GP hyperparameters
try
    [gp,fval,exitflag] = gp_optim(gp, theta_tr, log_post_tr(:), 'opt', gp_optim_opt);
catch
    % Optimisation failed, try again by starting from default settings
    [gp,gp_optim_opt] = init_gp(gp_opt, th_grid, log_post_tr, theta_tr, sigma_tr);
    [gp,fval,exitflag] = gp_optim(gp, theta_tr, log_post_tr(:), 'opt', gp_optim_opt);
    % If this still fails, let the error happen 
    
    %error('Optimisation of GP hypers failed due to code error.');
    %warning('Optimisation of GP hypers failed due to code error.');
    %keyboard;
end

%% debugging
if ~strcmp(gp_opt.display_type,'off')
    [w,s] = gp_pak(gp);
    disp('Hyperparameters:');
    disp(s);
    %disp('Values of log(hyperparameters):')
    %disp(w);
    disp('Values of hyperparameters:');
    disp(exp(w(1))); 
    if gp_opt.noise_model
        disp(exp(w(2:end)));
    else
        disp(exp(w(2:end-1)));
        disp(exp(w(end)));
    end
    disp(' ');
end
end


function [gp, gp_optim_opt] = init_gp(gp_opt, th_grid, log_post_tr, theta_tr, sigma_tr)
% Sets up initial GP structure and its default init values. 
%
% We set up rather noninformative priors for the GP model. The lengthscale is set 
% according to the bounds of the parameter space. The signal variance is set acording to
% 'typical' magnitudes of loglik. The noise variance is given a very wide prior (needed 
% only when gp_opt.noise_model is not used). 

% TODO: better default init values?
th_range = th_grid.range(:,2)-th_grid.range(:,1);

%% set up GP model etc.
pl = prior_t('s2',(th_range(1)/2)^2); % we suppose parameters are scaled similarly
ydata_var = 1000^2;
%ydata_var = 100^2;
%pm = prior_unif();
pm = prior_sqrtt('s2',ydata_var);
%pn = prior_unif();
pn = prior_sqrtt('s2',50^2);

gpcf = gpcf_sexp('lengthScale', th_range(:)'/3, 'magnSigma2', ydata_var, ...
    'lengthScale_prior', pl, 'magnSigma2_prior', pm);

switch gp_opt.noise_model        
    case 0
        sigma_n = 10;
        lik = lik_gaussian('sigma2', sigma_n^2, 'sigma2_prior', pn);
    case 1
        n_tr = length(log_post_tr);
        lik = lik_gaussiansmt('ndata', n_tr, 'sigma2', sigma_tr.^2);
        % Note: there is no sigma_n parameter to infer!
    case 2
        sigma_n = max(gp_opt.sigma_n_const,1e-3);
        lik = lik_gaussian('sigma2', sigma_n^2, 'sigma2_prior', prior_fixed()); % sigma fixed
end

if gp_opt.meanf
    % set different mean functions
    gpmf1 = gpmf_constant('prior_mean',0,'prior_cov',30^2);
    gpmf2 = gpmf_linear('prior_mean',0,'prior_cov',30^2);
    gpmf3 = gpmf_squared('prior_mean',0,'prior_cov',30^2);
    
    gp = gp_set('lik', lik, 'cf', gpcf, 'meanf', {gpmf1,gpmf2,gpmf3}, 'jitterSigma2', 1e-9);
    %gp = gp_set('lik', lik, 'cf', gpcf, 'meanf', {gpmf1}, 'jitterSigma2', 1e-9); % only const term
else
    % use the zero mean function
    gp = gp_set('lik', lik, 'cf', gpcf, 'jitterSigma2', 1e-9);
end
gp_optim_opt = optimset('TolFun',1e-3,'TolX',1e-3','display', 'off');

end



