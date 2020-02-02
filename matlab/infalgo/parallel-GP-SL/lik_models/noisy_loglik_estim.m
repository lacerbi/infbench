function [llik,bootvar] = noisy_loglik_estim(sim_model,opt,theta,method,use_boot)
% Computes a noisy log likelihood estimate by simulating N datasets and using SL 
% approximation. Also estimate of the variance of the lik value can be computed using the
% bootstrap. 
% TODO: LFIRE

if nargin < 4
    method = 'sl';
end
if nargin < 5
    use_boot = 1;
end
bootvar = NaN;

%% (noisy) log lik can be computed without simulation runs:
if strcmp(method,'exact')
    llik = sim_model.loglik_eval(theta);
    return;
elseif strcmp(method,'noisy_exact')
    [llik,llik_std] = sim_model.loglik_eval(theta);
    bootvar = llik_std.^2;
    return;
end

%% SL or LFIRE case:
N = opt.N;
summaries_th = NaN(sim_model.summary_dim,N); % each column one summary vector

% repeated sampling for SL (this could be parallelized in practice)
for i = 1:N
    data_i = sim_model.gen(theta, sim_model.n_data);
    summaries_th(:,i) = sim_model.comp_summaries(data_i,sim_model.data);
end

% compute SL or LFIRE log-likelihood value
if ~strcmp(method,'lfire')
    llik = eval_sl_loglik(sim_model.summary_true, summaries_th, opt.estimator);
    
    if use_boot && nargout > 1
        % uses bootstrap to compute estimate of the variance of the loglik estimate
        bootn = 2000;
        llik_boot = NaN(bootn,1);
        for i = 1:bootn
            boot_inds = randsample(N,N,1); % sample with replacement
            s_thi = summaries_th(:,boot_inds);
            llik_boot(i) = eval_sl_loglik(sim_model.summary_true, s_thi, opt.estimator);
        end
        bootvar = var(llik_boot); % use robust estimate here instead?
        
        if 0
            % show bootstrap samples
            figure(25);
            histogram(llik_boot,30);
            title('Bootstrap samples of SL loglik');
            xlabel('SL loglik values');
            xlim([min(llik_boot) max(llik_boot)]);
            pause;
        end
    end

else % lfire
    error('Not implemented.');
end


%% Some individual possibly very large loglik values occurring near some boundaries can 
%% make fitting the GP problematic so we deal with those cases by lower bounding the loglik
%trunclik = -Inf;
trunclik = -10^5;
llik = max(trunclik,llik);
bootvar(llik<=trunclik) = 100^2;
%llik = soft_loglik_lb(llik);

end







