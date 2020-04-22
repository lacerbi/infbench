function [] = sl_bootvar_test()
% Compute the logSL variance of the Ricker, Lorenz, G-and-k -simulation models at their
% true parameter values and compare to bootstrapped estimate.
% This is a simple investigation, more thoroughful analysis could be also carried out.
%
% Note: Computing these take some time and the results are not saved to file currently.

clc; close all; format('compact');

test_models = {'ricker','lorenz','gk_model'};
Ns = [100, 100, 100]; % same values as in actual simulations!

reps = 200; % how many repeated set of N summaries for each test model to approx variance
sl_method = 'sl'; % 'sl' or 'lfire' ('lfire' not implemented)
sl_opt.estimator = 'sl'; %'sl','ubsl','ublogsl'

nmodels = length(test_models);
for t = 1:nmodels
    % get sim model settings etc.
    [~,sim_model] = get_test_model(test_models{t},[],Ns(t));
    true_theta = sim_model.true_theta;
    
    sl_opt.N = Ns(t);
    
    % compute...
    disp('Computing...');
    llik = NaN(1,reps);
    bootvar = NaN(1,reps);
    for k = 1:reps
    	[llik(k),bootvar(k)] = noisy_loglik_estim(sim_model,sl_opt,true_theta,sl_method,1);
    end
    true_var = var(llik);
    disp('Done.');
    
    % print results
    model = test_models{t}
    true_var
    bootvar=sort(bootvar)
    mean_bootvar = mean(bootvar)
    med_bootvar = median(bootvar)
    disp('----------------------');
end
end




