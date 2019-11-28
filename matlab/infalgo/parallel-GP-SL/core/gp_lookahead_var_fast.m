function [d_var,new_var] = gp_lookahead_var_fast(gp,cur_var,x_tr,y_tr,sigma_tr,x_s,x_new,P,is_pend)
% Computes the GP 'lookahead' variance and the reduction of current variance. 
%
% INPUT:
% cur_var - current variance based on training data (x_tr,y_tr) evaluated at x_s
% x_tr - current training data (for which we already know corresponding y-values in y_tr)
% y_tr - data realisations for 'x_tr' evaluation locations
% x_s - where the variances are evaluated (can be matrix, each row is one location)
% x_new - new evaluation locations (can be matrix, each row is a location in the batch)
% P - some precomputed GP terms (to save computation time)
% is_pend - if == 1, then all the points in 'x_new' except the last one are treated as
% pending points
%
% OUTPUT:
% d_var - the reduction in the variance at x_s i.e. d_var = var_cur - var_new
% (can be matrix)
% new_var - new variance at x_s after x_new has been taken into account (can be matrix)

if is_pend ~= 1
    %% No pending points
    cov_new = gp_pred_cov_fast(gp,x_tr,y_tr,x_new,x_new,P,0);
    
    cov_new = cov_new + gp_noise_model_var(gp, x_tr, sigma_tr) * eye(size(cov_new));
    L_new = chol(cov_new,'lower');
    
    cov_snew = gp_pred_cov_fast(gp,x_tr,y_tr,x_s,x_new,P,1);
    cov_news = cov_snew';
    
    % compute only diagonals
    % naive computation:
    %d_var = diag(cov_snew*((L_new)'\(L_new\cov_news)));
    
    d_var = sum(cov_snew.*((L_new)'\(L_new\cov_news))',2);
    new_var = cur_var - d_var;
    
else
    %% Pending points, some precomputed values related to those should be provided in 'P'
    %x_new1 = x_new;
    x_pend = x_new(1:end-1,:); % all these are assumed to be pending points here
    x_new = x_new(end,:);
    
    cov_snew = gp_pred_cov_fast(gp,x_tr,y_tr,x_s,x_new,P,1);
    
    %cov_new = gp_pred_cov_fast(gp,x_tr,y_tr,x_new,x_new,P,0); % assumed 1 x 1 here
    [~,cov_new] = gp_pred_fast(gp,x_tr,y_tr,x_new,P);
    
    cov_new = cov_new + gp_noise_model_var(gp, x_tr, sigma_tr);
    
    cov_new_pend = gp_pred_cov_fast(gp,x_tr,y_tr,x_new,x_pend,P,0);
    cov_pend_new = cov_new_pend';
    
    cc = (P.Lpend)'\(P.Lpend\cov_pend_new);
    d_var_new = (cov_snew - P.cov_spend*cc).^2/(cov_new - cov_new_pend*cc);
    
    d_var = P.d_var_pend + d_var_new;
    new_var = cur_var - d_var;
end
end




