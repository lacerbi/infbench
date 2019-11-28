function d_cov = gp_lookahead_cov_fast(gp,x_tr,y_tr,sigma_tr,x_s1,x_new,P)
% Computes the full GP 'lookahead' covariance matrix at 'x_s1'. 
% This is a quickly made implementation needed only for testing purposes. 
%
% INPUT:
% x_tr - current training data (for which we already know corresponding y-values in y_tr)
% y_tr - data realisations for 'x_tr' evaluation locations
% x_s1 - where the covariance matrix is evaluated (can be matrix, each row is one location)
% x_new - new evaluation locations (can be matrix, each row is a location in the batch)
% P - some precomputed GP terms (to save computation time)
%
% OUTPUT:
% d_cov - the reduction of the full covariance at x_s1 i.e. d_cov = cov_cur - cov_new

cov_new = gp_pred_cov_fast(gp,x_tr,y_tr,x_new,x_new,P,0);

cov_new = cov_new + gp_noise_model_var(gp, x_tr, sigma_tr) * eye(size(cov_new));
L_new = chol(cov_new,'lower');

cov_s1new = gp_pred_cov_fast(gp,x_tr,y_tr,x_s1,x_new,P,1);
cov_news1 = cov_s1new';

d_cov = cov_s1new*((L_new)'\(L_new\cov_news1)); % computes now the whole cov-matrix
end




