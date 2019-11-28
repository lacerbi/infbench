function [d_var,new_var] = gp_lookahead_var_fast_pend(gp,cur_var,x_tr,y_tr,x_eval,x_pend,P)
% Computes the GP 'lookahead' variance and the reduction of current variance when there
% are fixed pending points 'x_pend'. 
% This is used for greedy acquisition where we evaluate where the current uncertainty is 
% highest, i.e. MAXV and MAXIQR.
%
% INPUT:
% cur_var - current variance based on training data (x_tr,y_tr) evaluated at x_s
% x_tr - current training data (for which we already know corresponding y-values in y_tr)
% y_tr - data realisations for 'x_tr' evaluation locations
% x_eval - where the variances are evaluated (can be matrix, each row is one location)
% x_pend - pending locations (can be matrix, each row is one pending point)
% P - some precomputed GP terms (to save computation time)
%
% OUTPUT:
% d_var - the reduction in the variance at x_eval i.e. d_var = var_cur - var_new
% (can be matrix)
% new_var - new variance at x_eval after x_pend has been taken into account (can be matrix)

L_new = P.L_new;

% TODO: x_pred-part of the following could be precomputed to save a little comput. time:
cov_snew = gp_pred_cov_fast(gp,x_tr,y_tr,x_eval,x_pend,P,0);
cov_news = cov_snew';

d_var = sum(cov_snew.*((L_new)'\(L_new\cov_news))',2);
new_var = cur_var - d_var;
end

