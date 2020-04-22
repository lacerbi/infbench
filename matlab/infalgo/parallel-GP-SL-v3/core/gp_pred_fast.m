function [eft,varft] = gp_pred_fast(gp,x_tr,y_tr,x_s,P)
% Computes GP mean and variance at input points 'x_s'. Unlike corresponding GPstuff 
% function gp_pred, this version assumes that we have standard GP model and that some 
% quantities in the GP formulas have been precomputed and are provided in the input struct
% P. No input checkings etc. are done. 

%% GPstuff gp_pred
% just calls the GPstuff prediction function, this is slow because cov-matrix is always
% inverted even if here it could be precomputed
%[eft, varft] = gp_pred(gp, x_tr, y_tr(:), x_s); 
%return;

%% fast GP prediction (thanks to precomputed values and lack of input checkings)
K_tr_s = gp_cov(gp, x_tr, x_s, []);
eft = K_tr_s' * P.a;

if isfield(gp,'meanf')
    % terms with non-zero mean -prior
    Hs = mean_prep_xs(gp,x_s);
    L = P.L;
    K_nf = K_tr_s;
    KsK = L'\(L\K_nf);  
    Rs = Hs - P.H*KsK;
    RB  = Rs'*P.Beta;  
    %[RB, RAR] = mean_predf(gp, x_tr ,x_s, K_tr_s, P.L, P.a);
    eft = eft + RB;
end

V = gp_trvar(gp,x_s,[]);
v = P.L\K_tr_s;
varft = V - sum(v'.*v',2);

% If there are specified mean functions
if isfield(gp,'meanf')
    B1 = P.B1;
    invARs = B1\Rs;
    RARapu = Rs'.*invARs'; 
    RAR = sum(RARapu,2); 
    varft = varft + RAR;
end
end




