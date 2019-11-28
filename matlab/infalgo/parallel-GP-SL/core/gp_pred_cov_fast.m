function covft = gp_pred_cov_fast(gp,x_tr,y_tr,x_s1,x_s2,P,K_precomput)
% Takes a GP structure together with matrix x_tr of training inputs and vector y_tr of 
% training targets, and evaluates the posterior covariance between test points x_s1 and
% x_s2. (This should be faster than computing the whole joint covariance matrix of input 
% [x_s1;x_s2] and then extracting the required covariances.)
%
% TODO: could handle separately the specific case where x_s1==x_s2 to save some more 
% comput. time

K_s1_s2 = gp_cov(gp, x_s1, x_s2, []);
K_tr_s2 = gp_cov(gp, x_tr, x_s2, []);

if K_precomput 
    % 'K_s1_trK' already computed and provided in 'P'
    %K_s1_trK = P.K_s1_trK;
    covft = P.K_s1_trK * K_tr_s2;
else
    % compute now
    K_s1_tr = gp_cov(gp, x_s1, x_tr, []); 
    covft = K_s1_tr*((P.L)'\(P.L\K_tr_s2));
end

covft = K_s1_s2 - covft;

% If there are specified mean functions
if isfield(gp,'meanf')
    Hs1 = mean_prep_xs(gp,x_s1);
    Hs2 = mean_prep_xs(gp,x_s2);
    
    L = P.L;
    if K_precomput 
        Rs1 = Hs1 - P.HK_s1_trKT;
    else
        Ks1K = L'\(L\K_s1_tr'); 
        Rs1 = Hs1 - P.H*Ks1K;
    end
    Ks2K = L'\(L\K_tr_s2);  
    Rs2 = Hs2 - P.H*Ks2K;
    
    invARs2 = P.B1\Rs2;
    RARapu = Rs1'*invARs2;
    RAR = RARapu; 
	covft = covft + RAR; 
end
end


