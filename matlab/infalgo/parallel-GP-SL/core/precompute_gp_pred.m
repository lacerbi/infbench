function P = precompute_gp_pred(gp, x_tr, y_tr, gp_opt)
% Precompute some values for gp_pred_fast.

%P = [];
%return;

% standard GP case:
[K, C] = gp_trcov(gp, x_tr);
L = chol(C,'lower');
a = L'\(L\y_tr);

P.L = L;
P.a = a;

if isfield(gp,'meanf')
    % some terms in mean function case can be precomputed since they don't depend on the
    % test points
    [H,b,B] = mean_prep1(gp,x_tr);
    KsH = L'\(L\H');
    invB = B\eye(size(B));
    Ksy = a;
    B1 = invB + H*KsH;
    B2 = H*Ksy + invB*b;
    Beta = B1\B2;
    
    P.H = H;
    P.B1 = B1;
    P.Beta = Beta;
end
end




