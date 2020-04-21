function [NEV] = expectedVarL(xs, s2hat, lambda, VV, lHatD, xxxScaled, invKxx, noise, bb, BB, aa)
% Active sampling loss function for wsabi_L

xs = xs';

distxsxx    = ...
pdist2_squared_fast(xs.*repmat(sqrt(1./VV),length(xs(:,1)),1),xxxScaled);

Kxsx        = lambda^2 * (1/(prod(2*pi*VV).^0.5)) * exp(-0.5*distxsxx);
          
Kxsxs = lambda^2 * (1/(prod(2*pi*VV).^0.5))*(1+noise)*ones(length(xs(:,1)),1);

varPred = (Kxsxs - diag(Kxsx*(invKxx * Kxsx')));

l_0 = (Kxsx*(invKxx*lHatD));
if any(varPred <= 0) || any(isnan(xs(:))) || any(~isreal(varPred))
    %keyboard;
    varPred = 0;
end

priorWeighting = mvnpdf(xs,bb,BB);

%Negative expected variance
if isscalar(s2hat)
    NEV(1,:) = -(l_0.^2 .* varPred) .* (1 - s2hat./(s2hat + varPred)) .* priorWeighting.^2;    
else
    % Estimate observation noise at test points from nearest neighbor
    [~,pos] = min(pdist2_squared_fast(bsxfun(@rdivide,xs,sqrt(VV)),xxxScaled),[],2);
    sn2 = s2hat(pos);
    NEV(1,:) = -(l_0.^2 .* varPred) .* (1 - sn2./(sn2 + varPred)) .* priorWeighting.^2;
end

end