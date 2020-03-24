function [fmu,fs2] = sqrtgp_pred(gp,xs)
% SQRTGP_PRED gp prediction from sqrt gp representation used in WSABI

% Luigi Acerbi 2019

compute_var = nargout > 1 || lower(gp.method(1)) == 'm';

VV = gp.VV;
lambda = gp.lambda;
noise = gp.noise;
alpha = gp.alpha;

% xs = xs';
xxxScaled   = gp.X .* repmat(sqrt(1./VV),[size(gp.X,1),1]);

distxsxx    = ...
pdist2_squared_fast(xs.*repmat(sqrt(1./VV),[size(xs,1),1]),xxxScaled);

Kxsx        = lambda^2 * (1/(prod(2*pi*VV).^0.5)) * exp(-0.5*distxsxx);

l_0 = (Kxsx*(gp.invKxx*gp.y));

if compute_var
    Kxsxs = lambda^2 * (1/(prod(2*pi*VV).^0.5))*(1+noise)*ones(length(xs(:,1)),1);
    varPred = (Kxsxs - diag(Kxsx*(gp.invKxx * Kxsx')));

    if any(varPred <= 0) || any(isnan(xs(:))) || any(~isreal(varPred))
        %keyboard;
        varPred = 0;
    end
end

switch lower(gp.method)
    case 'l'
        fmu = alpha + 0.5*l_0.^2;
        if nargout > 1
            fs2 = l_0.^2.*varPred;
        end
    case 'm'
        fmu = alpha + 0.5*(l_0.^2 + varPred);
        if nargout > 1
            fs2 = l_0.^2.*varPred + 0.5*varPred.^2;
        end
end

end