function nelbo = negelbo_mfvi(phi,fun,Ns)
%ELBO_MFVI

if nargin < 3 || isempty(Ns); Ns = 100; end

D = size(phi,2)/2;

mu = phi(1,1:D);
ln_sigma = (phi(1,(1:D)+D));
sigma = exp(ln_sigma);

if isscalar(Ns)
    xx = bsxfun(@plus,mu,bsxfun(@times,randn(Ns,D),sigma));
else
    xx = bsxfun(@plus,mu,bsxfun(@times,Ns,sigma));
    Ns = size(xx,1);
end

logp = zeros(Ns,1);
for i = 1:Ns
    logp(i) = fun(xx(i,:));    
end

H = D*0.5*(1+log(2*pi)) + sum(ln_sigma);

nelbo = -(sum(logp)/Ns + H);

end