function [ll,Pr_mat] = krajbich2010_loglike(params,data)
%KRAJBICH2010_GENDATA Log likelihood for Krajbich et al. model

% Assign parameters from parameter vector
sigma = params(1)*sqrt(data.dt);
d = params(2)*data.dt;
theta = params(3);
lambda = params(4);

[Ntrials,Nbins] = size(data.fixmat);

% Grid size
if isfield(data,'Ng') && ~isempty(data.Ng)
    Ng = data.Ng;
else
    Ng = 301;
end

grid = linspace(-1+0.5/Ng,1-0.5/Ng,Ng)';
dg = 2/Ng;

leftval = data.leftval;
rightval = data.rightval;
fixmat = data.fixmat;

Pr_left = zeros(1,Nbins);
Pr_right = zeros(1,Nbins);

ll = zeros(Ntrials,1);

for iTrial = 1:Ntrials

    p = zeros(Ng,1);
    p((Ng+1)/2,1) = 1/dg;     % All probability mass at 0    

    drift_vec = d*(bsxfun(@times,leftval(iTrial) - theta*rightval(iTrial),fixmat(iTrial,:) == 1) ...
    + bsxfun(@times,theta*leftval(iTrial) - rightval(iTrial),fixmat(iTrial,:) == 2));

    shift_mat_left = bsxfun_normpdf(grid+d*(leftval(iTrial)-theta*rightval(iTrial)),grid',sigma);
    shift_mat_right = bsxfun_normpdf(grid+d*(theta*leftval(iTrial)-rightval(iTrial)),grid',sigma);

    nf = 1;
    for iTime = 1:Nbins
        Pr_left(iTime) = sum(bsxfun(@times,p,1 - bsxfun_normcdf(1,grid+drift_vec(iTime),sigma)),1)*dg;
        Pr_right(iTime) = sum(bsxfun(@times,p,bsxfun_normcdf(-1,grid+drift_vec(iTime),sigma)),1)*dg;
        nf = nf - Pr_left(iTime) - Pr_right(iTime);
        if fixmat(iTrial,iTime) == 1
            p = sum(bsxfun(@times,p,shift_mat_left),1)'*dg;            
        else
            p = sum(bsxfun(@times,p,shift_mat_right),1)'*dg;
        end
        p = nf*p/(sum(p)*dg);
    end

    if data.choice(iTrial) == 1
        Presp = Pr_left(data.totfixdurbin(iTrial));
    else
        Presp = Pr_right(data.totfixdurbin(iTrial));
    end

    if nargout > 1
        Pr_mat(iTrial,:,1) = Pr_left;
        Pr_mat(iTrial,:,2) = Pr_right;        
    end
    
    Presp(~isfinite(Presp)) = 0;
    Presp = (1-lambda)*Presp + lambda/2/Nbins;    
    ll(iTrial) = log(Presp);    
end

ll = sum(ll);

if nargout > 1
    Pr_mat = (1-lambda)*Pr_mat + lambda/2/Nbins;
end


end

function p = bsxfun_normcdf(x,mu,sigma)
%BSXFUN_NORMCDF Vectorized normal cumulative distribution function (cdf).
%   P = BSXFUN_NORMCDF(X,MU,SIGMA) returns the cdf of the normal 
%   distribution with mean MU and standard deviation SIGMA, evaluated at 
%   the values in X. Dimensions of X, MU, and SIGMA must either match, or 
%   be equal to one. Computation of the pdf is performed with singleton
%   expansion enabled via BSXFUN. The size of Y is the size of the input 
%   arguments (expanded to non-singleton dimensions).
%
%   All elements of SIGMA are assumed to be non-negative (no checks).
%
%   See also BSXFUN, BSXFUN_NORMPDF, NORMCDF.

%   Author: Luigi Acerbi
%   Release date: 15/07/2015

if nargin<3
    error('bmp:bsxfun_normcdf:TooFewInputs','Input argument X, MU or SIGMA are undefined.');
end

try
    if isscalar(mu)
        z = bsxfun(@rdivide, x-mu, sigma);
    elseif isscalar(sigma)
        z = bsxfun(@minus, x, mu)./sigma;
    else
        z = bsxfun(@rdivide, bsxfun(@minus, x, mu), sigma);
    end
    
    % Use the complementary error function, rather than .5*(1+erf(z/sqrt(2))), 
    % to produce accurate near-zero results for large negative x.
    p = 0.5 * erfc(-z ./ sqrt(2));

catch
    error('bmp:bsxfun_normpdf:InputSizeMismatch',...
          'Non-singleton dimensions must match in size.');
end

end

function y = bsxfun_normpdf(x,mu,sigma)
%BSXFUN_NORMPDF Vectorized normal probability density function (pdf).
%   Y = BSXFUN_NORMPDF(X,MU,SIGMA) returns the pdf of the normal 
%   distribution with mean MU and standard deviation SIGMA, evaluated at 
%   the values in X. Dimensions of X, MU, and SIGMA must either match, or 
%   be equal to one. Computation of the pdf is performed with singleton
%   expansion enabled via BSXFUN. The size of Y is the size of the input 
%   arguments (expanded to non-singleton dimensions).
%
%   All elements of SIGMA are assumed to be non-negative (no checks).
%
%   See also BSXFUN, BSXFUN_NORMCDF, NORMPDF.

%   Author: Luigi Acerbi
%   Release date: 15/07/2015

if nargin<3
    error('bmp:bsxfun_normpdf:TooFewInputs','Input argument X, MU or SIGMA are undefined.');
end

try
    if isscalar(mu)
        y = bsxfun(@rdivide, exp(-0.5*bsxfun(@rdivide, x - mu, sigma).^2), sigma)/sqrt(2*pi);
    elseif isscalar(sigma)
        y = exp(-0.5*(bsxfun(@minus, x, mu)/sigma).^2)/(sigma*sqrt(2*pi));
    else
        y = bsxfun(@rdivide, exp(-0.5*bsxfun(@rdivide, bsxfun(@minus, x, mu), sigma).^2), sigma)/sqrt(2*pi);
    end
catch
    error('bmp:bsxfun_normpdf:InputSizeMismatch',...
          'Non-singleton dimensions must match in size.');
end

end
