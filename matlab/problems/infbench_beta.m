function [y,y_std] = infbench_beta(x,infprob,mcmc_params)
%INFBENCH_BETA Inference benchmark log pdf -- Beta distribution with different parameters.

if nargin < 3; mcmc_params = []; end

if isempty(x)
    
    % Initialization call -- define problem and set up data
    n = infprob(1);

    % Are we transforming the entire problem to unconstrained space?
    transform_to_unconstrained_coordinates = false;

    % Assign dimension D, alpha and beta parameters based on n
    % Lowest digit (units) a sets alpha = 2^(a-1)
    % Second-lowest digit (tens) b sets beta = 2^(b-1)
    % Top digit (hundreds) + 1 is D
    a = mod(n,10);
    b = mod(floor(n/10),10);
    alpha = 2^(a-1);
    beta = 2^(b-1);
    data.beta_params = [alpha, beta];
    D = 1 + floor(n/100);
        
    % Same Beta distribution duplicated in multiple dimensions D (just to 
    % make it more challenging)

    % Define parameter upper/lower bounds
    lb = zeros(1,D);
    ub = ones(1,D);
    plb = 0.1*ones(1,D);
    pub = 0.9*ones(1,D);
    noise = [];

    if transform_to_unconstrained_coordinates
        trinfo = warpvars(D,lb,ub,plb,pub);     % Transform to unconstrained space
        trinfo.mu = zeros(1,D);     % Necessary for retro-compatibility
        trinfo.delta = ones(1,D);
    else
        trinfo = [];
    end
    
    % Mode    
    if alpha > 1 && beta > 1
        xmin = (alpha - 1) / (alpha + beta - 2);
    elseif alpha <= 1 && beta > 1
        xmin = 0;
    elseif alpha > 1 && beta <= 1
        xmin = 1;
    else
        xmin = 0.5;
    end
    fval = ll_beta(xmin,data);
    xmin = xmin*ones(1,D);
    fval = fval*ones(1,D);
    xmin_post = xmin;
    fval_post = fval;
    
    xmin = warpvars(xmin,'d',trinfo);
    fval = fval + warpvars(xmin,'logp',trinfo);

    Mean = (alpha/(alpha + beta))*ones(1,D);
    Var = alpha*beta / ((alpha + beta)^2*(alpha + beta + 1));
    Cov = diag(Var*ones(1,D));
    Mode = xmin;

    y.D = D;
    y.LB = warpvars(lb,'d',trinfo);
    y.UB = warpvars(ub,'d',trinfo);
    y.PLB = warpvars(plb,'d',trinfo);
    y.PUB = warpvars(pub,'d',trinfo);

    y.lnZ = 0;              % Log normalization factor
    y.Mean = Mean;          % Distribution moments
    y.Cov = Cov;
    y.Mode = Mode;          % Mode of the pdf
    y.ModeFval = fval;

    priorMean = 0.5*(y.PUB + y.PLB);
    priorSigma2 = (0.5*(y.PUB - y.PLB)).^2;
    priorCov = diag(priorSigma2);
    y.Prior.Mean = priorMean;
    y.Prior.Cov = priorCov;

    % Compute each coordinate separately
    xmin_post = warpvars(xmin_post,'d',trinfo);
    fval_post = fval_post + warpvars(xmin_post,'logp',trinfo);

    y.Post.Mean = Mean;
    y.Post.Mode = xmin_post;          % Mode of the posterior
    y.Post.ModeFval = fval_post;        
    y.Post.lnZ = 0;
    y.Post.Cov = Cov;

    % Compute marginals
    nx = 2^13;
    xmesh = linspace(eps, 1-eps, nx);
    beta_pdf = exp((alpha-1)*log(xmesh) + (beta-1)*log1p(-xmesh));
    beta_pdf = beta_pdf/(qtrapz(beta_pdf)*(xmesh(2)-xmesh(1))); % Ensure normalization
    for i = 1:D
        y.Post.MarginalBounds(:,i) = [0; 1];
        y.Post.MarginalPdf(i,:) = beta_pdf;    
    end
    
    % Save data and coordinate transformation struct
    data.trinfo = trinfo;
    y.Data = data;
    y.DeterministicFlag = true;
    
else
    
    % Iteration call -- evaluate objective function
    
    % Transform unconstrained variables to original space
    x_orig = warpvars(x,'i',infprob.Data.trinfo);
    dy = warpvars(x,'logpdf',infprob.Data.trinfo);   % Jacobian correction
    
    % Compute log likelihood of data and possibly std of log likelihood
    LL = ll_beta(x_orig,infprob.Data);
    y_std = 0;
    y = LL + dy;
    
end

end

%--------------------------------------------------------------------------
function ll = ll_beta(theta,data)

D = numel(theta);
alpha = data.beta_params(1);
beta = data.beta_params(2);

ll = sum((alpha-1)*log(theta) + (beta-1)*log1p(-theta));
ll = ll - D*(gammaln(alpha) + gammaln(beta) - gammaln(alpha + beta));

end