function y = infbench_funnel(x,infprob)
%INFBENCH_FUNNEL Inference benchmark log pdf -- Neal's funnel density.

if isempty(x)
    D = infprob(1);         % Call with the number of dimensions

    sigma_prior = 3;
    sigma = 3;
        
    y.func = ['@(x,infprob) ' mfilename '(x,infprob)'];
    y.D = D;
    y.LB = -Inf(1,D);
    y.UB = Inf(1,D);
    y.PLB = -sigma_prior*ones(1,D);
    y.PUB = sigma_prior*ones(1,D);
        
    
    y.sigma2 = sigma^2;
        
    range = 5*(y.PUB-y.PLB);

    priorMean = 0.5*(y.PUB + y.PLB);
    priorCov = diag((0.5*(y.PUB - y.PLB)).^2);
    y.Prior.Mean = priorMean;
    y.Prior.Cov = priorCov;

    % Unused
    y.lnZ = 0;
    y.Mean = priorMean;
    y.Cov = priorCov;
    y.Mode = priorMean;
        
    % Compute normalization constant
    fun = @(x) f0(x,sigma,sigma_prior,D);
    Z = integral(fun,y.PLB(1)-range(1),y.PUB(1)+range(1));
    y.Post.lnZ = log(Z);
    
    y.Post.Mean = zeros(1,D);
    
    % Compute posterior variance (covariance is zero)
    post_var = zeros(1,D);
    
    fun = @(x) f2(x,sigma,sigma_prior,D);
    post_var(1) = integral(fun,y.PLB(1)-range(1),y.PUB(1)+range(1))/Z;
    
    fun = @(x) f2b(x,sigma,sigma_prior,D);
    post_var(2:D) = integral(fun,y.PLB(1)-range(1),y.PUB(1)+range(1))/Z;
    
    y.Post.Cov = diag(post_var);
    
    % Compute posterior mode analytically
    tau2 = sigma_prior^2*sigma^2/(sigma^2 + sigma_prior^2);
%     fun = @(x) 0.5*x(:,1).^2/y.sigma2 ...
%         + 0.5*sum(x(:,2:D).^2,2)./exp(x(:,1)) + 0.5*(D-1)*x(:,1) ...
%         + 0.5*sum(x.^2,2)/sigma_prior^2;
%     Mode = fminunc(fun,zeros(1,D));   % Numerical check
    y.Post.Mode = [-0.5*(D-1)*tau2,zeros(1,D-1)];
        
else
    sigma2 = infprob.sigma2;
    D = infprob.D;
    y = -0.5*x(:,1).^2/sigma2 -0.5*log(sigma2) ...
        - 0.5*sum(x(:,2:D).^2,2)./exp(x(:,1)) - 0.5*(D-1)*x(:,1) ...
        - 0.5*D*log(2*pi);
end

end

%--------------------------------------------------------------------------
function y = f0(x,sigma,sigma_prior,D)

tau = sigma_prior*sigma/sqrt(sigma^2 + sigma_prior^2);
y = normpdf(x,0,tau) ./ sqrt(2*pi*(sigma_prior^2 + exp(x))).^(D-1) ...
    ./ sqrt(2*pi*(sigma^2 + sigma_prior^2));

end

%--------------------------------------------------------------------------
function y = f2(x,sigma,sigma_prior,D)

tau = sigma_prior*sigma/sqrt(sigma^2 + sigma_prior^2);
y = x.^2 .* ...
    normpdf(x,0,tau) ./ sqrt(2*pi*(sigma_prior^2 + exp(x))).^(D-1) ...
    ./ sqrt(2*pi*(sigma^2 + sigma_prior^2));

end

function y = f2b(x,sigma,sigma_prior,D)

tau = sigma_prior*sigma/sqrt(sigma^2 + sigma_prior^2);
y = (sigma_prior^2.*exp(x))./(sigma_prior^2 + exp(x)) .* ...
    normpdf(x,0,tau) ./ sqrt(2*pi*(sigma_prior^2 + exp(x))).^(D-1) ...
    ./ sqrt(2*pi*(sigma^2 + sigma_prior^2));

end
