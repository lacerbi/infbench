function [y,y_std] = infbench_acerbi2012(x,infprob,mcmc_params)
%INFBENCH_ACERBI2012 Inference benchmark log pdf -- sensorimotor timing.

if nargin < 3; mcmc_params = []; end

problem_name = 'acerbi2012';
infbench_fun = str2func(['infbench_' problem_name]);

if isempty(x)
    if isempty(infprob) % Generate this document        
        fprintf('\n');

        % Add sampler directory to MATLAB path
        pathstr = fileparts(mfilename('fullpath'));
        addpath([pathstr,filesep(),'..',filesep(),'infalgo',filesep(),'parallel-GP-SL',filesep(),'utils',filesep(),'mcmcstat-master']);
                
        for n = 1     
            name = ['S' num2str(n)];
            
            infprob = infbench_fun([],n);
            infprob.DeterministicFlag = true;
            if isempty(mcmc_params); id = 0; else; id = mcmc_params(1); end

            D = infprob.D;
            trinfo = infprob.Data.trinfo;
            
            % Prior used for sampling
            infprob.PriorMean = infprob.Prior.Mean;
            infprob.PriorVar = diag(infprob.Prior.Cov)';
            infprob.PriorVolume = prod(infprob.UB - infprob.LB);
            infprob.PriorType = 'uniform';            

            LB = infprob.LB;
            UB = infprob.UB;
            PLB = infprob.PLB;
            PUB = infprob.PUB;
            
            if id == 0
            
                % Compute optimum
                Nopts = 5;
                
                opts = struct('Display','iter','MaxFunEvals',2e3);

                for iOpt = 1:Nopts                
                    x0 = rand(1,D).*(PUB - PLB) + PLB;
                    fun = @(x) -infbench_fun(x,infprob);
                    [xnew(iOpt,:),fvalnew(iOpt)] = bads(fun,x0,LB,UB,PLB,PUB,[],opts);
%                    [xnew(iOpt,:),fvalnew(iOpt)] = fmincon(fun,x0,[],[],[],[],LB,UB,[],opts);
                end
                
                [fvalnew,idx_best] = min(fvalnew);
                xnew = xnew(idx_best,:);
                
                fvalnew = -fvalnew;
                xmin = warpvars(xnew,'inv',trinfo);
                fval = fvalnew + warpvars(xnew,'logp',trinfo);

                x0 = xmin;
                x0 = warpvars(x0,'d',trinfo);   % Convert to unconstrained coordinates            
                fun = @(x) -logpost(x,infprob);
                [xnew,fvalnew] = bads(fun,x0,LB,UB,PLB,PUB,[],opts);
%                [xnew,fvalnew] = fmincon(fun,x0,[],[],[],[],LB,UB,[],opts);

                fvalnew = -fvalnew;
                xmin_post = xmin;
                xmin_post = warpvars(xnew,'inv',trinfo);
                fval_post = fvalnew + warpvars(xnew,'logp',trinfo);

                fprintf('\t\t\tcase %d\n',n);
                fprintf('\t\t\t\tname = ''%s'';\n\t\t\t\txmin = %s;\n\t\t\t\tfval = %s;\n',name,mat2str(xmin),mat2str(fval));
                fprintf('\t\t\t\txmin_post = %s;\n\t\t\t\tfval_post = %s;\n',mat2str(xmin_post),mat2str(fval_post));
                
            elseif id > 0 && n == mcmc_params(2)

                rng(id);
                widths = 0.5*(PUB - PLB);
                logpfun = @(x) logpost(x,infprob);
                
                % Number of samples
                if numel(mcmc_params) > 2
                    Ns = mcmc_params(3);
                else
                    Ns = 1e3;
                end
                
                W = 2*(infprob.D+1);    % Number of walkers
                
                sampleopts.Thin = 7;
                sampleopts.Burnin = Ns*sampleopts.Thin;
                sampleopts.Display = 'notify';
                sampleopts.Diagnostics = false;
                sampleopts.VarTransform = false;
                sampleopts.InversionSample = false;
                sampleopts.FitGMM = false;
                sampleopts.TolX = 1e-5;
                % sampleopts.TransitionOperators = {'transSliceSampleRD'};

                x0 = infprob.Post.Mode;
                x0 = warpvars(x0,'d',trinfo);
                
                [Xs,lls,exitflag,output] = eissample_lite(logpfun,x0,Ns,W,widths,LB,UB,sampleopts);
                
%                 % Re-evaluate final lls with higher precision
%                 fprintf('Re-evaluate log posteriors...\n');
%                 infprob.Ng = 501;
%                 lls = zeros(size(Xs,1),1);
%                 for i = 1:size(Xs,1)
%                     lls(i) = logpost(Xs(i,:),infprob);
%                     if mod(i,ceil(size(Xs,1)/10))==0; fprintf('%d/%d..',i,size(Xs,1)); end
%                 end
%                 fprintf('\nDone.\n');
                
                
                filename = [problem_name '_mcmc_n' num2str(n) '_id' num2str(id) '.mat'];
                save(filename,'Xs','lls','exitflag','output');                
            end
            
        end
        
    else
        % Initialization call -- define problem and set up data
        n = infprob(1);
        
        % Are we transforming the entire problem to unconstrained space?
        transform_to_unconstrained_coordinates = false;
        
        % Add problem directory to MATLAB path
        pathstr = fileparts(mfilename('fullpath'));
        addpath([pathstr,filesep(),problem_name]);
                
        D = 5;
        xmin = NaN(1,D);       fval = Inf;
        xmin_post = NaN(1,D);  fval_post = Inf;
        Mean_mcmc = NaN(1,D);  Cov_mcmc = NaN(D,D);    lnZ_mcmc = NaN;
        
        if transform_to_unconstrained_coordinates
            trinfo = warpvars(D,lb,ub,plb,pub);     % Transform to unconstrained space
            trinfo.mu = zeros(1,D);     % Necessary for retro-compatibility
            trinfo.delta = ones(1,D);
        else
            trinfo = [];
        end
        
        switch n
			case {1,101}
                nid = 1;
                subjid = 2;
				xmin = [0.0903024950646795 0.0308776679681614 0.70145987368378 0.13736980468675 0.0100000000005821];
				fval = -3839.17327065442;
				xmin_post = [0.0903022675098327 0.0308786661764316 0.701457336964017 0.137370632316879 0.0100000000000091];
				fval_post = -3836.24763972821;                
                % R_max = 1.009. Ntot = 100000. Neff_min = 97466.3. Total funccount = 3683564.
                Mean_mcmc = [0.0889507384236801 0.031037033631714 0.700860186274828 0.136268613484507 0.0112742636946608];
                Cov_mcmc = [5.69319324140584e-05 -7.21868242423165e-05 -5.71600545036279e-05 8.73572534343407e-05 -1.98445964866227e-07;-7.21868242423165e-05 0.000114752676423207 9.04981392015061e-05 -0.000122248013438833 1.71305359428043e-07;-5.71600545036279e-05 9.04981392015061e-05 0.000189266266554489 -0.000133295925616371 2.42345926662173e-07;8.73572534343407e-05 -0.000122248013438833 -0.000133295925616371 0.00016674960610084 -3.35907235005097e-07;-1.98445964866227e-07 1.71305359428043e-07 2.42345926662173e-07 -3.35907235005097e-07 1.51752670977801e-06];
                lnZ_mcmc = -3859.86804439948;                
        end

        % Read and preprocess data
        temp = load('acerbi2012internal.mat');
        data = timer_loaddata(temp.acerbi2012time.exp3{subjid});
        
        % Bin response data in 20 ms bins
        data.binsize = 0.02;
        dr = data.binsize;
        data.X(:,2) = round((data.X(:,2) - 0.5*dr)/dr)*dr + 0.5*dr;
        data.R = data.X(:,2);
        
        S = data.S;
        dS = S(2)-S(1);
        
        % Define parameter upper/lower bounds
        lb = [0.01,0.01,min(S)/2,(dS/2),0.01];
        ub = [0.5,0.5,2*max(S),(2*(S(end)-S(1))),0.2];
        plb = [0.05,0.05,min(S),(dS),0.02];
        pub = [0.25,0.25,max(S),(S(end)-S(1)),0.05];
        noise = [];        
                
        xmin = warpvars(xmin,'d',trinfo);
        fval = fval + warpvars(xmin,'logp',trinfo);
        
        Mean = zeros(1,D);
        Cov = eye(D);
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
                
        y.Post.Mean = Mean_mcmc;
        y.Post.Mode = xmin_post;          % Mode of the posterior
        y.Post.ModeFval = fval_post;        
        y.Post.lnZ = lnZ_mcmc;
        y.Post.Cov = Cov_mcmc;
                
        if n > 100
            data.IBSNreps = 0; % Deterministic problems            
        else
            data.IBSNreps = 150; % Reps used for IBS estimator
        end
        
        % Read marginals from file
        marginals = load([problem_name '_marginals.mat']);
        y.Post.MarginalBounds = marginals.MarginalBounds{nid};
        y.Post.MarginalPdf = marginals.MarginalPdf{nid};
        
        % Save data and coordinate transformation struct
        data.trinfo = trinfo;
        y.Data = data;
        y.DeterministicFlag = (data.IBSNreps == 0);
                
    end
    
else
    
    % Iteration call -- evaluate objective function
    
    % Transform unconstrained variables to original space
    x_orig = warpvars(x,'i',infprob.Data.trinfo);
    dy = warpvars(x,'logpdf',infprob.Data.trinfo);   % Jacobian correction
    
    % Compute log likelihood of data and possibly std of log likelihood
    if infprob.DeterministicFlag
        LL = timer_loglike(x_orig,infprob.Data);
        y_std = 0;
    else        
        Nibs = infprob.Data.IBSNreps/5;
        IBSNreps = 5; % Split to avoid memory errors (should fix ibslike)
        ibs_opts = struct('Nreps',IBSNreps,...
            'ReturnPositive',true,'ReturnStd',false);        
        for iter = 1:Nibs
            [LL(iter),y_var(iter)] = ibslike(@timer_gendata,x_orig,infprob.Data.R,[],ibs_opts,infprob.Data);
        end
        LL = mean(LL);
        y_std = sqrt(mean(y_var)/Nibs);
    end
    y = LL + dy;
    
end

end

%--------------------------------------------------------------------------
function y = logpost(x,infprob)    
    y = infbench_acerbi2012(x,infprob);
    lnp = infbench_lnprior(x,infprob);
    y = y + lnp;
end

function timer_test(infprob)

theta = rand(1,numel(infprob.PLB)).*(infprob.PUB - infprob.PLB) + infprob.PLB;
LL = timer_loglike(theta,infprob.Data);
ibs_opts = struct('Nreps',infprob.Data.IBSNreps,...
    'ReturnPositive',true,'ReturnStd',true);
[LL_ibs,y_std] = ibslike(@timer_gendata,theta,infprob.Data.R,[],ibs_opts,infprob.Data);

[LL,LL_ibs,y_std]
[LL - LL_ibs, y_std]

end

function data = timer_loaddata(dataset)

data.X = [];
for i = 1:numel(dataset.X)
    data.X = [data.X; dataset.X{i}];    
end

% Remove missed trials
idx_nan = any(isnan(data.X),2);
data.X(idx_nan,:) = [];

data.S = unique(dataset.xrange);
data.prior = dataset.prior;

end

function ll = timer_loglike(theta,data)

MAXSD = 5;
Ns = 101;
Nx = 401;

ws = theta(1);
wm = theta(2);
mu_prior = theta(3);
sigma_prior = (theta(4));
if numel(theta) < 5
    lambda = 0.01; % Lapse rate
else
    lambda = theta(5);
end

if isfield(data,'binsize') && ~isempty(data.binsize)
    dr = data.binsize;
else
    dr = 0;
end

srange = linspace(0,2,Ns)';
ds = srange(2) - srange(1);

ll = zeros(size(data.X,1),1);
Nstim = numel(data.S);

for iStim = 1:Nstim
    mu_s = data.S(iStim);
    sigma_s = ws*mu_s;
    xrange = linspace(max(0,mu_s-MAXSD*sigma_s), mu_s+MAXSD*sigma_s, Nx);
    dx = xrange(2)-xrange(1);
    
    xpdf = normpdf(xrange,mu_s,sigma_s);
    xpdf = xpdf / (qtrapz(xpdf)*dx);
    
    like = bsxfun_normpdf(xrange,srange,ws*srange + eps);
    prior = normpdf(srange,mu_prior,sigma_prior);
    
    post = bsxfun(@times,like,prior);
    post = bsxfun(@rdivide,post,qtrapz(post,1)*ds);
    
    post_mean = qtrapz(bsxfun(@times,post,srange),1)*ds;
    s_hat = post_mean / (1 + wm^2);
    
    idx = data.X(:,3) == iStim;
    
    sigma_m = wm*s_hat;
    if dr > 0
        pr = bsxfun_normcdf(data.R(idx)+0.5*dr,s_hat,sigma_m) ...
            - bsxfun_normcdf(data.R(idx)-0.5*dr,s_hat,sigma_m);
    else 
        pr = bsxfun_normpdf(data.R(idx),s_hat,sigma_m);
    end
    ll(idx) = qtrapz(bsxfun(@times,xpdf,pr),2)*dx;
end

if dr > 0
    ll = log(ll*(1-lambda) + lambda/((srange(end)-srange(1))/dr));    
else
    ll = log(ll*(1-lambda) + lambda/(srange(end)-srange(1)));
end

ll = sum(ll);

end


function R = timer_gendata(theta,idx,data)

MAXSD = 5;
Nx = 401;
Ns = 101;

ws = theta(1);
wm = theta(2);
mu_prior = theta(3);
sigma_prior = (theta(4));
if numel(theta) < 5
    lambda = 0.01; % Lapse rate
else
    lambda = theta(5);
end

if isfield(data,'binsize') && ~isempty(data.binsize)
    dr = data.binsize;
else
    dr = 0;
end

S_vec = data.X(idx,3);
srange = linspace(0,2,Ns)';
ds = srange(2) - srange(1);

R = zeros(numel(idx),1);
Nstim = numel(data.S);

% Compute s_hat for a given x

sigma_s_max = ws*max(data.S);
xrange = linspace(max(0,min(data.S)-MAXSD*sigma_s_max), max(data.S)+MAXSD*sigma_s_max, Nx);
like = bsxfun_normpdf(xrange,srange,ws*srange + eps);
prior = normpdf(srange,mu_prior,sigma_prior);
post = bsxfun(@times,like,prior);
post = bsxfun(@rdivide,post,qtrapz(post,1)*ds);
post_mean = qtrapz(bsxfun(@times,post,srange),1)*ds;
s_hat = post_mean / (1 + wm^2);

for iStim = 1:Nstim
    mu_s = data.S(iStim);
    S_idx = S_vec == iStim;
    Nt = sum(S_idx);
    if Nt == 0; continue; end
        
    % Generate noisy sensory measurements
    sigma_s = ws*mu_s;
    xx = max(0,mu_s + sigma_s*randn(Nt,1));
            
%     like = bsxfun_normpdf(xx,srange,ws*srange + eps);
%     prior = normpdf(srange,mu_prior,sigma_prior);
%     
%     post = bsxfun(@times,like,prior);
%     post = bsxfun(@rdivide,post,qtrapz(post,1)*ds);    
%     post_mean = qtrapz(bsxfun(@times,post,srange),1)'*ds;
    s_hat_x = lininterp1(xrange,s_hat,xx);    
    
    sigma_m = wm*s_hat_x;
    
    % Add motor noise
    rr = max(0,s_hat_x + sigma_m.*randn(Nt,1));
    
    if dr > 0
        rr = round((rr - 0.5*dr)/dr)*dr + 0.5*dr;
    end
    R(S_idx) = rr;
end

lapse_idx = rand(size(R,1),1) < lambda;

if sum(lapse_idx) > 0
    if dr > 0
        Nr_bins = (srange(end)-srange(1))/dr;
        R_lapse = randi(Nr_bins,[sum(lapse_idx),1])*dr - 0.5*dr + srange(1);
    else
        L = srange(end)-srange(1);
        R_lapse = rand([sum(lapse_idx),1])*L + srange(1);
    end
    R(lapse_idx) = R_lapse;
end

end

function y = bsxfun_normpdf(x,mu,sigma)
%BSXFUN_NORMPDF Vectorized normal probability density function (pdf).

if isscalar(mu)
    y = bsxfun(@rdivide, exp(-0.5*bsxfun(@rdivide, x - mu, sigma).^2), sigma)/sqrt(2*pi);
elseif isscalar(sigma)
    y = exp(-0.5*(bsxfun(@minus, x, mu)/sigma).^2)/(sigma*sqrt(2*pi));
else
    y = bsxfun(@rdivide, exp(-0.5*bsxfun(@rdivide, bsxfun(@minus, x, mu), sigma).^2), sigma)/sqrt(2*pi);
end

end


function p = bsxfun_normcdf(x,mu,sigma)
%BSXFUN_NORMCDF Vectorized normal cumulative distribution function (cdf).

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

end

function Vout = lininterp1(X,V,Xq)

if isvector(X); X = X(:); end
if isvector(V); V = V(:); end

% Gets the x spacing
dx = 1./(X(2,:)-X(1,:));  % one over to perform divide only once
Xq = bsxfun(@minus,Xq,X(1,:));      % subtract minimum of x

Xqi = bsxfun(@times,Xq,dx)+1;          % indices of nearest-lower-neighbors
flag = Xqi<1 | Xqi>size(X,1) | isnan(Xqi);
                                % finds indices out of bounds
Xqi = floor(Xqi);
                                
V = [V; NaN(1, size(V,2))];
Xqi(flag) = size(V,1)-1;

delta = Xqi - bsxfun(@times,Xq,dx);

linind1 = bsxfun(@plus, Xqi, size(V,1)*(0:size(V,2)-1));
linind2 = bsxfun(@plus, Xqi + 1, size(V,1)*(0:size(V,2)-1));

out1 = reshape(V(linind1),size(Xq));
out2 = reshape(V(linind2),size(Xq));

out1(isnan(out1)) = 0;
out2(isnan(out2)) = 0;

Vout = bsxfun(@times,delta,out1) + bsxfun(@times,1-delta,out2);
Vout(flag) = NaN;

end