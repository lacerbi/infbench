function [log_mu,log_Var,output,samples] = wsabiplus( ...
            loglikfun,      ... 1) Handle to log-likelihood function
            PLB,            ... 2) Plausible lower bound
            PUB,            ... 3) Plausible upper bound
            LB,             ... 4) Hard lower bound
            UB,             ... 5) Hard upper bound
            x0,             ... 6) Starting point (optional)
            options         ...
            )
        
% Output structures:
% log_mu:   log of the integral posterior mean.
% log_var:  log of the integral posterior variance.
% clktime:  vector of times per iteration, may want to cumulative sum.
% xxIter:   numSamples x D array of sample locations used to build model.
% loglIter: numSamples x 1 vector of log likelihood fcn evaluations.
% hyp:      integral hyperparameters.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4; LB = []; end
if nargin < 5; UB = []; end
if nargin < 6; x0 = []; end
if nargin < 7; options = []; end

D = numel(PLB);

% Default options
defopts.Method          = 'L';              % Method ('L' or 'M')
defopts.MaxFunEvals     = (D+2)*50;         % Maximum number of fcn evals
defopts.AlphaFactor     = 0.8;              % Alpha offset fraction, as in paper
defopts.Display         = 'iter';           % Print output
defopts.GPOutputLength  = 1;                % GP ouput scale
defopts.GPInputLengths  = [];               % Vector of GP input lengths
defopts.GPOutputHypVar  = log(10)^2;        % Variance of prior over GP output length
defopts.GPInputHypVar   = log(20)^2;        % Variance of prior over GP input lengths (normalized)
defopts.Nsearch         = 2^10;             % Starting search points for acquisition fcn
defopts.AcqFun          = [];               % Acquisition function
defopts.Plot            = 0;                % Make diagnostic plots at each iteration
defopts.SpecifyTargetNoise = false;         % Loglik function returns noise estimate (SD) as second output
defopts.TolBoundX       = 1e-5;             % Tolerance on closeness to bound constraints (fraction of total range)
defopts.PriorMean       = [];               % Gaussian prior mean
defopts.PriorCov        = [];               % Gaussian prior covariance
defopts.LCBFactor       = 0;                % Lower Confidence Bound parameter
defopts.PreciseSearch   = true;             % Assume zero noise during active search?
defopts.BoundedTransform = 'logit';         % Input transform for bounded variables
defopts.NMCMCsamples    = 5e3;              % Posterior samples (via MCMC)
defopts.SavePosterior   = Inf;              % Array of # evaluations when to save posteriors (Inf = at the end)

for f = fields(defopts)'
    if ~isfield(options,f{:}) || isempty(options.(f{:}))
        options.(f{:}) = defopts.(f{:});
    end
end


add2path(); % Add folders to path

method = upper(options.Method);
if method(1) ~= 'L' && method(1) ~= 'M'
    error('wsabi:UnknownMethod', ...
        'Allowed values for METHOD are (L)inearized and (M)oment-matching.');
end

numSamples = options.MaxFunEvals;
alpha = options.AlphaFactor;
if ischar(options.Display)
    switch lower(options.Display)
        case {'yes','on','iter','notify','all'}; prnt = 1;
        case {'no','off','none'}; prnt = 0;
    end
else
    prnt = options.Display;
end
Nsearch = options.Nsearch;

% Initialize optimization structure
optimState = setupvars_wsabi(x0,LB,UB,PLB,PUB,[],options,prnt);
range = [optimState.LB_search; optimState.UB_search];
prange = optimState.PUB-optimState.PLB;
if isempty(options.GPInputLengths)
    options.GPInputLengths = 0.5*prange/sqrt(3);
end
lambda = options.GPOutputLength;

% Prior over GP hyperparameters    

hypMu = [0,log(0.1*prange)];
hypVar = [options.GPOutputHypVar,options.GPInputHypVar*ones(1,D)];


% Relabel prior mean and covariance for brevity of code
bb          = optimState.priorMu;
BB          = optimState.priorVar;

% Relabel input hyperparameters squared for brevity of code
VV          = (options.GPInputLengths(:)').^2;

jitterNoise = 1e-6;     % Jitter on the GP model over the log likelihood
hypOptEvery = 1;        % 1 => optimise hyperparameters every iteration

% Limit absolute range of likelihood model hyperparameters for stability
hypLims     = 30*ones(1,D+1); 

% Allocate Storage
mu              = zeros(numSamples-1,1);
logscaling      = zeros(numSamples-1,1);
Var             = zeros(numSamples-1,1);
clktime         = zeros(numSamples-1,1);
loglHatD_0_tmp  = zeros(numSamples,1);
if options.SpecifyTargetNoise
    loglsdhat_tmp   = zeros(numSamples,1);
end
hyp             = zeros(1,1+D);

% Minimiser options (fmincon for hyperparameters)
options1                        = optimset('fmincon');
options1.Display                = 'none';
options1.GradObj                = 'off';
options1.Algorithm              = 'active-set';
options1.TolX                   = 1e-5;
options1.TolFun                 = 1e-5;
options1.MaxTime                = 5;
options1.MaxFunEvals            = 100*D;
%options1.UseParallel           = 'always';
options1.AlwaysHonorConstraints = 'true';

% Minimiser options (CMAES for active sampling)
opts                            = cmaes_modded('defaults');
opts.LBounds                    = range(1,:)';
opts.UBounds                    = range(2,:)';
opts.DispModulo                 = Inf;
opts.DispFinal                  = 'off';
opts.SaveVariables              = 'off';
opts.LogModulo                  = 0;
opts.CMA.active                 = 1;      % Use Active CMA (generally better)
opts.EvalParallel               = 'on';
%opts.PopSize                   = 100;
%opts.Restarts                  = 1;

% MCMC sampling from posterior
NMCMCsamples = options.NMCMCsamples;
mcmc_flag = (nargout > 3) && NMCMCsamples > 0;
if mcmc_flag
    mcmc_iters = options.SavePosterior;
    mcmc_iters(mcmc_iters == Inf) = numSamples;
    mcmc_iters = unique(mcmc_iters);
end

% Initial Sample
optimState.X = zeros(numSamples,D);
optimState.X(1,:) = optimState.x0;

for t = 1:numSamples
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Pre-process new samples -- i.e. convert to log space etc.
    
    % Get batch of samples & variables from stack
    tmpT            = cputime; 
    xxIter          = optimState.X(1:t,:);
        
    % Convert back to original space
    xx = xxIter(end,:);
    xx_orig = warpvars_wsabi(xx,'i',optimState.trinfo);    
    
    % Call loglik handle for latest sample
    if options.SpecifyTargetNoise
        [fval_orig,loglsdhat_tmp(t)] = loglikfun(xx_orig);
    else
        fval_orig = loglikfun(xx_orig);
    end
    
    fval = fval_orig + warpvars_wsabi(xx,'logp',optimState.trinfo);
    if optimState.removePrior
        fval = fval + 0.5*sum(((xx - bb).^2)./BB,2) + optimState.priorLogNorm;
    end
    loglHatD_0_tmp(t) = fval;
    
    % Rescaling factor
    if options.SpecifyTargetNoise
        % Find max LCB
        logscaling(t)   = max(loglHatD_0_tmp(1:t) - options.LCBFactor*loglsdhat_tmp(1:t));
    else
        % Simply find the max in log space                             
        logscaling(t)   = max(loglHatD_0_tmp(1:t));        
    end
    
    % Scale batch by max, and exponentiate
    lHatD_0_tmp = exp(loglHatD_0_tmp(1:t) - logscaling(t));    
    
    % Evaluate the offset, alpha fraction of minimum value seen
    aa      = alpha * min(lHatD_0_tmp);
    
    % Transform into sqrt space
    lHatD = sqrt(2*abs(lHatD_0_tmp - aa));
    
    if options.SpecifyTargetNoise
        % Noise warping via unscented transform
        u = 0.6745; % norminv(0.75)
        sigma_tmp = loglsdhat_tmp(1:t);
        ucb_tmp = sqrt(abs(exp(loglHatD_0_tmp(1:t) + u*sigma_tmp - logscaling(t)) - aa)*2);
        lcb_tmp = sqrt(max(0,exp(loglHatD_0_tmp(1:t) - u*sigma_tmp - logscaling(t)) - aa)*2);
        s2hat = (0.5*(ucb_tmp - lcb_tmp)/u).^2;
        s2hat = max(s2hat, max(s2hat)*jitterNoise);
    else
        s2hat = jitterNoise;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% ML-II On GP Likelihood model hyperparameters
    
    hyp(1)          = log(lambda);
    hyp(2:end)      = log(VV);
        
    if t > 3 && mod(t,hypOptEvery) == 0
        hypLogLik = @(x) logLikGPDim(xxIter,lHatD,s2hat,x,hypMu,hypVar);
        [hyp,nll] = fmincon(hypLogLik, ...
                         hyp,[],[],[],[],-hypLims,hypLims,[],options1);

        % Attempt restart of hyperparameters
        if mod(t,hypOptEvery*10) == 0
            hyp0 = randn(100,numel(hyp));
            nll0 = zeros(size(hyp0,1),1);
            for iHyp = 1:size(hyp0,1); nll0(iHyp) = hypLogLik(hyp0(iHyp,:)); end
            [~,idx] = min(nll0);
            [hyp2,nll2] = fmincon(hypLogLik, ...
                             hyp0(idx,:),[],[],[],[],-hypLims,hypLims,[],options1);
            if nll2 < nll; hyp = hyp2; end
        end
    end
    
    lambda          = exp(hyp(1));
    VV              = exp(hyp(2:end));
    
    % Scale samples by input length scales.
    xxIterScaled    = xxIter .* repmat(sqrt(1./VV),t,1);
    
    % Squared distance matrix
    dist2           = pdist2_squared_fast(xxIterScaled, xxIterScaled);
    
    % Evaluate Gram matrix
    Kxx = lambda.^2 * (1/(prod(2*pi*VV).^0.5)) * exp(-0.5*dist2);
    if isscalar(s2hat)
        Kxx = Kxx + s2hat*eye(size(Kxx));
    else
        Kxx = Kxx + diag(s2hat);
    end
    Kxx = Kxx/2 + Kxx'/2; % Make sure symmetric for stability.
    
    % Invert Gram matrix.
    invKxx = Kxx \ eye(size(Kxx));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %% Expected value of integral
    
    ww              = invKxx * lHatD;
    
    Yvar            = (VV.*VV + 2*VV.*BB)./VV; 
    postProd        = ww*ww';
    
    xx2sq           = xxIter .* repmat(sqrt(1./Yvar),t,1);
    bbch            = bb .* sqrt(1./Yvar);
    
    xx2sqFin        = pdist2_squared_fast(xx2sq,bbch);
    xxIterScaled2   = xxIter .* repmat(sqrt(BB./(VV.*Yvar)),t,1);
    
    dist4           = pdist2_squared_fast(xxIterScaled2,xxIterScaled2);
    
    if method(1) == 'M'
        % Sigma^2 term:
        sig2t = - ... 
                lambda^4 * (1 / prod(4*pi^2*((VV.*VV + 2*VV.*BB)))^0.5) * ...
                 exp(-0.5 * (pdist2(xx2sqFin,-xx2sqFin) + dist4)) .* invKxx;
    else
        sig2t = 0;
    end
    
    YY              = lambda^4 * ... 
                    (1 / prod(4*pi^2*((VV.*VV + 2*VV.*BB)))^0.5) * ...
                    exp(-0.5 * (pdist2(xx2sqFin,-xx2sqFin) + dist4)) .* ...
                    postProd + sig2t;

    if 0
        xx1 = xxIter ./ sqrt(2*VV);    
        xx2 = 0.5*(xxIter - bb) ./ sqrt(0.5*VV + BB);
        YY2 = lambda^4 * ...
            (1 / prod(4*pi^2*((VV.*VV + 2*VV.*BB)))^0.5) * ...
            exp(-0.5 * pdist2_squared_fast(xx1,xx1) - 0.5*pdist2_squared_fast(xx2,-xx2)) .* ...
            postProd;
        mu2 = aa + 0.5*sum(YY2(:));
    end
    
    % Mean of the integral at iteration 't', before scaling back up:
    mu(t) = aa + 0.5*sum(YY(:));
    if method(1) == 'M'
        mu(t) = mu(t) + 0.5*lambda.^2 * (1/(prod(2*pi*VV)^0.5));
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %% Variance of the integral, before scaling back up
    
    %---------------- Tmp Vars to calculate first term -------------------%
        
    GG_coeff = lambda^6 * 1/prod(16*pi^4*(VV.*VV + 2*VV.*BB).*VV.*BB.* ...
               (((VV.*VV + 2*VV.*BB)./BB) + ...
               ((VV.*VV + 2*VV.*BB)./VV) + VV + BB))^0.5 * ...
               prod(2*pi*(VV.*VV + 2*VV.*BB))^0.5;
    
    tmpVar   = ((VV.*VV + 2*VV.*BB).*(Yvar + Yvar.*(VV./BB + VV + BB)));
    
    xx3sq    = xxIter .* ...
               repmat(sqrt((Yvar.*VV)./tmpVar),t,1);
    bb3ch    = bb .* sqrt((Yvar.*VV)./tmpVar);
    
    xx2sqGGlh = xx2sqFin + pdist2_squared_fast(xx3sq,bb3ch);
    
    xx4sq    = xxIter .* ...
             repmat(sqrt((VV.*Yvar)./tmpVar),t,1);
    bb4ch    = bb .* sqrt((VV.*Yvar)./tmpVar);
    
    xx5sq    = xxIter .* ...
             repmat(sqrt((VV.*Yvar.^2)./(tmpVar.*BB)),t,1);
    bb5ch    = bb .* sqrt((VV.*Yvar.^2)./(tmpVar.*BB));
    
    xx2sqGGrh       = pdist2_squared_fast(xx4sq,bb4ch) + ...
                      pdist2_squared_fast(xx5sq,bb5ch);
                  
    xxIterScaled3   = xxIter .* ...
                      repmat(sqrt(Yvar.*BB./tmpVar),t,1);
                  
    dist4           = pdist2_squared_fast(xxIterScaled3,xxIterScaled3);
    xxIterScaled4   = xxIter .* ...
                      repmat(sqrt(BB./(Yvar.*VV)),t,1);
                  
    dist5           = pdist2_squared_fast(xxIterScaled4,xxIterScaled4);
    
    
    GG      = GG_coeff * postProd .* ...
            exp(-0.5*(pdist2(xx2sqGGlh, -xx2sqGGrh) + dist4));
     
    YY2     = lambda^4 * (1 / prod(4*pi^2*((VV.*VV + 2*VV.*BB)))^0.5) * ...
            exp(-0.5 * (pdist2(xx2sqFin,-xx2sqFin) + dist5)) .* ...
            repmat(ww',length(ww),1);
    
    Var(t)  = (sum(sum(GG)) - sum(YY2,2)'*(invKxx * sum(YY2,2)));
    
    if method(1) == 'M'
        %---------------- Tmp Vars to calculate second term ------------------%

        tmp_2_mainvar = (BB.*(VV.*VV + 2*VV.*BB) + (((VV.*VV)./BB) + 2*VV) +...
                        (((VV.*VV+2*VV.*BB).*(VV.*VV+2*VV.*BB))./(VV)) + ...
                        VV.*(VV.*VV + 2*VV.*BB));
        tmp_2_coeff   = lambda^8 * 1/prod(8*pi^3*(VV.*VV+2*VV.*BB))^0.5 * ...
                        prod((VV+2*BB).*(VV.*VV./BB+2*VV))^0.5 * ...
                        1/prod(tmp_2_mainvar)^0.5;

        ScaledVar2_0       = ((VV)./(VV.*VV+2*VV.*BB));
        xxIterScaledVar2_0 = xxIter .* repmat(sqrt(ScaledVar2_0),t,1);
        bbScaledVar2_0     = bb .* sqrt(ScaledVar2_0);
        distVar2_0         = pdist2_squared_fast(xxIterScaledVar2_0, bbScaledVar2_0);     

        ScaledVar2_1       = ((VV.*VV + 2*VV.*BB)./(tmp_2_mainvar));
        xxIterScaledVar2_1 = xxIter .* repmat(sqrt(ScaledVar2_1),t,1);
        bbScaledVar2_1     = bb .* sqrt(ScaledVar2_1);
        distVar2_1         = pdist2_squared_fast(xxIterScaledVar2_1, bbScaledVar2_1);

        ScaledVar2_2       = ((VV.*VV + 2*VV.*BB).*(VV.*VV + 2*VV.*BB))./(tmp_2_mainvar.*(VV.*BB));
        xxIterScaledVar2_2 = xxIter .* repmat(sqrt(ScaledVar2_2),t,1);
        bbScaledVar2_2     = bb .* sqrt(ScaledVar2_2);
        distVar2_2         = pdist2_squared_fast(xxIterScaledVar2_2, bbScaledVar2_2);

        ScaledVar2_3       = ScaledVar2_1;
        xxIterScaledVar2_3 = xxIter .* repmat(sqrt(ScaledVar2_3),t,1);
        bbScaledVar2_3     = bb .* sqrt(ScaledVar2_3);
        distVar2_3         = pdist2_squared_fast(xxIterScaledVar2_3, bbScaledVar2_3);

        ScaledVar2_4       = ((VV.*BB)./(tmp_2_mainvar));
        xxIterScaledVar2_4 = xxIter .* repmat(sqrt(ScaledVar2_4),t,1);
        bbScaledVar2_4     = bb .* sqrt(ScaledVar2_4);
        distVar2_4         = pdist2_squared_fast(xxIterScaledVar2_4, bbScaledVar2_4);

        distVar2lh         = distVar2_1 + distVar2_2;
        distVar2rh         = distVar2_0 + distVar2_3 + distVar2_4;
        distVar2           = pdist2(distVar2lh,-distVar2rh);  

        %---------------- Tmp Vars to calculate third term -------------------%

        tmp_3_mainvar      = VV.*VV + 2*VV.*BB;
        tmp_3_coeff        = lambda^12 * 1/prod(4*pi^2*tmp_3_mainvar);

        ScaledVar3_0       = (VV./((VV.*VV+2*VV.*BB)));
        xxIterScaledVar3_0 = xxIter .* repmat(sqrt(ScaledVar3_0),t,1);
        bbScaledVar3_0     = bb .* sqrt(ScaledVar3_0);
        distVar3_0         = pdist2_squared_fast(xxIterScaledVar3_0, bbScaledVar3_0);    

        ScaledVar3_1       = (VV./(VV.*VV+2*VV.*BB));
        xxIterScaledVar3_1 = xxIter .* repmat(sqrt(ScaledVar3_1),t,1);
        bbScaledVar3_1     = bb .* sqrt(ScaledVar3_1);
        distVar3_1         = pdist2_squared_fast(xxIterScaledVar3_1, bbScaledVar3_1);

        ScaledVar3_2       = (BB./(VV.*VV+2*VV.*BB));
        xxIterScaledVar3_2 = xxIter .* repmat(sqrt(ScaledVar3_2),t,1);
        distVar3_2         = pdist2_squared_fast(xxIterScaledVar3_2, xxIterScaledVar3_2);

        ScaledVar3_3       = (BB./(VV.*VV+2*VV.*BB));
        xxIterScaledVar3_3 = xxIter .* repmat(sqrt(ScaledVar3_3),t,1);
        distVar3_3         = pdist2_squared_fast(xxIterScaledVar3_3, xxIterScaledVar3_3);

        ScaledVar3_4       = (VV./(VV.*VV+2*VV.*BB));
        xxIterScaledVar3_4 = xxIter .* repmat(sqrt(ScaledVar3_4),t,1);
        bbScaledVar3_4     = bb .* sqrt(ScaledVar3_4);
        distVar3_4         = pdist2_squared_fast(xxIterScaledVar3_4, bbScaledVar3_4);

        ScaledVar3_5       = (VV./(VV.*VV+2*VV.*BB));
        xxIterScaledVar3_5 = xxIter .* repmat(sqrt(ScaledVar3_5),t,1);
        bbScaledVar3_5     = bb .* sqrt(ScaledVar3_5);
        distVar3_5         = pdist2_squared_fast(xxIterScaledVar3_5, bbScaledVar3_5);

        %----------------- Combine terms to get total var --------------------%

        tmp_1 = lambda.^4/prod(8*pi^2*VV.*(0.5*(VV+2*BB)+BB))^0.5;

        tmp_2 = invKxx .* (tmp_2_coeff * exp(-0.5*distVar2));

        tmp_3 = tmp_3_coeff * exp(-0.5*distVar3_0)' * ...
                (invKxx.*exp(-0.5*distVar3_2))* exp(-0.5*distVar3_1) * ...
                exp(-0.5*distVar3_4)'*(invKxx.*exp(-0.5*distVar3_3)) * ...
                exp(-0.5*distVar3_5);

        Var(t) = Var(t) + 0.5*(tmp_1 - 2*sum(tmp_2(:)) + sum(tmp_3(:)));  
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %% Printout every few iterations
    if prnt
        if ~mod(t,10)
            fprintf('Iter %d. Log Current Mean Integral: %g +/- %g.\n', ...
                            t, log(mu(t-1)) + logscaling(t-1), ...
                            sqrt(Var(t-1))/mu(t-1));
        end
        
        if prnt == 2 && ~mod(t,10)
            gp.method = method(1);
            gp.lambda = lambda;
            gp.VV = VV;
            gp.noise = jitterNoise;
            gp.alpha = aa;
            gp.X = xxIter;
            gp.y = lHatD;
            gp.invKxx = invKxx;
            sqrtgp_plot(gp,[],100);
            drawnow;
            % pause;
        end        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %% If requested, sample from GP
    if mcmc_flag && any(t == mcmc_iters)
    
        % Create GP struct
        gp.method = method(1);
        gp.lambda = lambda;
        gp.VV = VV;
        gp.noise = jitterNoise;
        gp.alpha = aa;
        gp.X = xxIter;
        gp.y = lHatD;
        gp.invKxx = invKxx;
        
        logprior = @(x_) -0.5*sum(bsxfun(@rdivide,bsxfun(@minus,x_, bb).^2,BB),2) - optimState.priorLogNorm;

        logpostfun = @(x_) log(sqrtgp_pred(gp,x_)) + logprior(x_);
        %ss = stratais(gp.X,gp.y,ifun,2e4);            
        %ss = warpvars_wsabi(ss,'inv',optimState.trinfo);
        %cornerplot(ss);

        sampleopts.Burnin = ceil(NMCMCsamples/2);
        sampleopts.Thin = 1;
        sampleopts.Display = 'off';
        sampleopts.Diagnostics = false;
        sampleopts.VarTransform = false;
        sampleopts.InversionSample = false;
        sampleopts.FitGMM = false;

        widths = std(gp.X,[],1);
        W = 2*(D+1);
        % Take starting points from high posterior density region
        hpd_frac = 0.2;
        N = numel(gp.y);
        N_hpd = min(N,max(W,round(hpd_frac*N)));
        yy = logpostfun(gp.X);        
        [yy,ord] = sort(yy,'descend');                
        X_hpd = gp.X(ord(1:N_hpd),:);
        
        % Take only points with non-zero likelihoods
        yy = yy(1:N_hpd);
        X_hpd = X_hpd(isfinite(yy), :);                
        x0 = X_hpd(randperm(size(X_hpd,1),min(W,size(X_hpd,1))),:);

        try
            ss = eissample_lite(logpostfun,x0,NMCMCsamples,W,widths,-Inf,Inf,sampleopts);
        catch
            ss = NaN(NMCMCsamples,D);
        end
        ss = warpvars_wsabi(ss,'inv',optimState.trinfo);
        
        if numel(mcmc_iters) == 1
            samples = ss;            
        else
            idx = find(t == unique(mcmc_iters),1);
            samples{idx} = ss;
        end
        
        % cornerplot(ss);
    end    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %% Actively select next sample point
    if t < numSamples

        xxIterScaled = bsxfun(@rdivide,xxIter,sqrt(VV));
        
        % Define acquisition function
        acqfun_name = str2func(['expectedVar' method]);
        if options.PreciseSearch
            s2hat_search = 0;
        else
            s2hat_search = s2hat;
        end
        acqfun = @(x) acqfun_name( x', s2hat_search, lambda, VV, ...
            lHatD, xxIterScaled, invKxx, jitterNoise, bb, BB, aa );

        if Nsearch > 0
            % Perform shotgun search

            ll = loglHatD_0_tmp(1:t) ...
                - 0.5*sum(bsxfun(@rdivide,bsxfun(@minus,xxIter,bb).^2,BB),2) ...
                - 0.5*sum(log(BB));

            if size(xxIter,1) > max(4,D)
                [cma_Sigma,cma_mu] = covcma(xxIter,ll,[],'descend');
            else
                cma_Sigma = []; cma_mu = [];            
            end        

            % 1) Draw a small number of points from prior
            Nrnd = ceil(Nsearch/10);
            murnd = bb;
            sigmarnd = sqrt(BB);
            Xsearch = bsxfun(@plus,bsxfun(@times,randn(Nrnd,D),sigmarnd),murnd);

            % 2) Draw points around current points
            Nrnd = ceil((Nsearch - size(Xsearch,1))/2);
            if Nrnd > 0 && t > 3
                cma_SS = 0.01*covcma(xxIter,ll,xxIter(1,:),'descend');
                sqlik = exp(0.5*ll-0.5*max(ll));            
                weights = sqlik./sum(sqlik);
                iirnd = randsample((1:size(xxIter,1))', Nrnd, true, weights);
                murnd = xxIter(iirnd,:);

                Xsearch = [Xsearch; ...
                    mvnrnd(murnd,cma_SS,Nrnd)];
            end

            % 3) Draw remaining points
            Nrnd = Nsearch - size(Xsearch,1);
            if Nrnd > 0
                if ~isempty(cma_Sigma)  % Draw from multivariate normal ~ hpd region            
                    Xsearch = [Xsearch; ...
                        mvnrnd(cma_mu,cma_Sigma,Nrnd)];
                else            % Uniform draw inside search box
                    Xsearch = [Xsearch; ...
                        bsxfun(@plus, range(1,:), ...
                        bsxfun(@times,range(2,:)-range(1,:),rand(Nrnd,D)))];                    
                end
            end

            % Evaluate acquisition function on all candidate search points
            aval = acqfun(Xsearch);

            % Take best point
            [strtFval,idx] = min(aval);
            strtSamp = Xsearch(idx,:);
            
            % Fix starting INSIGMA
            if ~isempty(cma_Sigma)
                diag_sigma = sqrt(diag(cma_Sigma));
                diag_sigma = min(diag_sigma,10*sqrt(BB(:)));
                diag_sigma = max(diag_sigma,max(diag_sigma)/1e5);
                insigma = 0.1*diag_sigma;
            else
                insigma = [];
            end
        else
            if rand < 1.1 % Sample starting location for search from prior.
                strtSamp = mvnrnd(bb,diag(BB),1);
            else
                strtSamp = 2*range(2,:).*rand(1,D) - 50;
            end
            strtFval = acqfun(strtSamp);
            insigma = [];
        end

        % Using global optimiser (cmaes):
        % insigma = exp((log(VV) + log(BB))/4);
        % if size(xxIter,1) > dim; insigma = std(xxIter); end
        
        [newX,cmaesFval] = cmaes_modded( ['expectedVar' method(1)], strtSamp', insigma(:), opts, s2hat, lambda, VV, ...
                      lHatD, xxIterScaled, invKxx, jitterNoise, bb, BB, aa);
        newX = newX';

        % If CMA-ES somehow does not improve from starting point, just use that
        if strtFval < cmaesFval; newX = strtSamp; end

        optimState.X(t+1,:) = newX;
    end
    
    clktime(t) = cputime - tmpT;

end

fprintf('\n done. \n');
log_mu  = log(mu) + logscaling;
log_Var = log(Var) + 2*logscaling;

% Output structure
if nargout > 2

    % Transformed space
    output.X_transformed = xxIter;
    output.fval_transformed = loglHatD_0_tmp(1:t);
    if optimState.removePrior
        output.fval_transformed = output.fval_transformed - 0.5*sum(((xxIter - bb).^2)./BB,2) - optimState.priorLogNorm;
    end    
    
    % Transform back to original space
    loglIter = loglHatD_0_tmp(1:t);
    loglIter = loglIter - warpvars_wsabi(xxIter,'logp',optimState.trinfo);
    if optimState.removePrior
        loglIter = loglIter - 0.5*sum(((xxIter - bb).^2)./BB,2) - optimState.priorLogNorm;
    end
    xxIter = warpvars_wsabi(xxIter,'inv',optimState.trinfo);
    if options.SpecifyTargetNoise
        lsdIter = loglsdhat_tmp(1:t);
    else
        lsdIter = [];
    end
    
    output.X = xxIter;
    output.fval = loglIter;
    output.fval_sd = lsdIter;
    output.clktime = clktime;    
    output.gp_hyp = hyp;
    output.trinfo = optimState.trinfo;    
    if mcmc_flag; output.mcmc_iters = mcmc_iters; end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function add2path()
%ADD2PATH Adds WSABI+ subfolders to MATLAB path.

subfolders = {'utils'};
pathCell = regexp(path, pathsep, 'split');
baseFolder = fileparts(mfilename('fullpath'));

onPath = true;
for iFolder = 1:numel(subfolders)
    folder = [baseFolder,filesep,subfolders{iFolder}];    
    if ispc  % Windows is not case-sensitive
      onPath = onPath & any(strcmpi(folder, pathCell));
    else
      onPath = onPath & any(strcmp(folder, pathCell));
    end
end

% ADDPATH is slow, call it only if folders are not on path
if ~onPath
    addpath(genpath(fileparts(mfilename('fullpath'))));
end

end
