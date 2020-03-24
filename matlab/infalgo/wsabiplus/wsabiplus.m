function [ log_mu, log_Var, clktime, xxIter, loglIter, hyp ] = wsabiplus( ...
            loglikfun,      ... 1) Handle to log-likelihood function
            priorMu,        ... 2) Gaussian prior mean, 1 x D
            priorVar,       ... 3) Gaussian prior covariance, D x D
            LB,             ... 4) Lower bound for search / uniform prior
            UB,             ... 5) Upper bound for search / uniform prior
            x0,             ... 6) Use as starting point
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

D = numel(priorMu);

% Default options
defopts.Method          = 'L';              % Method ('L' or 'M')
defopts.MaxFunEvals     = (D+2)*50;         % Maximum number of fcn evals
defopts.AlphaFactor     = 0.8;              % Alpha offset fraction, as in paper
defopts.Display         = 'iter';           % Print output
defopts.GPInputLengths  = [];               % Vector of GP input lengths
defopts.GPOutputLength  = 1;                % GP ouput scale
defopts.GPHypVar        = 100^2;            % Variance of prior over GP hyperparams
defopts.Nsearch         = 2^13;             % Starting search points for acquisition fcn
defopts.AcqFun          = [];               % Acquisition function
defopts.Plot            = 0;                % Make diagnostic plots at each iteration
defopts.SpecifyTargetNoise = false;         % Loglik function returns noise estimate (SD) as second output

for f = fields(defopts)'
    if ~isfield(options,f{:}) || isempty(options.(f{:}))
        options.(f{:}) = defopts.(f{:});
    end
end

method = upper(options.Method);
if method(1) ~= 'L' && method(1) ~= 'M'
    error('wsabi:UnknownMethod', ...
        'Allowed values for METHOD are (L)inearized and (M)oment-matching.');
end

numSamples = options.MaxFunEvals+1;
alpha = options.AlphaFactor;
if ischar(options.Display)
    switch lower(options.Display)
        case {'yes','on','iter','notify','all'}; printing = 1;
        case {'no','off','none'}; printing = 0;
    end
else
    printing = options.Display;
end
Nsearch = options.Nsearch;

ss = sqrt(diag(priorVar))';
if isempty(LB); LB = priorMu - 6*ss; end
if isempty(UB); UB = priorMu + 6*ss; end
range = [LB; UB];

if isempty(options.GPInputLengths)
    options.GPInputLengths = (0.5*(range(2,:)-range(1,:))/sqrt(3));
end
lambda = options.GPOutputLength;
hypVar = options.GPHypVar;


% Relabel prior mean and covariance for brevity of code.
bb          = priorMu;
BB          = diag(priorVar)';

% Relabel input hyperparameters squared for brevity of code.
VV          = (options.GPInputLengths(:)').^2;

jitterNoise = 1e-6;     % Jitter on the GP model over the log likelihood.
numEigs     = inf;      % If trying to use Nystrom ** NOT RECOMMENDED **
hypOptEvery = 1;        % 1 => optimise hyperparameters every iteration.

dim         = length(bb);   % Dimensionality of integral.

% Limit absolute range of likelihood model hyperparameters for stability.
hypLims     = 30*ones(1,dim+1); 

% Allocate Storage
mu              = zeros(numSamples-1,1);
logscaling      = zeros(numSamples-1,1);
Var             = zeros(numSamples-1,1);
clktime         = zeros(numSamples-1,1);
lHatD_0_tmp     = zeros(numSamples,1);
loglHatD_0_tmp  = zeros(size(lHatD_0_tmp));
loglsdhat_tmp   = zeros(size(lHatD_0_tmp));
hyp             = zeros(1,1+dim);

% Variance of prior over GP hyperparameters
if isscalar(hypVar); hypVar = hypVar*ones(size(hyp)); end

% Minimiser options (fmincon for hyperparameters)
options1                        = optimset('fmincon');
options1.Display                = 'none';
options1.GradObj                = 'off';
options1.Algorithm              = 'active-set';
options1.TolX                   = 1e-5;
options1.TolFun                 = 1e-5;
options1.MaxTime                = 5;
options1.MaxFunEvals            = 100*dim;
%options1.UseParallel           = 'always';
options1.AlwaysHonorConstraints = 'true';

% Minimiser options (fmincon if desired for active sampling)
options2                        = optimset('fmincon');
options2.Display                = 'none';
options2.GradObj                = 'on';
%options2.DerivativeCheck       = 'on';
options2.TolX                   = 1e-5;
options2.TolFun                 = 1e-5;
%options2.MaxTime               = 0.5;
%options2.MaxFunEvals           = 75;
options2.UseParallel            = 'always';
options2.AlwaysHonorConstraints = 'true';

% Minimiser options (CMAES - advised for active sampling)
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

% Initial Sample:
xx = zeros(numSamples,dim);
if isempty(x0)
    xx(end,:) = bb+1e-6;    % Prior mean
else
    xx(end,:) = x0+1e-6;
end
currNumSamples = 1;

for t = 1:numSamples - 1
    if printing
        if ~mod(t,10)
            fprintf('Iter %d. Log Current Mean Integral: %g +/- %g.\n', ...
                            t, log(mu(t-1)) + logscaling(t-1), ...
                            sqrt(Var(t-1))/mu(t-1));
        end
        
        if printing == 2 && ~mod(t,10)
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
    % Pre-process new samples -- i.e. convert to log space etc.
    
    % Get batch of samples & variables from stack.
    tmpT            = cputime; 
    xxIter          = xx(numSamples-currNumSamples+1:end,:); % Curr samples
    
    % Call loglik handle for latest sample.
    if options.SpecifyTargetNoise
        [loglHatD_0_tmp(numSamples-currNumSamples+1), ...
            loglsdhat_tmp(numSamples-currNumSamples+1)] ...
            = loglikfun( xxIter(1,:) );        
    else
        loglHatD_0_tmp(numSamples-currNumSamples+1) = loglikfun( xxIter(1,:) );
    end
                                           
    % Find the max in log space.                                       
    logscaling(t)   = max(loglHatD_0_tmp(numSamples-currNumSamples+1:end));
    
    % Scale batch by max, and exponentiate.
    lHatD_0_tmp(numSamples-currNumSamples+1:end) = ... 
      exp(loglHatD_0_tmp(numSamples-currNumSamples+1:end) - logscaling(t));    
  
    % Evaluate the offset, alpha fraction of minimum value seen.
    aa      = alpha * min( lHatD_0_tmp(numSamples-currNumSamples+1:end) );
    
    % Transform into sqrt space.
    lHatD = sqrt(abs(lHatD_0_tmp(numSamples-currNumSamples+1:end)- aa)*2);
         
    if options.SpecifyTargetNoise
        sigma_tmp = loglsdhat_tmp(numSamples-currNumSamples+1:end);
        ucb_tmp = sqrt(abs(exp(loglHatD_0_tmp(numSamples-currNumSamples+1:end) + sigma_tmp - logscaling(t)) - aa)*2);
        lcb_tmp = sqrt(max(0,exp(loglHatD_0_tmp(numSamples-currNumSamples+1:end) - sigma_tmp - logscaling(t)) - aa)*2);
        s2hat = (0.5*(ucb_tmp - lcb_tmp)).^2;
    else
        s2hat = 1e-5;       
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ML-II On GP Likelihood model hyperparameters
    
    hyp(1)          = log(lambda);
    hyp(2:end)      = log(VV);
    
    if currNumSamples > 3 &&  ~mod(currNumSamples,hypOptEvery)
        hypLogLik = @(x) logLikGPDim(xxIter, lHatD, s2hat, x, hypVar);
        [hyp,nll] = fmincon(hypLogLik, ...
                         hyp,[],[],[],[],-hypLims,hypLims,[],options1);

        % Attempt restart of hyperparameters
        if ~mod(currNumSamples,hypOptEvery*10)
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
    xxIterScaled    = xxIter .* repmat(sqrt(1./VV),currNumSamples,1);
    
    % Squared distance matrix
    dist2           = pdist2_squared_fast(xxIterScaled, xxIterScaled);
    
    % Evaluate Gram matrix
    Kxx = lambda.^2 * (1/(prod(2*pi*VV).^0.5)) * exp(-0.5*dist2);
%     Kxx = Kxx + ...
%               lambda.^2*(1/(prod(2*pi*VV).^0.5))*jitterNoise*eye(size(Kxx));
    if isscalar(s2hat)
        Kxx = Kxx + s2hat*eye(size(Kxx));
    else
        Kxx = Kxx + diag(s2hat);
    end
    Kxx = Kxx/2 + Kxx'/2; % Make sure symmetric for stability.
    
    % Invert Gram matrix.
    invKxx = Kxx \ eye(size(Kxx));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Expected value of integral:
    ww              = invKxx * lHatD;
    
    Yvar            = (VV.*VV + 2*VV.*BB)./VV; 
    postProd        = ww*ww';
    
    xx2sq           = xxIter .* repmat(sqrt(1./Yvar),currNumSamples,1);
    bbch            = bb .* sqrt(1./Yvar);
    
    xx2sqFin        = pdist2_squared_fast(xx2sq,bbch);
    xxIterScaled2   = xxIter .* repmat(sqrt(BB./(VV.*Yvar)),currNumSamples,1);
    
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
    
    % Variance of the integral, before scaling back up:
    
    %---------------- Tmp Vars to calculate first term -------------------%
        
    GG_coeff = lambda^6 * 1/prod(16*pi^4*(VV.*VV + 2*VV.*BB).*VV.*BB.* ...
               (((VV.*VV + 2*VV.*BB)./BB) + ...
               ((VV.*VV + 2*VV.*BB)./VV) + VV + BB))^0.5 * ...
               prod(2*pi*(VV.*VV + 2*VV.*BB))^0.5;
    
    tmpVar   = ((VV.*VV + 2*VV.*BB).*(Yvar + Yvar.*(VV./BB + VV + BB)));
    
    xx3sq    = xxIter .* ...
               repmat(sqrt((Yvar.*VV)./tmpVar),currNumSamples,1);
    bb3ch    = bb .* sqrt((Yvar.*VV)./tmpVar);
    
    xx2sqGGlh = xx2sqFin + pdist2_squared_fast(xx3sq,bb3ch);
    
    xx4sq    = xxIter .* ...
             repmat(sqrt((VV.*Yvar)./tmpVar),currNumSamples,1);
    bb4ch    = bb .* sqrt((VV.*Yvar)./tmpVar);
    
    xx5sq    = xxIter .* ...
             repmat(sqrt((VV.*Yvar.^2)./(tmpVar.*BB)),currNumSamples,1);
    bb5ch    = bb .* sqrt((VV.*Yvar.^2)./(tmpVar.*BB));
    
    xx2sqGGrh       = pdist2_squared_fast(xx4sq,bb4ch) + ...
                      pdist2_squared_fast(xx5sq,bb5ch);
                  
    xxIterScaled3   = xxIter .* ...
                      repmat(sqrt(Yvar.*BB./tmpVar),currNumSamples,1);
                  
    dist4           = pdist2_squared_fast(xxIterScaled3,xxIterScaled3);
    xxIterScaled4   = xxIter .* ...
                      repmat(sqrt(BB./(Yvar.*VV)),currNumSamples,1);
                  
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
        xxIterScaledVar2_0 = xxIter .* repmat(sqrt(ScaledVar2_0),currNumSamples,1);
        bbScaledVar2_0     = bb .* sqrt(ScaledVar2_0);
        distVar2_0         = pdist2_squared_fast(xxIterScaledVar2_0, bbScaledVar2_0);     

        ScaledVar2_1       = ((VV.*VV + 2*VV.*BB)./(tmp_2_mainvar));
        xxIterScaledVar2_1 = xxIter .* repmat(sqrt(ScaledVar2_1),currNumSamples,1);
        bbScaledVar2_1     = bb .* sqrt(ScaledVar2_1);
        distVar2_1         = pdist2_squared_fast(xxIterScaledVar2_1, bbScaledVar2_1);

        ScaledVar2_2       = ((VV.*VV + 2*VV.*BB).*(VV.*VV + 2*VV.*BB))./(tmp_2_mainvar.*(VV.*BB));
        xxIterScaledVar2_2 = xxIter .* repmat(sqrt(ScaledVar2_2),currNumSamples,1);
        bbScaledVar2_2     = bb .* sqrt(ScaledVar2_2);
        distVar2_2         = pdist2_squared_fast(xxIterScaledVar2_2, bbScaledVar2_2);

        ScaledVar2_3       = ScaledVar2_1;
        xxIterScaledVar2_3 = xxIter .* repmat(sqrt(ScaledVar2_3),currNumSamples,1);
        bbScaledVar2_3     = bb .* sqrt(ScaledVar2_3);
        distVar2_3         = pdist2_squared_fast(xxIterScaledVar2_3, bbScaledVar2_3);

        ScaledVar2_4       = ((VV.*BB)./(tmp_2_mainvar));
        xxIterScaledVar2_4 = xxIter .* repmat(sqrt(ScaledVar2_4),currNumSamples,1);
        bbScaledVar2_4     = bb .* sqrt(ScaledVar2_4);
        distVar2_4         = pdist2_squared_fast(xxIterScaledVar2_4, bbScaledVar2_4);

        distVar2lh         = distVar2_1 + distVar2_2;
        distVar2rh         = distVar2_0 + distVar2_3 + distVar2_4;
        distVar2           = pdist2(distVar2lh,-distVar2rh);  

        %---------------- Tmp Vars to calculate third term -------------------%

        tmp_3_mainvar      = VV.*VV + 2*VV.*BB;
        tmp_3_coeff        = lambda^12 * 1/prod(4*pi^2*tmp_3_mainvar);

        ScaledVar3_0       = (VV./((VV.*VV+2*VV.*BB)));
        xxIterScaledVar3_0 = xxIter .* repmat(sqrt(ScaledVar3_0),currNumSamples,1);
        bbScaledVar3_0     = bb .* sqrt(ScaledVar3_0);
        distVar3_0         = pdist2_squared_fast(xxIterScaledVar3_0, bbScaledVar3_0);    

        ScaledVar3_1       = (VV./(VV.*VV+2*VV.*BB));
        xxIterScaledVar3_1 = xxIter .* repmat(sqrt(ScaledVar3_1),currNumSamples,1);
        bbScaledVar3_1     = bb .* sqrt(ScaledVar3_1);
        distVar3_1         = pdist2_squared_fast(xxIterScaledVar3_1, bbScaledVar3_1);

        ScaledVar3_2       = (BB./(VV.*VV+2*VV.*BB));
        xxIterScaledVar3_2 = xxIter .* repmat(sqrt(ScaledVar3_2),currNumSamples,1);
        distVar3_2         = pdist2_squared_fast(xxIterScaledVar3_2, xxIterScaledVar3_2);

        ScaledVar3_3       = (BB./(VV.*VV+2*VV.*BB));
        xxIterScaledVar3_3 = xxIter .* repmat(sqrt(ScaledVar3_3),currNumSamples,1);
        distVar3_3         = pdist2_squared_fast(xxIterScaledVar3_3, xxIterScaledVar3_3);

        ScaledVar3_4       = (VV./(VV.*VV+2*VV.*BB));
        xxIterScaledVar3_4 = xxIter .* repmat(sqrt(ScaledVar3_4),currNumSamples,1);
        bbScaledVar3_4     = bb .* sqrt(ScaledVar3_4);
        distVar3_4         = pdist2_squared_fast(xxIterScaledVar3_4, bbScaledVar3_4);

        ScaledVar3_5       = (VV./(VV.*VV+2*VV.*BB));
        xxIterScaledVar3_5 = xxIter .* repmat(sqrt(ScaledVar3_5),currNumSamples,1);
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
    % Actively select next sample point:
     
    % Define acquisition function
    acqfun_name = str2func(['expectedVar' method]);
    acqfun = @(x) acqfun_name( transp(x), lambda, VV, ...
        lHatD, xxIter, invKxx, jitterNoise, bb, BB, aa );
        
    if Nsearch > 0
        % Perform shotgun search

        ll = loglHatD_0_tmp(numSamples-currNumSamples+1:end) ...
            - 0.5*sum(bsxfun(@rdivide,bsxfun(@minus,xxIter,bb).^2,BB),2) ...
            - 0.5*sum(log(BB));
        
        if size(xxIter,1) > max(4,dim)
            [cma_Sigma,cma_mu] = covcma(xxIter,ll,[],'descend');
        else
            cma_Sigma = []; cma_mu = [];            
        end        
        
        % 1) Draw a small number of points from prior
        Nrnd = ceil(Nsearch/10);
        murnd = bb;
        sigmarnd = sqrt(BB);
        Xsearch = bsxfun(@plus,bsxfun(@times,randn(Nrnd,dim),sigmarnd),murnd);
                
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
                    bsxfun(@times,range(2,:)-range(1,:),rand(Nrnd,dim)))];                    
            end
        end
        
        % Evaluate acquisition function on all candidate search points
        aval = acqfun(Xsearch);
        
        % Take best point
        [strtFval,idx] = min(aval);
        strtSamp = Xsearch(idx,:);
    
        insigma = 0.1*sqrt(diag(cma_Sigma));
    else
        if rand < 1.1 % Sample starting location for search from prior.
            strtSamp = mvnrnd(bb,diag(BB),1);
        else
            strtSamp = 2*range(2,:).*rand(1,dim) - 50;
        end
        strtFval = acqfun(strtSamp);
        insigma = [];
    end
    
    % Using global optimiser (cmaes):
    % insigma = exp((log(VV) + log(BB))/4);
    % if size(xxIter,1) > dim; insigma = std(xxIter); end
    [newX,cmaesFval] = cmaes_modded( ['expectedVar' method(1)], strtSamp', insigma(:), opts, lambda, VV, ...
                  lHatD, xxIter, invKxx, jitterNoise, bb, BB, aa);
    newX = newX';
    
    % If CMA-ES somehow does not improve from starting point, just use that
    if strtFval < cmaesFval; newX = strtSamp; end
    
    xx(numSamples-currNumSamples,:) = newX;
    
    clktime(t) = cputime - tmpT;
    
    currNumSamples = currNumSamples + 1;

end

fprintf('\n done. \n');
log_mu  = log(mu) + logscaling;
log_Var = log(Var) + 2*logscaling;
loglIter = loglHatD_0_tmp(numSamples-currNumSamples+2:end,:);

end
