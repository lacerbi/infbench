function optimState = setupvars_wsabi(x0,LB,UB,PLB,PUB,optimState,options,prnt)
%SETUPVARS Initialize and transform variables for WSABI+.

nvars = size(PLB,2);

if isempty(LB); LB = -Inf(1,nvars); end
if isempty(UB); UB = Inf(1,nvars); end
if isscalar(LB); LB = LB*ones(1,nvars); end
if isscalar(UB); UB = UB*ones(1,nvars); end
if isscalar(x0); x0 = x0*ones(1,nvars); end

% Starting point
if isempty(x0) || any(~isfinite(x0))   % Invalid/not provided starting point
    if prnt > 0
        fprintf('Initial starting point is invalid or not provided. Starting from center of plausible region.\n');
    end
    x0 = 0.5*(PLB + PUB);       % Midpoint
end

optimState.LB_orig = LB;
optimState.UB_orig = UB;
optimState.PLB_orig = PLB;
optimState.PUB_orig = PUB;
optimState.LBeps_orig = LB + (UB - LB)*options.TolBoundX;
optimState.UBeps_orig = UB - (UB - LB)*options.TolBoundX;


%% Transform variables
trinfo = warpvars_wsabi(nvars,LB,UB,PLB,PUB);
trinfo.x0_orig = x0;
if ~isfield(trinfo,'R_mat'); trinfo.R_mat = []; end
if ~isfield(trinfo,'scale'); trinfo.scale = []; end

optimState.LB = warpvars_wsabi(LB,'dir',trinfo);
optimState.UB = warpvars_wsabi(UB,'dir',trinfo);
optimState.PLB = warpvars_wsabi(PLB,'dir',trinfo);
optimState.PUB = warpvars_wsabi(PUB,'dir',trinfo);
optimState.x0 = warpvars_wsabi(x0,'dir',trinfo);

optimState.trinfo = trinfo;


%% Setup Gaussian prior

if ~isempty(options.PriorMean)
    if isempty(options.PriorCov)
        error('Need to specify both a prior mean and prior covariance.');
    end
    if any(isfinite(optimState.LB_orig)) || any(isfinite(optimState.UB_orig))
        error('Cannot specify Gaussian prior with finite bounds.');
    end
    
    priorMu = options.PriorMean;
    if isscalar(priorMu); priorMu = priorMu*ones(1,nvars); end
    
    priorCov = options.PriorCov;
    if isscalar(priorCov); priorCov = priorCov*ones(1,nvars); end
    
    if size(priorCov,1) == size(priorCov,2)
        priorCov = diag(priorCov)';
    end
    
    optimState.priorMu = priorMu;
    optimState.priorVar = priorCov;
    optimState.removePrior = false;
    
else
    
    optimState.priorMu = 0.5*(optimState.PLB + optimState.PUB);
    optimState.priorVar = (0.5*(optimState.PUB - optimState.PLB)).^2;
    optimState.removePrior = true;
    
end

optimState.priorLogNorm = 0.5*sum(log(2*pi*optimState.priorVar));

optimState.LB_search = optimState.PLB - 6*sqrt(optimState.priorVar);
optimState.UB_search = optimState.PUB + 6*sqrt(optimState.priorVar);


%% Get warnings state

optimState.DefaultWarnings.singularMatrix = warning('query','MATLAB:singularMatrix');
warning('off',optimState.DefaultWarnings.singularMatrix.identifier);