function [mu,Sigma,lnZ,output] = inflaplace(fun,x0,PLB,PUB,LB,UB,options)
%INFLAPLACE Approximate Bayesian inference via Laplace approximation.

% Default options
defopts.Display             = 'iter';       % Print output
defopts.BoundedTransform    = 'logit';      % Input transform for bounded variables
defopts.Optimization        = 'fminunc';    % Optimization algorithm

for f = fields(defopts)'
    if ~isfield(options,f{:}) || isempty(options.(f{:}))
        options.(f{:}) = defopts.(f{:});
    end
end

add2path(); % Add folders to path

if ischar(options.Display)
    switch lower(options.Display)
        case {'yes','on','iter','notify','all'}; prnt = 1;
        case {'no','off','none'}; prnt = 0;
    end
else
    prnt = options.Display;
end

% Initialize optimization structure
optimState = setupvars_laplace(x0,LB,UB,PLB,PUB,[],options,prnt);

% Define function in transformed space
fun_warped = @(x) fun(warpvars_laplace(x,'i',optimState.trinfo)) + warpvars_laplace(x,'logp',optimState.trinfo);

% Optimization
D = size(x0,2);
x0 = optimState.x0;

% Optimize from starting point to find MAP estimate
switch lower(options.Optimization)
    case 'fminunc'
        mu = fminunc(@(x) -fun_warped(x),x0);
    case 'bads'
        mu = bads(@(x) -fun_warped(x),x0,optimState.LB,optimState.UB,optimState.PLB,optimState.PUB);
    otherwise
        error('Unknown optimization algorithm.');
end

% Hessian matrix
A = -hessian(fun_warped,mu);
[C,p] = chol(A);
nugget = 1e-6;
while p > 0
    if prnt
        fprintf('Cholesky decomposition failed. Retry with nugget = %.2g.\n', nugget);
    end
    [C,p] = chol(A + nugget*eye(D));
    nugget = nugget*1.2;
end
lnZ = fun_warped(mu) + 0.5*D*log(2*pi) - sum(log(diag(C)));
Sigma = C'\(C\eye(D));

% OUTPUT struct
output.trinfo = optimState.trinfo;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function add2path()
%ADD2PATH Adds LAPLACE subfolders to MATLAB path.

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function optimState = setupvars_laplace(x0,LB,UB,PLB,PUB,optimState,options,prnt)
%SETUPVARS Initialize and transform variables for WSABI+.

nvars = size(x0,2);

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
%optimState.LBeps_orig = LB + (UB - LB)*options.TolBoundX;
%optimState.UBeps_orig = UB - (UB - LB)*options.TolBoundX;


%% Transform variables
% Transform variables
switch lower(options.BoundedTransform)
    case 'logit'
        trinfo = warpvars_laplace(nvars,LB,UB,PLB,PUB,3);
    case {'norminv','probit'}
        trinfo = warpvars_laplace(nvars,LB,UB,PLB,PUB,12);
    case 'student4'
        trinfo = warpvars_laplace(nvars,LB,UB,PLB,PUB,13);
    otherwise
        error('wsabi:UnknwonBoundedTransform','Unknown bounded transform.');
end
trinfo.x0_orig = x0;
if ~isfield(trinfo,'R_mat'); trinfo.R_mat = []; end
if ~isfield(trinfo,'scale'); trinfo.scale = []; end

optimState.LB = warpvars_laplace(LB,'dir',trinfo);
optimState.UB = warpvars_laplace(UB,'dir',trinfo);
optimState.PLB = warpvars_laplace(PLB,'dir',trinfo);
optimState.PUB = warpvars_laplace(PUB,'dir',trinfo);
optimState.x0 = warpvars_laplace(x0,'dir',trinfo);

% Fix infinities
optimState.x0(optimState.x0 == -Inf) = -10;
optimState.x0(optimState.x0 == Inf) = 10;

optimState.trinfo = trinfo;

%% Get warnings state

optimState.DefaultWarnings.singularMatrix = warning('query','MATLAB:singularMatrix');
warning('off',optimState.DefaultWarnings.singularMatrix.identifier);

end
