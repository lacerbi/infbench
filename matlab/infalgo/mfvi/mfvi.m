function [mu,sigma,elbo,output] = mfvi(fun,x0,PLB,PUB,LB,UB,options)
%MFVI Black-box mean field variational inference.

% Default options
defopts.Display             = 'iter';       % Print output
defopts.BoundedTransform    = 'logit';      % Input transform for bounded variables
defopts.Optimization        = 'bads';       % Optimization algorithm
defopts.SamplesPerIter      = [10,50,200,500]; % Number of samples per iteration
defopts.FinalSamples        = 1e4;          % Final samples to evaluate ELBO

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
optimState = setupvars_mfvi(x0,LB,UB,PLB,PUB,[],options,prnt);

% Define function in transformed space
fun_warped = @(x) fun(warpvars_mfvi(x,'i',optimState.trinfo)) + warpvars_mfvi(x,'logp',optimState.trinfo);

% Optimization
D = size(x0,2);

sigmas0 = 0.05*(optimState.PUB - optimState.PLB);
phi = [optimState.x0, log(sigmas0)];

lb = [optimState.LB, -Inf(1,D)];
ub = [optimState.UB, Inf(1,D)];
plb = [optimState.PLB, log(1e-3*(optimState.PUB - optimState.PLB))];
pub = [optimState.PUB, log(optimState.PUB - optimState.PLB)];

Ns_iter = options.SamplesPerIter;
xx = [];

% Run variational optimization
for iter = 1:numel(Ns_iter)
    
    if prnt; fprintf('Iteration %d.\n', iter); end
    
    Ns = Ns_iter(iter);
    xx = [xx; randn(Ns - size(xx,1),D)];
    target = @(phi_) negelbo_mfvi(phi_,fun_warped,xx);        

    phi0 = phi;
    
    switch lower(options.Optimization)
        case 'fminunc'
            opt_options = optimoptions('fminunc');
            opt_options.Display = options.Display;
            if iter < numel(Ns_iter)
                opt_options.MaxFunEvals = 50*D*iter;
            end            
            phi = fminunc(target,phi0,opt_options);
        case 'bads'
            opt_options = bads('defaults');
            opt_options.Display = options.Display;
            if iter < numel(Ns_iter)
                opt_options.MaxFunEvals = 50*D*iter;
            else
                opt_options.MaxFunEvals = 50*D;
                opt_options.AccelerateMeshSteps = 1;
            end
            opt_options.UncertaintyHandling = false;
            opt_options.TolFun = 0.02;
            opt_options.TolMesh = 1e-4;
            phi = bads(target,phi0,lb,ub,plb,pub,opt_options);
        otherwise
            error('Unknown optimization algorithm.');
    end    
end

% Final optimization (always fminunc)
% if prnt; fprintf('Iteration %d (FINAL).\n', numel(Ns_iter)); end
% Ns = Ns_iter(end);
% xx = [xx; randn(Ns - size(xx,1),D)];
% target = @(phi_) negelbo_mfvi(phi_,fun_warped,xx);        
% phi0 = phi;
% opt_options = optimoptions('fminunc');
% opt_options.Display = options.Display;
% opt_options.MaxFunEvals = 20*D;
% phi = fminunc(target,phi0,opt_options);

mu = phi(1:D);
sigma = exp(phi(D+(1:D)));

Ns_final = options.FinalSamples;
elbo = -negelbo_mfvi(phi,fun_warped,Ns_final);

% OUTPUT struct
output.trinfo = optimState.trinfo;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function add2path()
%ADD2PATH Adds MFVI subfolders to MATLAB path.

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
function optimState = setupvars_mfvi(x0,LB,UB,PLB,PUB,optimState,options,prnt)
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
        trinfo = warpvars_mfvi(nvars,LB,UB,PLB,PUB,3);
    case {'norminv','probit'}
        trinfo = warpvars_mfvi(nvars,LB,UB,PLB,PUB,12);
    case 'student4'
        trinfo = warpvars_mfvi(nvars,LB,UB,PLB,PUB,13);
    otherwise
        error('mfvi:UnknwonBoundedTransform','Unknown bounded transform.');
end
trinfo.x0_orig = x0;
if ~isfield(trinfo,'R_mat'); trinfo.R_mat = []; end
if ~isfield(trinfo,'scale'); trinfo.scale = []; end

optimState.LB = warpvars_mfvi(LB,'dir',trinfo);
optimState.UB = warpvars_mfvi(UB,'dir',trinfo);
optimState.PLB = warpvars_mfvi(PLB,'dir',trinfo);
optimState.PUB = warpvars_mfvi(PUB,'dir',trinfo);
optimState.x0 = warpvars_mfvi(x0,'dir',trinfo);

% Fix infinities
optimState.x0(optimState.x0 == -Inf) = -10;
optimState.x0(optimState.x0 == Inf) = 10;

optimState.trinfo = trinfo;

%% Get warnings state

optimState.DefaultWarnings.singularMatrix = warning('query','MATLAB:singularMatrix');
warning('off',optimState.DefaultWarnings.singularMatrix.identifier);

end
