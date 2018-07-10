function [vbmodel,exitflag,output] = agp_lite(fun,x0,LB,UB,PLB,PUB,options)
%AGP_LITE Light implementation of AGP method for Bayesian inference.

% This is a variant of the AGP algorithm proposed in:
% Wang, H., & Li, J. (2017). Adaptive Gaussian process approximation for 
% Bayesian inference with expensive likelihood functions. 
% arXiv preprint arXiv:1703.09930.

% Code by Luigi Acerbi (2018).

% Add variational Gaussian mixture model toolbox to path
mypath = fileparts(mfilename('fullpath'));
addpath([mypath filesep 'vbgmm']);

% Check existence of GPlite toolbox in path
if ~exist('gplite_train.m','file')
    error('agp_lite:NoGPliteToolbox','The GPlite toolbox needs to be in the MATLAB path to run AGP_LITE.');
end

exitflag = 0;   % To be used in the future
D = size(x0,2);

% Check hard and plausible bounds
if isempty(LB); LB = -Inf; end
if isempty(UB); UB = Inf; end
if isscalar(LB); LB = LB*ones(1,D); end
if isscalar(UB); UB = UB*ones(1,D); end
if isempty(PLB); PLB = LB; end
if isempty(PUB); PUB = UB; end
if isscalar(PLB); PLB = PLB*ones(1,D); end
if isscalar(PUB); PUB = PUB*ones(1,D); end
if any(~isfinite(PLB)) || any(~isfinite(PUB))
    error('agp_lite:NotFinitePB','Plausible lower/upper bounds need to be finite.');
end

%% Algorithm options and defaults

defopts.MaxFunEvals = D*100;            % Maximum number of fcn evals
defopts.Ns = 2e4;                       % Number of MCMC samples to get approximate posterior
defopts.NsMax_gp = 0;                   % Max GP hyperparameter samples (0 = optimize)
defopts.StopGPSampling = 200 + 10*D;    % Training set size for switching to GP hyperparameter optimization (if sampling)
defopts.AcqFun = @acqagp;               % AGP acquisition function
defopts.SamplingMethod = 'parallel';    % MCMC sampler for approximate posterior
defopts.Plot = 0;                       % Make diagnostic plots at each iteration
defopts.NcompMax = 30;                  % Maximum number of mixture components
defopts.FracExpand = 0.1;               % Expand search box by this amount
defopts.ProposalFcn = @(x) agp_proposal(x,PLB,PUB); % Proposal fcn based on PLB and PUB (unused)

for f = fields(defopts)'
    if ~isfield(options,f{:}) || isempty(options.(f{:}))
        options.(f{:}) = defopts.(f{:});
    end
end

% Convert AcqFun to function handle if passed as a string
if ischar(options.AcqFun); options.AcqFun = str2func(options.AcqFun); end

Ninit = 20;     % Initial design
Nstep = 10;     % Training points at each iteration
Nsearch = 2^13; % Starting search points for acquisition fcn

% Variational Bayesian Gaussian mixture options
vbopts.Display     = 'off';        % No display
vbopts.TolBound    = 1e-8;         % Minimum relative improvement on variational lower bound
vbopts.Niter       = 2000;         % Maximum number of iterations
vbopts.Nstarts     = 1;            % Number of runs
vbopts.TolResponsibility = 0.5;    % Remove components with less than this total responsibility
vbopts.ClusterInit = 'kmeans';     % Initialization of VB (methods are 'rand' and 'kmeans')

% GPLITE model options
gpopts.Nopts = 1;       % Number of hyperparameter optimization runs
gpopts.Ninit = 2^10;    % Initial design size for hyperparameter optimization
gpopts.Thin = 5;        % Thinning for hyperparameter sampling (if sampling)
gp_meanfun = 'zero';    % Constant-zero mean function

% Setup options for CMA-ES optimization
cmaes_opts = cmaes_modded('defaults');
cmaes_opts.EvalParallel = 'yes';
cmaes_opts.DispFinal = 'off';
cmaes_opts.SaveVariables = 'off';
cmaes_opts.DispModulo = Inf;
cmaes_opts.LogModulo = 0;
cmaes_opts.LBounds = LB(:);
cmaes_opts.UBounds = UB(:);
cmaes_opts.CMA.active = 1;      % Use Active CMA (generally better)

%% Initialization

% Evaluate fcn on random starting grid
Nrnd = Ninit - size(x0,1);
Xrnd = bsxfun(@plus,PLB,bsxfun(@times,PUB-PLB,rand(Nrnd,D)));
X = [x0;Xrnd];
X = bsxfun(@min,bsxfun(@max,LB,X),UB);  % Force X inside hard bounds
y = zeros(Ninit,1);
for i = 1:Ninit; y(i) = fun(X(i,:)); end

% Draw initial samples
mu0 = 0.5*(PLB + PUB);
width = 0.5*(PUB - PLB);
sigma0 = width;
Xs = bsxfun(@plus,bsxfun(@times,sigma0,randn(options.Ns,D)),mu0);
Xs = bsxfun(@min,bsxfun(@max,LB,Xs),UB);  % Force Xs inside hard bounds

% Fit single Gaussian to initial pdf as a VBGMM object
vbmodel = vbgmmfit(Xs',1,[],vbopts);

% Starting GP hyperparameter vector
hyp = [log(width(:));log(std(y));log(1e-3)];

%% Main loop
iter = 1;
while 1
    fprintf('Iter %d...', iter);
    N = size(X,1);
    
    % Build GP approximation
    fprintf(' Building GP approximation...');
    Ns_gp = min(round(options.NsMax_gp/10),round(options.NsMax_gp / sqrt(N)));
    if N >= options.StopGPSampling; Ns_gp = 0; end
    py = vbgmmpdf(vbmodel,X');      % Evaluate approximation at X    
    y_gp = y - log(py(:));          % Log difference
    
    % Train GP
    hypprior = getGPhypprior(X);    % Get prior over GP hyperparameters
    [gp,hyp] = gplite_train(hyp,Ns_gp,X,y_gp,gp_meanfun,hypprior,[],gpopts);
    
    % Sample from GP plus log approximation
    fprintf(' Sampling from GP...');
    lnpfun = @(x) lnprior(x,vbmodel,LB,UB);
    Xs = gplite_sample(gp,options.Ns,[],options.SamplingMethod,lnpfun);
     
    % Plot current approximate posterior and training points
    if options.Plot
        %Xrnd = vbgmmrnd(vbmodel,1e4)';
        %cornerplot(Xrnd);    
        for i = 1:D; names{i} = ['x_{' num2str(i) '}']; end
        [~,ax] = cornerplot(Xs,names);
        for i = 1:D-1
            for j = i+1:D
                axes(ax(j,i));  hold on;
                scatter(X(:,i),X(:,j),'ok');
            end
        end
        drawnow;
    end
    
    % Refit vbGMM
    fprintf(' Refit vbGMM...\n');
    vbmodel = vbgmmfit(Xs',options.NcompMax,[],vbopts);

    %Xrnd = vbgmmrnd(vbmodel,1e5)';
    %Mean = mean(Xrnd,1);
    %Cov = cov(Xrnd);
    
    % Estimate normalization constant in HPD region
    [lnZ,lnZ_var] = estimate_lnZ(X,y,vbmodel);
        
    fprintf('Estimate of lnZ = %f +/- %f.\n',lnZ,sqrt(lnZ_var));
    
    % Record stats
    stats(iter).N = N;
    % stats(iter).Mean = Mean;
    % stats(iter).Cov = Cov;
    stats(iter).lnZ = lnZ;
    stats(iter).lnZ_var = lnZ_var;
    stats(iter).vbmodel = vbmodel;
    stats(iter).gp = gplite_clean(gp);
    
    % Find max of approximation among GP samples and record approximate mode
    ys = vbgmmpdf(vbmodel,Xs')';
    [~,idx] = max(ys);
    stats(iter).mode = Xs(idx,:);
    
    if N >= options.MaxFunEvals; break; end
    
    % Select new points
    fprintf(' Active sampling...');
    for iNew = 1:Nstep
        fprintf(' %d..',iNew);
        % Random uniform search inside search box
        width = max(X) - min(X);
        lb = min(X) - width*options.FracExpand; ub = max(X) + width*options.FracExpand;
        lb = max(LB,min(lb,PLB)); ub = min(UB,max(ub,PUB));
        [xnew,fval] = fminfill(@(x) options.AcqFun(x,vbmodel,gp,options),[],[],[],lb,ub,[],struct('FunEvals',floor(Nsearch/2)));
        
        % Random search sample from vbGMM
        xrnd = vbgmmrnd(vbmodel,ceil(Nsearch/2))';
        xrnd = bsxfun(@min,bsxfun(@max,LB,xrnd),UB);  % Force Xs inside hard bounds
        frnd = options.AcqFun(xrnd,vbmodel,gp,options);
        [frnd_min,idx] = min(frnd);        
        if frnd_min < fval; xnew = xrnd(idx,:); fval = frnd_min; end

        % Optimize from best point with CMA-ES
        insigma = (max(X) - min(X))'/sqrt(3);
        [xnew_cmaes,fval_cmaes] = cmaes_modded(func2str(options.AcqFun),xnew',insigma,cmaes_opts,vbmodel,gp,options,1);
        if fval_cmaes < fval; xnew = xnew_cmaes'; end
        
        % Add point
        ynew = fun(xnew);
        X = [X; xnew];
        y = [y; ynew];
        
        py = vbgmmpdf(vbmodel,xnew');   % Evaluate approximation at X    
        ynew_gp = ynew - log(py);       % Log difference
        gp = gplite_post(gp,xnew,ynew_gp,[],1);   % Rank-1 update
        
        % Plot new candidate points
        if options.Plot
            for i = 1:D-1
                for j = i+1:D
                    axes(ax(j,i));  hold on;
                    scatter(xnew(:,i),xnew(:,j),'or','MarkerFaceColor','r');
                end
            end
            drawnow;
        end

    end
    fprintf('\n');
    
    iter = iter + 1;
end

output.X = X;
output.y = y;
output.stats = stats;

end

%--------------------------------------------------------------------------
function lp = lnprior(x,vbmodel,LB,UB)
%LNPRIOR Log prior and base function for approximate posterior.

lp = log(vbgmmpdf(vbmodel,x'))';
if any(isfinite(LB)) || any(isfinite(UB))
    idx = any(bsxfun(@gt,x,UB),2) | any(bsxfun(@lt,x,LB),2);
    lp(idx) = -Inf;
end

end

%--------------------------------------------------------------------------
function [lnZ,lnZ_var] = estimate_lnZ(X,y,vbmodel)
%ESTIMATE_LNZ Rough approximation of normalization constant

hpd_frac = 0.2;     % Top 20% HPD
N = size(X,1);

lp = log(vbgmmpdf(vbmodel,X')');

% Take HPD points according to both fcn samples and model
[~,ord] = sort(lp + y,'descend');

idx_hpd = ord(1:ceil(N*hpd_frac));
lp_hpd = lp(idx_hpd);
y_hpd = y(idx_hpd);

delta = -(lp_hpd - y_hpd);

lnZ = mean(delta);
lnZ_var = var(delta)/numel(delta);

end

%--------------------------------------------------------------------------
function hypprior = getGPhypprior(X)
%GETGPHYPPRIOR Define empirical Bayes prior over GP hyperparameters.

D = size(X,2);
hypprior = [];
Nhyp = D+2;
hypprior.mu = NaN(1,Nhyp);
hypprior.sigma = NaN(1,Nhyp);
hypprior.df = 3*ones(1,Nhyp);    % Broad Student's t prior
hypprior.mu(1:D) = log(std(X));
hypprior.sigma(1:D) = max(2,log(max(X)-min(X)) - log(std(X)));
hypprior.mu(D+2) = log(1e-2);
hypprior.sigma(D+2) = 0.5;

end