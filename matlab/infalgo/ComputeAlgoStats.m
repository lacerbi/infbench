function [gsKL,Mean,Cov,lnZ,lnZ_var,Mode,MTV] = ComputeAlgoStats(X,y,probstruct,compute_lnZ,Ns_moments,gp,s2)
%COMPUTEALGOSTATS Compute GP model-based statistics from given training set.
    
if nargin < 4 || isempty(compute_lnZ); compute_lnZ = false; end
if nargin < 5 || isempty(Ns_moments); Ns_moments = 2e4; end
if nargin < 6; gp = []; end
if nargin < 7; s2 = []; end

try

    % Add prior to y if not previously added
    if ~probstruct.AddLogPrior
        lnp = infbench_lnprior(X,probstruct);
        y = y + lnp;
    end

    if isempty(gp)
        % Compute Gaussianized symmetric KL-divergence with ground truth
        gp.X = X;
        gp.y = y;
        if ~isempty(s2); gp.s2 = s2; end
        gp.meanfun = 4;         % Negative quadratic mean fcn
        gp.covfun = 1;          % SE ard covariance
        if ~isempty(s2)
            gp.noisefun = [1 1];  % Provided noise            
        else
            gp.noisefun = [1 0];  % Gaussian observation noise
        end
    elseif isstruct(gp)
        gp = gplite_post(gp);   % Fill in GP
    end    

    if isstruct(gp)
        gp.Nopts = 2;   % Perform two restarts
        xx = gplite_sample(gp,Ns_moments,[],'parallel');
    else
        xx = gp;    % GP was passed as samples
    end
    Mean = mean(xx,1);
    Cov = cov(xx);
    [kl1,kl2] = mvnkl(Mean,Cov,probstruct.Post.Mean,probstruct.Post.Cov);
    gsKL = 0.5*(kl1 + kl2);

    if compute_lnZ
        NcompMax = 30;  % Max number of mixture components
        Nsvb_max = 5e3; % Max number of samples for VBGMM

        % Add variational Gaussian mixture model toolbox to path
        mypath = fileparts(mfilename('fullpath'));
        addpath([mypath filesep 'vbgmm']);       

        % Variational Bayesian Gaussian mixture options
        vbopts.Display     = 'off';        % No display
        vbopts.TolBound    = 1e-8;         % Minimum relative improvement on variational lower bound
        vbopts.Niter       = 2000;         % Maximum number of iterations
        vbopts.Nstarts     = 1;            % Number of runs
        vbopts.TolResponsibility = 0.5;    % Remove components with less than this total responsibility
        vbopts.ClusterInit = 'kmeans';     % Initialization of VB (methods are 'rand' and 'kmeans')

        idx = round(linspace(1,size(xx,1),min(Nsvb_max,size(xx,1))));    
        Xs = xx(idx,:);

        vbmodel = vbgmmfit(Xs',NcompMax,[],vbopts);
        [lnZ,lnZ_var] = estimate_lnZ(X,y,vbmodel);
    else
        lnZ = [];   lnZ_var = [];
    end

    if nargout > 5
        if isstruct(gp)
            Mode = gplite_fmin(gp,[],1);    % Max flag - finds maximum
        else
            Mode = NaN(1,size(xx,2));       % No mode for now with samples
        end
    end
    
    % Compute Marginal Total Variation (MTV)
    if nargout > 6
        MTV = ComputeMarginalTotalVariation(xx,probstruct);
    end
    
catch ME    
    msgText = getReport(ME);
    warning(['Cannot compute GP-based statistics: ' msgText]);
    
    % Try continuing, return NaNs
    D = size(X,2);
    gsKL = NaN;
    Mean = NaN(1,D);
    Cov = NaN(D,D);
    lnZ = NaN;
    lnZ_var = NaN;
    Mode = NaN(1,D);
    MTV = NaN(1,D);
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

idx = isfinite(delta);

if any(idx)
    lnZ = mean(delta(idx));
    lnZ_var = var(delta(idx))/numel(delta(idx));
else
    lnZ = NaN;
    lnZ_var = NaN;
end

end
