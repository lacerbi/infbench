function [y,y_std] = infbench_wood2010(x,infprob,mcmc_params)
%INFBENCH_GORIS2015 Inference benchmark log pdf -- neuronal model from Goris et al. (2015).

if nargin < 3; mcmc_params = []; end

if isempty(x)
    if isempty(infprob) % Generate this document        
        fprintf('\n');

        % Add sampler directory to MATLAB path
        pathstr = fileparts(mfilename('fullpath'));
        addpath([pathstr,filesep(),'..',filesep(),'infalgo',filesep(),'parallel-GP-SL',filesep(),'utils',filesep(),'mcmcstat-master']);
                
        for n = 1
            switch n
                case 1
                    % xmin = [138.8505216444501 2.38317564009549 0.682321320237262 1.1613095607596 1 0.231748337632257 -0.272638945416596 3.10117864852662 72.8822298534178 0.00789002312857097 0.101380347749067 0.693895739234024];
                    % fval = -2594.08310420223;
            end
            
            infprob = infbench_wood2010([],n);
            if isempty(mcmc_params); id = 0; else; id = mcmc_params(1); end

            sim_model = infprob.SimModel;
            D = sim_model.dim;            
            trinfo = sim_model.trinfo;
            
            % Prior used for sampling
            infprob.PriorMean = infprob.Prior.Mean;
            infprob.PriorVar = diag(infprob.Prior.Cov)';
            infprob.PriorVolume = prod(infprob.UB - infprob.LB);
            infprob.PriorType = 'uniform';            
            
            if id == 0
            
                % First, check optimum            
                opts = optimoptions('fminunc','Display','off','MaxFunEvals',700);

                x0 = xmin;
                x0 = warpvars(x0,'d',trinfo);   % Convert to unconstrained coordinates            
                [xnew,fvalnew] = fminunc(@(x) -infbench_goris2015(x,infprob), x0, opts);

                fvalnew = -fvalnew;
                xmin(infprob.idxParams) = warpvars(xnew,'inv',trinfo);
                fval = fvalnew + warpvars(xnew,'logp',trinfo);

                x0 = xmin(infprob.idxParams);
                x0 = warpvars(x0,'d',trinfo);   % Convert to unconstrained coordinates            
                [xnew,fvalnew] = fminunc(@(x) nlogpost(x,infprob), x0, opts);

                fvalnew = -fvalnew;
                xmin_post = xmin;
                xmin_post(infprob.idxParams) = warpvars(xnew,'inv',trinfo);
                fval_post = fvalnew + warpvars(xnew,'logp',trinfo);

                fprintf('\t\t\tcase %d\n',n);
                fprintf('\t\t\t\tname = ''%s'';\n\t\t\t\txmin = %s;\n\t\t\t\tfval = %s;\n',name,mat2str(xmin),mat2str(fval));
                fprintf('\t\t\t\txmin_post = %s;\n\t\t\t\tfval_post = %s;\n',mat2str(xmin_post),mat2str(fval_post));
                
            elseif id > 0 && n == mcmc_params(2)

                rng(id);    % Randomize seed
                
                % Define MCMC options
                if numel(mcmc_params) > 2
                    mcmc_opt.nsimu = mcmc_params(3);
                else
                    mcmc_opt.nsimu = 1e5;
                end
                mcmc_opt.nfinal = min(mcmc_opt.nsimu,1e3);
                mcmc_opt.display_type = 'on';
                
                infprob.Ns = 400;   % High precision
                
                model.ssfun = @(x,data) -2*logpost(x,infprob);
                model.N = 1;
                data = [];

                % Starting point
                x0 = infprob.Post.Mode;
                params = cell(D,1);
                for i = 1:D
                    params{i} = {sprintf('\\theta_{%d}',i),x0(i),infprob.LB(i),infprob.UB(i)};
                end

                % Additional MCMC settings
                options.nsimu = mcmc_opt.nsimu;
                widths = 0.1*(infprob.PUB - infprob.PLB);
                options.qcov = diag(widths.^2);
                options.method = 'am';
                options.updatesigma = 0;
                options.verbosity = 0; % no printing from mcmc
                options.waitbar = 0;

                % Run MCMC chains
                fprintf('Running MCMC...\n');
                tic
                [output,Xs,~,s2chain] = mcmcrun(model,data,params,options);
                toc
                exitflag = 0;
                
                % Remove burn-in
                Ns = size(Xs,1);
                Xs = Xs(ceil(Ns/2):Ns,:);
                
                % Thin remaining samples
                idx = round(linspace(1,size(Xs,1),mcmc_opt.nfinal))';                
                Xs = Xs(idx,:);
                
                % Re-evaluate final lls with very high precision
                fprintf('Re-evaluate log posteriors...\n');
                infprob.Ns = 8e3;
                lls = zeros(size(Xs,1));
                for i = 1:size(Xs,1)
                    lls(i) = logpost(Xs(i,:),infprob);
                    if mod(i,ceil(size(Xs,1)/10))==0; fprintf('%d/%d..',i,size(Xs,1)); end
                end
                fprintf('\nDone.\n');
                                
                filename = ['wood2010_mcmc_n' num2str(n) '_id' num2str(id) '.mat'];
                save(filename,'Xs','lls','exitflag','output');                
            end
            
        end
           
        
        
    else
        % Initialization call -- define problem and set up data
        n = infprob(1);
        
        % Are we transforming the entire problem to unconstrained space?
        transform_to_unconstrained_coordinates = false;
        
        % Add problem directory to MATLAB path
        % pathstr = fileparts(mfilename('fullpath'));
        % addpath([pathstr,filesep(),'synthlike']);
                
        xmin = NaN(1,3);       fval = Inf;
        xmin_post = NaN(1,3);  fval_post = Inf;
        Mean_mcmc = NaN(1,3);  Cov_mcmc = NaN(3,3);    lnZ_mcmc = NaN;

        % Define simulation model
        sim_model.name = 'Ricker';
        sim_model.theta_names = {'log(r)','\phi','\sigma_e'};
        sim_model.dim = 3;
        sim_model.gen = @(theta,n) simulate_ricker(theta,1,n); % generate one data set
        sim_model.n_data = 50; % T
        sim_model.summary_dim = 13;
        sim_model.comp_summaries = @(data,obs_data) ricker_summstats(data,obs_data);
                
        % Define parameter upper/lower bounds
        lb = [3,4,0];
        ub = [5,20,0.8];
        plb = (ub-lb)*0.1 + lb;
        pub = (ub-lb)*0.9 + lb;        
        noise = [];
        
        D = sim_model.dim;
        if transform_to_unconstrained_coordinates
            trinfo = warpvars(D,lb,ub,plb,pub);     % Transform to unconstrained space
            trinfo.mu = zeros(1,D);     % Necessary for retro-compatibility
            trinfo.delta = ones(1,D);
        else
            trinfo = [];
        end
        
        switch n
            case 1  % 
                theta_orig = [3.8,10,0.3];
                sim_model.data = [186 0 0 0 9 100 0 9 214 0 0 0 0 26 99 0 8 164 0 0 0 14 218 0 0 0 0 10 252 0 0 0 0 0 9 145 0 0 34 49 25 112 0 8 166 0 1 8 145 0]';
                % sim_model.gen(theta_orig, sim_model.n_data);
                xmin = [];
                fval = [];
                xmin_post = [];
                fval_post = [];
%                Mean_laplace = [-0.465520991260055 -1.69627543340411 -1.62447989201296 -0.767119415827957 0.470491515145298 -0.484169040138803 -2.60814931957517];
%                Cov_laplace = [0.000229133759218794 1.82326721309945e-05 0.000410425618637872 -0.000282902243846033 1.37333646747451e-05 2.55742161512486e-05 -7.49219827856342e-06;1.82326721309945e-05 0.00860578855990285 0.0113987299239016 -0.0046236785662946 -0.000242809250998626 -0.000356073871350522 -0.000143256002177701;0.000410425618637872 0.0113987299239016 0.162017238465403 -0.103918933733934 0.000992228291153627 0.00350094756853489 0.000154939078874722;-0.000282902243846033 -0.0046236785662946 -0.103918933733934 0.089505239660401 0.000464757180894073 -0.00336613410893957 -0.000435454797885497;1.37333646747451e-05 -0.000242809250998626 0.000992228291153627 0.000464757180894073 0.000982265536813824 0.00228454532591384 -3.74450494783302e-05;2.55742161512486e-05 -0.000356073871350522 0.00350094756853489 -0.00336613410893957 0.00228454532591384 0.00701913099415825 -8.12472287376342e-05;-7.49219827856342e-06 -0.000143256002177701 0.000154939078874722 -0.000435454797885497 -3.74450494783301e-05 -8.12472287376342e-05 0.00888756642652754];
%                lnZ_laplace = -2620.215196872;
                % R_max = 1.011. Ntot = 400000. Neff_min = 8013.5. Total funccount = 3124820.
        end
        
        sim_model.true_theta = theta_orig;
        
        if isempty(xmin); xmin = sim_model.true_theta; end
        if isempty(xmin_post); xmin_post = xmin; end       

        % data from the model and true posterior
        % sim_model.data = sim_model.gen(theta_orig, sim_model.n_data);
        sim_model.summary_true = ricker_summstats(sim_model.data,sim_model.data);
        sim_model.prior_eval = @(theta) ones(size(theta,1),1);
        sim_model.true_post_pdf = []; % sl_baseline_get(name)';
        
        xmin = warpvars(xmin,'d',trinfo);
        fval = fval + warpvars(xmin,'logp',trinfo);
        
        Mean = zeros(1,D);
        Cov = eye(D);
        Mode = xmin;
        
        y.Ns = 80;          % # samples for synthetic likelihood
                
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
        
        % Read marginals from file
        % marginals = load('goris2015_marginals.mat');
        % y.Post.MarginalBounds = marginals.MarginalBounds{n};
        % y.Post.MarginalPdf = marginals.MarginalPdf{n};
        
        % Save data and coordinate transformation struct
        sim_model.trinfo = trinfo;
        y.SimModel = sim_model;
        
        if isempty(fval) || isempty(fval_post)
            infprob_temp = y;
            infprob_temp.Ns = 2e3;            
            fval = infbench_wood2010(xmin,infprob_temp);
            fval = fval + warpvars(xmin,'logp',trinfo);
            fval_post = infbench_wood2010(xmin_post,infprob_temp);            
            fval_post = fval_post + warpvars(xmin_post,'logp',trinfo);
            y.ModeFval = fval;
            y.Post.ModeFval = fval_post;
        end
        
    end
    
else
    
    % Iteration call -- evaluate objective function
    
    % Transform unconstrained variables to original space
    x_orig = warpvars(x,'i',infprob.SimModel.trinfo);
    dy = warpvars(x,'logpdf',infprob.SimModel.trinfo);   % Jacobian correction
    
    sl_opt.N = infprob.Ns;
    sl_opt.estimator = 'sl';
    
    % Compute log likelihood of data and possibly std of log likelihood
    if nargout > 1
        [LL,y_std] = infbench_synthetic_loglik_estim(infprob.SimModel,sl_opt,x_orig);
    else
        LL = infbench_synthetic_loglik_estim(infprob.SimModel,sl_opt,x_orig);
    end    
    y = LL - dy;
    
end

end

%--------------------------------------------------------------------------
function y = logpost(x,infprob)
    y = infbench_wood2010(x,infprob);
    lnp = infbench_lnprior(x,infprob);
    y = y + lnp;
end

%--------------------------------------------------------------------------
function y = simulate_ricker(theta,N,T)
% Simulates one data set with length T from the Ricker model.
%
% INPUT:
% N - the starting population (equal to 1 in our application)
% T - the length of the data set
%
% OUTPUT:
% y - the simulated data set

r = exp(theta(1)); % the parameter we sample over is log(r)
phi = theta(2);
sigma_e = theta(3);

Ns = [N; zeros(T,1)];
srn = sigma_e*randn(T,1);
for t = 1:T
    Ns(t+1) = r*Ns(t)*exp(-Ns(t) + srn(t)); % population size
end
y = my_poissrnd(phi*Ns(2:end)); % the actual observed random variable

%y = y(51:end);
end


function r = my_poissrnd(lambda)
% Modified MATLAB poisson generator function: Removed input checkings etc. to speed up
% computations slightly. 
% TODO: This could make slightly faster by e.g. removing input checking from binornd.m
%
% Input lambda must be scalar or column vector with finite and non-negative values.
%
%   POISSRND uses a waiting time method for small values of LAMBDA,
%   and Ahrens' and Dieter's method for larger values of LAMBDA.

%   References:
%      [1]  Devroye, L. (1986) Non-Uniform Random Variate Generation,
%           Springer-Verlag.

%   Copyright 1993-2015 The MathWorks, Inc.

lambda(lambda < 0) = NaN;

%Initialize r to zero.
r = zeros(size(lambda));

r(isinf(lambda)) = Inf;

% For large lambda, use the method of Ahrens and Dieter as
% described in Knuth, Volume 2, 1998 edition.
k = find(15 <= lambda & lambda < Inf);
if ~isempty(k)
   alpha = 7/8;
   lk = lambda(k);
   m = floor(alpha * lk);

   % Generate m waiting times, all at once
   x = randg(m);
   t = (x <= lk);

   % If we did not overshoot, then the number of additional times
   % has a Poisson distribution with a smaller mean.
   r(k(t)) = m(t) + my_poissrnd(lk(t)-x(t));

   % If we did overshoot, then the times up to m-1 are uniformly
   % distributed on the interval to x, so the count of times less
   % than lambda has a binomial distribution.
   if ~all(t)
       r(k(~t)) = binornd(m(~t)-1, lk(~t)./x(~t));
   end
end

% For small lambda, generate and count waiting times.
j = find(lambda < 15);
p = zeros(numel(j),1,'like',lambda);
while ~isempty(j)
    p = p - log(rand(numel(j),1,'like',lambda));
    t = (p < lambda(j));
    j = j(t);
    p = p(t);
    r(j) = r(j) + 1;
end

% Return NaN if LAMBDA is negative.
r(isnan(lambda)) = NaN;
end

%--------------------------------------------------------------------------

function ss_x = ricker_summstats(x,y)
% Compute the summary statistics used in the Ricker model (as in Wood 2010)
%
% INPUT:
% x - simulated data
% y - observed data
%
% OUTPUT:
% ss_x - summary statistics for the simulated data x


% Autocovariances up to lag 5 
% [6 summaries]
[ss1x,lagsx] = xcov(x,5,'unbiased');
ind_ss1x = ((length(ss1x)+1)/2);
ss1x = ss1x(ind_ss1x: end);

% Coefficients of the cubic regression of order differences on their
% observed values
% [6 summaries]
order_diff = 1;
x_diff = diff(x, order_diff); y_diff = diff(y, order_diff);
x_diff = sort(x_diff); y_diff = sort(y_diff);
x_diff2 = x_diff - repmat(mean(x_diff),length(x_diff),1);
y_diff2 = y_diff - repmat(mean(y_diff),length(y_diff),1);
ss2x = regress(x_diff2, [y_diff2, y_diff2.^2, y_diff2.^3]);

%  Coefficients of the autoregression
%  y_{t+1}^{0.3} = beta_{1} y_{t}^{0.3} + beta_{2} y_{t}^{0.6} + epsilon_{t}
%  [2 summaries - beta_{1}, and beta_{2}]
x_mod = x.^0.3;
x_mod = x_mod - repmat(mean(x_mod),length(x_mod),1);
x_mod2 = x_mod(2:end); x_pred2 = x_mod(1:end-1);
ss3x = regress(x_mod2, [x_pred2, x_pred2.^2]);

% mean population (of y) [1 summary]
ss4x = mean(x); % Fine for mean of column [DEFAULT]

% no. of zeros observed [1 summary]
ss5x = sum(x == 0);

% Combining all summary statistics into a vector
ss_x = [ss1x;ss2x;ss3x;ss4x;ss5x];

% TESTING: what if one ignores some summaries...
%ss_x = [ss2x;ss3x;ss4x;ss5x];
end







