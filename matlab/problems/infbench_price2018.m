function [y,y_std] = infbench_price2018(x,infprob,mcmc_params)
%INFBENCH_PRICE2018 Inference benchmark log pdf -- g-and-k model from Price et al (2018).

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
            
            infprob = infbench_price2018([],n);
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
                    mcmc_opt.nsimu = 5e4;
                end
                mcmc_opt.nfinal = min(mcmc_opt.nsimu,1e3);
                mcmc_opt.display_type = 'on';
                
                infprob.Ns = 100;   % High precision
                
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
                infprob.Ns = 2e3;
                lls = zeros(size(Xs,1),1);
                for i = 1:size(Xs,1)
                    lls(i) = logpost(Xs(i,:),infprob);
                    if mod(i,ceil(size(Xs,1)/10))==0; fprintf('%d/%d..',i,size(Xs,1)); end
                end
                fprintf('\nDone.\n');
                                
                filename = ['price2018_mcmc_n' num2str(n) '_id' num2str(id) '.mat'];
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
                
        D = 4;
        xmin = NaN(1,D);       fval = Inf;
        xmin_post = NaN(1,D);  fval_post = Inf;
        Mean_mcmc = NaN(1,D);  Cov_mcmc = NaN(D,D);    lnZ_mcmc = NaN;

        % true data
        % data = load('data_gandk.mat'); % y
        ny = 1e4;                
        % [a,b,mut,sigmas] skew t MLEs based on observed data (Tskew_to_GandK_MLEs.mat)            
        Tskew = [2.39022814580771,0.56964234011095,2.30230171065523,0.363528274958011];
        
        % Define simulation model
        sim_model.name = 'g-and-k';
        sim_model.theta_names = {'a','b','g','k'};
        sim_model.dim = 4;
        sim_model.gen = @(theta,n) simulate_gandk(ny,theta); % generate one data set
        sim_model.n_data = 1; % T
        sim_model.summary_dim = 4;
        sim_model.comp_summaries = @(data,obs_data) Scores(Tskew(1),Tskew(2),Tskew(3),Tskew(4),data);
        sim_model.Tskew_to_GandK_MLEs = Tskew;
        sim_model.ny = ny;
                        
        % Define parameter upper/lower bounds
        lb = [2.5,0.5,1.5,0.3];
        ub = [3.5,1.5,2.5,0.7];
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
                theta_orig = [2.99,1,2.08,0.5];
                xmin = [];
                fval = [];
                xmin_post = [];
                fval_post = [];
                % R_max = 1.000. Ntot = 100000. Neff_min = 96172.7. Total funccount = 10000000.
            	%Mean_mcmc = [4.12423365181945 9.29185965281371 0.148350792260269];
            	%Cov_mcmc = [0.013081294550676 -0.0351307243157069 -0.00634174005917681;-0.0351307243157069 0.173607223018992 0.0151707836408846;-0.00634174005917681 0.0151707836408846 0.0104764488100004];
                %lnZ_mcmc = -23.501766415079;
        end
        
        sim_model.true_theta = theta_orig;
        
        if isempty(xmin); xmin = sim_model.true_theta; end
        if isempty(xmin_post); xmin_post = xmin; end       

        % data from the model and true posterior
        sim_model.data = []; % NOTE: true data is always used and loaded from file
        sim_model.summary_true = [0,0,0,0];
        sim_model.prior_eval = @(theta) ones(size(theta,1),1);
        sim_model.true_post_pdf = []; % sl_baseline_get(name)';
        
        xmin = warpvars(xmin,'d',trinfo);
        fval = fval + warpvars(xmin,'logp',trinfo);
        
        Mean = zeros(1,D);
        Cov = eye(D);
        Mode = xmin;
        
        y.Ns = 20;          % # samples for synthetic likelihood
                
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
        %marginals = load('price2018_marginals.mat');
        %y.Post.MarginalBounds = marginals.MarginalBounds{n};
        %y.Post.MarginalPdf = marginals.MarginalPdf{n};
        
        % Save data and coordinate transformation struct
        sim_model.trinfo = trinfo;
        y.SimModel = sim_model;
        
        if isempty(fval) || isempty(fval_post)
            infprob_temp = y;
            infprob_temp.Ns = 2e3;            
            fval = infbench_price2018(xmin,infprob_temp);
            fval = fval + warpvars(xmin,'logp',trinfo);
            fval_post = infbench_price2018(xmin_post,infprob_temp);            
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
    y = infbench_price2018(x,infprob);
    lnp = infbench_lnprior(x,infprob);
    y = y + lnp;
end

%--------------------------------------------------------------------------
% Following functions from Price et al. (2018)

function q = simulate_gandk(n,parms)
% simulate_gandk simulates a data set of length n from the g-and-k distribution
%
% INPUT:
% n - number of observations
% parms - vector [a b g k] of g-and-k distribution
%
% OUTPUT:
% q - n draws from g-and-k distribution with parameters parms.

zu = randn(n,1);
q = fun_gandk(parms,zu); %the g-and-k quantile function
end

function f = fun_gandk(parms,zu)
% fun_gandk is the g-and-k quantile function
%
% INPUT:
% parms - vector [a b g k] of g-and-k distribution
% zu - quantiles of the standard normal distribution (vector)
%
% OUTPUT:
% f - the calculated quantiles

a = parms(1); b = parms(2); c = 0.8; g = parms(3); k = parms(4);
f = a + b*(1 + c*(1-exp(-g*zu))./(1 + exp(-g*zu))).*(1 + zu.^2).^k.*zu;

end

function [score_results] = Scores(a,b,mut,sigmas,x)
% Scores computes the score for the simulated data, with the skew t
% parameters set to the MLEs for the observed data.
%
% INPUT:
% a - MLE estimate for skew t parameter 'a' based on observed data
% b - MLE estimate for skew t parameter 'b' based on observed data
% mut - MLE estimate for skew t parameter 'mu' based on observed data
% sigmas - MLE estimate for skew t parameter 'sigma' based on observed data
% x - simulated data
%
% OUTPUT:
% score_results - the vector of scores

m = length(x);

%dlogL_da = -m*log(2)-m*psi(a)+m*psi(a+b)-m/2/(a+b)+sum(d1_da+d2_da);
dlogL_da = -m*log(2)-m*psi(a)+m*psi(a+b)-m/2/(a+b)+sum(...
    log(1 - (mut - x)./(sigmas*(a + b + (mut - x).^2./sigmas^2).^(1/2))) - ((mut - x).*(a + 1/2))./(2*sigmas*((mut - x)./(sigmas*(a + b + (mut - x).^2./sigmas^2).^(1/2)) - 1).*(a + b + (mut - x).^2/sigmas^2).^(3/2))...
    -((mut - x).*(b + 1/2))./(2*sigmas*((mut - x)./(sigmas*(a + b + (mut - x).^2./sigmas^2).^(1/2)) + 1).*(a + b + (mut - x).^2./sigmas^2).^(3/2)));

%dlogL_db = -m*log(2)-m*psi(b)+m*psi(a+b)-m/2/(a+b)+sum(d1_db+d2_db);
dlogL_db = -m*log(2)-m*psi(b)+m*psi(a+b)-m/2/(a+b)+sum(...
    -((mut - x).*(a + 1/2))./(2*sigmas*((mut - x)./(sigmas*(a + b + (mut - x).^2./sigmas^2).^(1/2)) - 1).*(a + b + (mut - x).^2./sigmas^2).^(3/2))...
    +log((mut - x)./(sigmas*(a + b + (mut - x).^2./sigmas^2).^(1/2)) + 1) - ((mut - x).*(b + 1/2))./(2*sigmas*((mut - x)./(sigmas*(a + b + (mut - x).^2./sigmas^2).^(1/2)) + 1).*(a + b + (mut - x).^2./sigmas^2).^(3/2)));

%dlogL_dmu = sum(d1_dmu+d2_dmu);
dlogL_dmu = sum(((a + 1/2)*(1./(sigmas*(a + b + (mut - x).^2./sigmas^2).^(1/2)) - ((2*mut - 2.*x).*(mut - x))./(2*sigmas^3*(a + b + (mut - x).^2./sigmas^2).^(3/2))))./((mut - x)./(sigmas*(a + b + (mut - x).^2./sigmas^2).^(1/2)) - 1)...
    +((b + 1/2)*(1./(sigmas*(a + b + (mut - x).^2/sigmas^2).^(1/2)) - ((2*mut - 2.*x).*(mut - x))./(2*sigmas^3*(a + b + (mut - x).^2/sigmas^2).^(3/2))))./((mut - x)./(sigmas*(a + b + (mut - x).^2./sigmas^2).^(1/2)) + 1));

%dlogL_dsigma = -m/sigmas+sum(d1_dsigma+d2_dsigma);
dlogL_dsigma = -m/sigmas+sum(-((a + 1/2)*((mut - x)./(sigmas^2*(a + b + (mut - x).^2./sigmas^2).^(1/2)) - (mut - x).^3./(sigmas^4*(a + b + (mut - x).^2./sigmas^2).^(3/2))))./((mut - x)./(sigmas*(a + b + (mut - x).^2./sigmas^2).^(1/2)) - 1)...
    -((b + 1/2)*((mut - x)./(sigmas^2*(a + b + (mut - x).^2./sigmas^2).^(1/2)) - (mut - x).^3./(sigmas^4*(a + b + (mut - x).^2./sigmas^2).^(3/2))))./((mut - x)./(sigmas*(a + b + (mut - x).^2./sigmas^2).^(1/2)) + 1));

score_results = [dlogL_da dlogL_db dlogL_dmu dlogL_dsigma];

end