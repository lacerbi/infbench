function [y,y_std] = infbench_wood2010(x,infprob,mcmc_params)
%INFBENCH_GORIS2015 Inference benchmark log pdf -- neuronal model from Goris et al. (2015).

if nargin < 3; mcmc_params = []; end

if isempty(x)
    if isempty(infprob) % Generate this document        
        fprintf('\n');

        for n = 1
            switch n
                case 1
                    xmin = [138.8505216444501 2.38317564009549 0.682321320237262 1.1613095607596 1 0.231748337632257 -0.272638945416596 3.10117864852662 72.8822298534178 0.00789002312857097 0.101380347749067 0.693895739234024];
                    fval = -2594.08310420223;
            end
            
            infprob = infbench_goris2015([],n);
            if isempty(mcmc_params); id = 0; else; id = mcmc_params(1); end
            
            trinfo = infprob.Data.trinfo;
            
            if id == 0
            
                % First, check optimum            
                opts = optimoptions('fminunc','Display','off','MaxFunEvals',700);

                x0 = xmin(infprob.idxParams);
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
                
                rng(id);
                widths = 0.5*(infprob.PUB - infprob.PLB);
                logpfun = @(x) -nlogpost(x,infprob);
                
                % Number of samples
                if numel(mcmc_params) > 2
                    W_mult = mcmc_params(3);
                else
                    W_mult = 200;
                end
                
                W = 2*(infprob.D+1);    % Number of walkers
                Ns = W*W_mult;             % Number of samples
                
                sampleopts.Burnin = Ns;
                sampleopts.Thin = 1;
                sampleopts.Display = 'iter';
                sampleopts.Diagnostics = false;
                sampleopts.VarTransform = false;
                sampleopts.InversionSample = false;
                sampleopts.FitGMM = false;
                sampleopts.TolX = 1e-5;
                % sampleopts.TransitionOperators = {'transSliceSampleRD'};

                x0 = xmin(infprob.idxParams);
                x0 = warpvars(x0,'d',trinfo);   % Convert to unconstrained coordinates
                LB = infprob.PLB - 10*widths;
                UB = infprob.PUB + 10*widths;
                
                [Xs,lls,exitflag,output] = eissample_lite(logpfun,x0,Ns,W,widths,LB,UB,sampleopts);
                
                filename = ['goris2015_mcmc_n' num2str(n) '_id' num2str(id) '.mat'];
                save(filename,'Xs','lls','exitflag','output');                
            end
            
        end
           
        
        
    else
        % Initialization call -- define problem and set up data
        n = infprob(1);
        
        % Add problem directory to MATLAB path
        pathstr = fileparts(mfilename('fullpath'));
        addpath([pathstr,filesep(),'synthlike']);
                
        xmin = NaN(1,3);       fval = Inf;
        xmin_post = NaN(1,3);  fval_post = Inf;
        Mean_mcmc = NaN(1,3);  Cov_mcmc = NaN(3,3);    lnZ_mcmc = NaN;

        % Define simulation model
        sim_model.name = 'Ricker';
        sim_model.theta_names = {'log(r)','\phi','\sigma_e'};
        sim_model.dim = 3;
        sim_model.gen = @(theta,n) simulate_ricker(theta,1,50); % generate one data set
        sim_model.n_data = 50; % T
        sim_model.summary_dim = 13;
        sim_model.comp_summaries = @(data,obs_data) ricker_summstats(data,obs_data);
                
        switch n
            case 1  % 
                sim_model.true_theta = [3.8,10,0.3];                
                xmin = [];
                fval = [];
                xmin_post = [];
                fval_post = [];
%                Mean_laplace = [-0.465520991260055 -1.69627543340411 -1.62447989201296 -0.767119415827957 0.470491515145298 -0.484169040138803 -2.60814931957517];
%                Cov_laplace = [0.000229133759218794 1.82326721309945e-05 0.000410425618637872 -0.000282902243846033 1.37333646747451e-05 2.55742161512486e-05 -7.49219827856342e-06;1.82326721309945e-05 0.00860578855990285 0.0113987299239016 -0.0046236785662946 -0.000242809250998626 -0.000356073871350522 -0.000143256002177701;0.000410425618637872 0.0113987299239016 0.162017238465403 -0.103918933733934 0.000992228291153627 0.00350094756853489 0.000154939078874722;-0.000282902243846033 -0.0046236785662946 -0.103918933733934 0.089505239660401 0.000464757180894073 -0.00336613410893957 -0.000435454797885497;1.37333646747451e-05 -0.000242809250998626 0.000992228291153627 0.000464757180894073 0.000982265536813824 0.00228454532591384 -3.74450494783302e-05;2.55742161512486e-05 -0.000356073871350522 0.00350094756853489 -0.00336613410893957 0.00228454532591384 0.00701913099415825 -8.12472287376342e-05;-7.49219827856342e-06 -0.000143256002177701 0.000154939078874722 -0.000435454797885497 -3.74450494783301e-05 -8.12472287376342e-05 0.00888756642652754];
%                lnZ_laplace = -2620.215196872;
                % R_max = 1.011. Ntot = 400000. Neff_min = 8013.5. Total funccount = 3124820.
        end
        
        if isempty(xmin); xmin = sim_model.true_theta; end
        if isempty(xmin_post); xmin_post = xmin; end       

        % data from the model and true posterior
        sim_model.data = sim_model.gen(sim_model.true_theta, sim_model.n_data);
        sim_model.summary_true = ricker_summstats(sim_model.data,sim_model.data);
        sim_model.prior_eval = @(theta) ones(size(theta,1),1);
        sim_model.true_post_pdf = []; % sl_baseline_get(name)';
        
        % Define parameter upper/lower bounds
        lb = [3,4,0];
        ub = [5,20,0.8];
        plb = (ub-lb)*0.1 + lb;
        pub = (ub-lb)*0.9 + lb;        
        noise = [];
        
        D = sim_model.dim;
        trinfo = warpvars(D,lb,ub,plb,pub);     % Transform to unconstrained space
        trinfo.mu = zeros(1,D);     % Necessary for retro-compatibility
        trinfo.delta = ones(1,D);
        y.xBaseFull = xmin;
        xmin = warpvars(xmin,'d',trinfo);
        fval = fval + warpvars(xmin,'logp',trinfo);
        
        Mean = zeros(1,D);
        Cov = eye(D);
        Mode = xmin;
        
        y.Ns = 80;          % # samples for synthetic likelihood
                
        y.D = D;
        y.LB = -Inf(1,D);   % Using unconstrained space
        y.UB = Inf(1,D);
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
        y.Post.xBaseFull = xmin_post;
        xmin_post = warpvars(xmin_post(idx_params),'d',trinfo);
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
            infprob_temp = infprob;
            infprob_temp.Ns = 2e3;            
            fval = infbench_wood2010(xmin,infprob_temp);
            fval = fval + warpvars(xmin,'logp',trinfo);
            fval_post = infbench_wood2010(xmpost,infprob_temp);            
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
        [LL,y_std] = noisy_loglik_estim(infprob.SimModel,sl_opt,x_orig);
    else
        LL = noisy_loglik_estim(infprob.SimModel,sl_opt,x_orig);
    end    
    y = LL - dy;
    
end

end

%--------------------------------------------------------------------------
function y = nlogpost(x,infprob)
    y = -infbench_goris2015(x,infprob);
    infprob.PriorMean = infprob.Prior.Mean;
    infprob.PriorVar = diag(infprob.Prior.Cov)';
    lnp = infbench_lnprior(x,infprob);
    y = y - lnp;
end