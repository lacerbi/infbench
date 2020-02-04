function [y,y_std] = infbench_acerbidokka2018(x,infprob,mcmc_params)
%INFBENCH_ACERBIDOKKA2018 Inference benchmark log pdf -- simple causal inference model from Acerbi, Dokka et al. (2018).

if nargin < 3; mcmc_params = []; end

if isempty(x)
    if isempty(infprob) % Generate this document        
        fprintf('\n');

        % Add sampler directory to MATLAB path
        pathstr = fileparts(mfilename('fullpath'));
        addpath([pathstr,filesep(),'..',filesep(),'infalgo',filesep(),'parallel-GP-SL',filesep(),'utils',filesep(),'mcmcstat-master']);
                
        for n = 1:3
            switch n
                case 1
                    % xmin = [138.8505216444501 2.38317564009549 0.682321320237262 1.1613095607596 1 0.231748337632257 -0.272638945416596 3.10117864852662 72.8822298534178 0.00789002312857097 0.101380347749067 0.693895739234024];
                    % fval = -2594.08310420223;
            end
            
            name = ['S' num2str(n)];
            
            infprob = infbench_acerbidokka2018([],n);
            infprob.DeterministicFlag = true;
            if isempty(mcmc_params); id = 0; else; id = mcmc_params(1); end

            D = infprob.D;
            trinfo = infprob.Data.trinfo;
            
            % Prior used for sampling
            infprob.PriorMean = infprob.Prior.Mean;
            infprob.PriorVar = diag(infprob.Prior.Cov)';
            infprob.PriorVolume = prod(infprob.UB - infprob.LB);
            infprob.PriorType = 'uniform';            

            LB = infprob.LB;
            UB = infprob.UB;
            PLB = infprob.PLB;
            PUB = infprob.PUB;
            
            if id == 0
            
                % Compute optimum
                Nopts = 5;
                
                opts = optimoptions('fmincon','Display','off','MaxFunEvals',5e4);

                for iOpt = 1:Nopts                
                    x0 = rand(1,D).*(PUB - PLB) + PLB;
                    fun = @(x) -infbench_acerbidokka2018(x,infprob);                    
                    [xnew(iOpt,:),fvalnew(iOpt)] = fmincon(fun,x0,[],[],[],[],LB,UB,[],opts);
                end
                
                [fvalnew,idx_best] = min(fvalnew);
                xnew = xnew(idx_best,:);
                
                fvalnew = -fvalnew;
                xmin = warpvars(xnew,'inv',trinfo);
                fval = fvalnew + warpvars(xnew,'logp',trinfo);

                x0 = xmin;
                x0 = warpvars(x0,'d',trinfo);   % Convert to unconstrained coordinates            
                fun = @(x) -logpost(x,infprob);
                [xnew,fvalnew] = fmincon(fun,x0,[],[],[],[],LB,UB,[],opts);

                fvalnew = -fvalnew;
                xmin_post = xmin;
                xmin_post = warpvars(xnew,'inv',trinfo);
                fval_post = fvalnew + warpvars(xnew,'logp',trinfo);

                fprintf('\t\t\tcase %d\n',n);
                fprintf('\t\t\t\tname = ''%s'';\n\t\t\t\txmin = %s;\n\t\t\t\tfval = %s;\n',name,mat2str(xmin),mat2str(fval));
                fprintf('\t\t\t\txmin_post = %s;\n\t\t\t\tfval_post = %s;\n',mat2str(xmin_post),mat2str(fval_post));
                
            elseif id > 0

                rng(id);
                widths = 0.5*(PUB - PLB);
                logpfun = @(x) logpost(x,infprob);
                
                % Number of samples
                if numel(mcmc_params) > 1
                    Ns = mcmc_params(2);
                else
                    Ns = 1e3;
                end
                
                W = 2*(infprob.D+1);    % Number of walkers
                
                sampleopts.Thin = 23;
                sampleopts.Burnin = Ns*sampleopts.Thin;
                sampleopts.Display = 'notify';
                sampleopts.Diagnostics = false;
                sampleopts.VarTransform = false;
                sampleopts.InversionSample = false;
                sampleopts.FitGMM = false;
                sampleopts.TolX = 1e-5;
                % sampleopts.TransitionOperators = {'transSliceSampleRD'};

                x0 = infprob.Post.Mode;
                x0 = warpvars(x0,'d',trinfo);
                
                [Xs,lls,exitflag,output] = eissample_lite(logpfun,x0,Ns,W,widths,LB,UB,sampleopts);
                
                filename = ['acerbidokka2018_mcmc_n' num2str(n) '_id' num2str(id) '.mat'];
                save(filename,'Xs','lls','exitflag','output');                
            end
            
        end
        
    else
        % Initialization call -- define problem and set up data
        n = infprob(1);
        
        % Are we transforming the entire problem to unconstrained space?
        transform_to_unconstrained_coordinates = false;
        
        % Add problem directory to MATLAB path
        pathstr = fileparts(mfilename('fullpath'));
        addpath([pathstr,filesep(),'acerbidokka2018']);
                
        D = 6;        
        xmin = NaN(1,D);       fval = Inf;
        xmin_post = NaN(1,D);  fval_post = Inf;
        Mean_mcmc = NaN(1,D);  Cov_mcmc = NaN(D,D);    lnZ_mcmc = NaN;
                
        % Define parameter upper/lower bounds
        lb = [log(0.5)*ones(1,4),5e-3,log(0.25)];
        ub = [log(80)*ones(1,4),0.5,log(180)];
        plb = [log(1)*ones(1,4),0.01,log(1)];
        pub = [log(40)*ones(1,4),0.2,log(45)];
        noise = [];
        
        if transform_to_unconstrained_coordinates
            trinfo = warpvars(D,lb,ub,plb,pub);     % Transform to unconstrained space
            trinfo.mu = zeros(1,D);     % Necessary for retro-compatibility
            trinfo.delta = ones(1,D);
        else
            trinfo = [];
        end
        
        switch n
			case 1
				name = 'S1';
				xmin = [1.64206809048607 1.56382540104707 2.26143803248302 1.66865229743062 0.0252736737711433 2.31634928106453];
				fval = -483.513343605144;
				xmin_post = [1.6421811436931 1.56395771711076 2.26147071487172 1.66854514580804 0.0252735950695666 2.31634927169251];
				fval_post = -491.191510124604;
			case 2
				name = 'S2';
				xmin = [1.38726609391213 2.60712631868221 3.08537558206964 1.5047752731606 0.143318676649389 2.99303342142661];
				fval = -429.775221706332;
				xmin_post = [1.38726676853333 2.60712639983091 3.08537554925837 1.50477233430562 0.143318797696291 2.99303327788359];
				fval_post = -437.453388225819;
			case 3
				name = 'S3';
				xmin = [1.19318104591285 1.98890187717717 2.6580800074951 2.02743456496656 0.0867816136434034 2.88917816646634];
				fval = -583.757699194538;
				xmin_post = [1.19394210042755 1.98905711125483 2.65812076676506 2.02729160630505 0.086781453413982 2.88917822194199];
				fval_post = -591.435865713995;                
                % R_max = 1.000. Ntot = 100000. Neff_min = 96172.7. Total funccount = 10000000.
%            	Mean_mcmc = [4.12423365181945 9.29185965281371 0.148350792260269];
%            	Cov_mcmc = [0.013081294550676 -0.0351307243157069 -0.00634174005917681;-0.0351307243157069 0.173607223018992 0.0151707836408846;-0.00634174005917681 0.0151707836408846 0.0104764488100004];
%                lnZ_mcmc = -23.501766415079;
        end

        temp = load('acerbidokka2018_data.mat');
        data.X = temp.unity_data{n};  % Get subject's dataset
                
        xmin = warpvars(xmin,'d',trinfo);
        fval = fval + warpvars(xmin,'logp',trinfo);
        
        Mean = zeros(1,D);
        Cov = eye(D);
        Mode = xmin;
                
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
        
        
        data.IBSNreps = 10; % Reps used for IBS estimator
        
        % Read marginals from file
%        marginals = load('acerbidokka2018_marginals.mat');
%        y.Post.MarginalBounds = marginals.MarginalBounds{n};
%        y.Post.MarginalPdf = marginals.MarginalPdf{n};
        
        % Save data and coordinate transformation struct
        data.trinfo = trinfo;
        y.Data = data;
        y.DeterministicFlag = false;
                
    end
    
else
    
    % Iteration call -- evaluate objective function
    
    % Transform unconstrained variables to original space
    x_orig = warpvars(x,'i',infprob.Data.trinfo);
    dy = warpvars(x,'logpdf',infprob.Data.trinfo);   % Jacobian correction
    
    % Compute log likelihood of data and possibly std of log likelihood
    if infprob.DeterministicFlag
        LL = ll_deterministic(x_orig,infprob.Data);
        y_std = 0;
    else
        [LL,y_std] = ll_stochastic(x_orig,infprob.Data);
    end
    y = LL - dy;
    
end

end

%--------------------------------------------------------------------------
function y = logpost(x,infprob)    
    y = infbench_acerbidokka2018(x,infprob);
    lnp = infbench_lnprior(x,infprob);
    y = y + lnp;
end

%--------------------------------------------------------------------------
function ll = ll_deterministic(theta,data)

sigma_vis = exp(theta(1:3));
sigma_vest = exp(theta(4));
lambda = theta(5);
kappa = exp(theta(6));

ll = 0;

for iNoise = 1:3    
    X = data.X{iNoise};
    
    s_vest = X(:,3);
    s_vis = X(:,4);
    
    a_plus = (kappa + s_vest - s_vis)/sigma_vis(iNoise);
    a_minus = (-kappa + s_vest - s_vis)/sigma_vis(iNoise);
    b = sigma_vest/sigma_vis(iNoise);
    
    p_resp = normcdf(a_plus./sqrt(1+b.^2)) - normcdf(a_minus./sqrt(1+b.^2));
    p_resp = lambda/2 + (1-lambda)*p_resp;
    
    p = p_resp.*(X(:,5) == 1) + (1-p_resp).*(X(:,5) == 2);
    ll = ll + sum(log(p));    
end

end


%--------------------------------------------------------------------------
function [ll,ll_std] = ll_stochastic(theta,data)

sigma_vis = exp(theta(1:3));
sigma_vest = exp(theta(4));
lambda = theta(5);
kappa = exp(theta(6));

ll = 0; ll_var = 0;

for iNoise = 1:3    
    X = data.X{iNoise}(:,3:4);
    R = data.X{iNoise}(:,5);
        
    theta = [log(sigma_vis(iNoise)),log(sigma_vest),lambda,log(kappa)];
    
    options.Nreps = data.IBSNreps;
    fun = @(theta_,X_) gendata(theta_,X_);
    [nll,nlogl_var] = ibslike(fun,theta,R,X,options);
    
    ll = ll - nll;
    ll_var = ll_var + nlogl_var;
end

ll_std = sqrt(ll_var);

end

%--------------------------------------------------------------------------
function R = gendata(theta,X)

sigma_vis = exp(theta(1));
sigma_vest = exp(theta(2));
lambda = theta(3);
kappa = exp(theta(4));

Nt = size(X,1);

x_vis = randn(Nt,1)*sigma_vis + X(:,1);
x_vest = randn(Nt,1)*sigma_vest + X(:,2);

R = 2 - (abs(x_vis - x_vest) < kappa);

lapse_idx = rand(Nt,1) < lambda;
R(lapse_idx) = (rand(sum(lapse_idx),1) < 0.5) + 1;


end




