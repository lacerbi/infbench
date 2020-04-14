function [y,y_std] = infbench_akrami2018b(x,infprob,mcmc_params)
%INFBENCH_AKRAMI2018B Inference benchmark log pdf -- history-dependent model for mice data.

if nargin < 3; mcmc_params = []; end

problem_name = 'akrami2018b';
infbench_fun = str2func(['infbench_' problem_name]);


if isempty(x)
    if isempty(infprob) % Generate this document        
        fprintf('\n');

        % Add sampler directory to MATLAB path
        pathstr = fileparts(mfilename('fullpath'));
        addpath([pathstr,filesep(),'..',filesep(),'infalgo',filesep(),'parallel-GP-SL',filesep(),'utils',filesep(),'mcmcstat-master']);
                
        for n = 1     
            name = ['S' num2str(n)];
            
            infprob = infbench_fun([],n);
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
                
                opts = struct('Display','iter','MaxFunEvals',2e3);

                for iOpt = 1:Nopts                
                    x0 = rand(1,D).*(PUB - PLB) + PLB;
                    fun = @(x) -infbench_fun(x,infprob);
%                    [xnew(iOpt,:),fvalnew(iOpt)] = bads(fun,x0,LB,UB,PLB,PUB,[],opts);
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
%                [xnew,fvalnew] = bads(fun,x0,LB,UB,PLB,PUB,[],opts);
                [xnew,fvalnew] = fmincon(fun,x0,[],[],[],[],LB,UB,[],opts);

                fvalnew = -fvalnew;
                xmin_post = xmin;
                xmin_post = warpvars(xnew,'inv',trinfo);
                fval_post = fvalnew + warpvars(xnew,'logp',trinfo);

                fprintf('\t\t\tcase %d\n',n);
                fprintf('\t\t\t\tname = ''%s'';\n\t\t\t\txmin = %s;\n\t\t\t\tfval = %s;\n',name,mat2str(xmin),mat2str(fval));
                fprintf('\t\t\t\txmin_post = %s;\n\t\t\t\tfval_post = %s;\n',mat2str(xmin_post),mat2str(fval_post));
                
            elseif id > 0 && n == mcmc_params(2)

                rng(id);
                widths = 0.5*(PUB - PLB);
                logpfun = @(x) logpost(x,infprob);
                
                % Number of samples
                if numel(mcmc_params) > 2
                    Ns = mcmc_params(3);
                else
                    Ns = 1e3;
                end
                
                W = 2*(infprob.D+1);    % Number of walkers
                
                sampleopts.Thin = 17;
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
                
%                 % Re-evaluate final lls with higher precision
%                 fprintf('Re-evaluate log posteriors...\n');
%                 infprob.Ng = 501;
%                 lls = zeros(size(Xs,1),1);
%                 for i = 1:size(Xs,1)
%                     lls(i) = logpost(Xs(i,:),infprob);
%                     if mod(i,ceil(size(Xs,1)/10))==0; fprintf('%d/%d..',i,size(Xs,1)); end
%                 end
%                 fprintf('\nDone.\n');
                
                
                filename = [problem_name '_mcmc_n' num2str(n) '_id' num2str(id) '.mat'];
                save(filename,'Xs','lls','exitflag','output');                
            end
            
        end
        
    else
        % Initialization call -- define problem and set up data
        n = infprob(1);
        
        % Are we transforming the entire problem to unconstrained space?
        transform_to_unconstrained_coordinates = false;
        
        % Add problem directory to MATLAB path
%        pathstr = fileparts(mfilename('fullpath'));
%        addpath([pathstr,filesep(),problem_name]);
        
        % Define parameter upper/lower bounds
        Nr = 9;
        lb = [-3*ones(1,Nr)];
        ub = [3*ones(1,Nr)];
        plb = [-1*ones(1,Nr)];
        pub = [1*ones(1,Nr)];
        noise = [];
        
        D = numel(lb);
        xmin = NaN(1,D);       fval = Inf;
        xmin_post = NaN(1,D);  fval_post = Inf;
        Mean_mcmc = NaN(1,D);  Cov_mcmc = NaN(D,D);    lnZ_mcmc = NaN;
        
        if transform_to_unconstrained_coordinates
            trinfo = warpvars(D,lb,ub,plb,pub);     % Transform to unconstrained space
            trinfo.mu = zeros(1,D);     % Necessary for retro-compatibility
            trinfo.delta = ones(1,D);
        else
            trinfo = [];
        end
        
        switch n
			case {1,101}
                nid = 1;
				xmin = [0.825598323512873 -1.18015471271915 -0.0201916543320723 -2.7582078856892 2.99912100423627 0.197392890134053 -0.141297986743053 2.1721911172423 -0.816990721695071];
				fval = -6199.18585061038;
				xmin_post = [0.825594488578359 -1.18015106882503 -0.0201876987245857 -2.75460461939717 2.99551781922876 0.197387460895249 -0.141294051089077 2.169488644509 -0.816971359470748];
				fval_post = -6215.31250581914;
                % R_max = 1.000. Ntot = 100000. Neff_min = 93935.5. Total funccount = 10992360.
                Mean_mcmc = [0.601893215681988 0.929478095523648 -0.903846947948588 -1.25510425422428 0.532878582500004 -0.215684629220838];
                Cov_mcmc = [0.00177938856793166 -0.000907866511738251 -0.00147981117322351 0.000727299297923634 0.000423994956108903 -8.0698112932624e-05;-0.000907866511738251 0.00182625143109427 0.000720435936430774 -0.00141199493801528 -7.99342138443656e-05 -8.5564801234589e-05;-0.00147981117322351 0.000720435936430774 0.00238230151468306 -0.00118870882858125 -0.000324204738509761 7.21378524355696e-05;0.000727299297923634 -0.00141199493801528 -0.00118870882858125 0.0023687101586604 7.5678494290554e-05 6.38438035732974e-05;0.000423994956108903 -7.99342138443656e-05 -0.000324204738509761 7.5678494290554e-05 0.0010319947305803 -0.00049979629889294;-8.0698112932624e-05 -8.5564801234589e-05 7.21378524355696e-05 6.38438035732974e-05 -0.00049979629889294 0.000962327214581166];
                lnZ_mcmc = -12666.9953224803;                
        end

        % Read and preprocess data
        temp = csvread('akrami2018_data.csv',1);
        idx = size(temp,1)/2+1:size(temp,1);    % Take 2nd half of the data
        x = temp(idx,1:2);
        
        % Build matrix of regressors;
        x1 = circshift(x,1);
        x1(1,:) = 0;
        x2 = circshift(x,2);
        x2(1:2,:) = 0;
        r = double(x1(:,1) > x1(:,2));
        r(r == 0) = -1;
        xbar = mean(x1,2);
        tau = 20;
        ww = exp(-(1:200)/tau);
        ff = filter(ww/sum(ww),1,xbar);        
        data.X = [x, ones(size(x,1),1), x1, x2, r, ff];
        data.y = temp(idx,3);
        data.Ntrials = size(x,1);
                
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
                
        if n > 100
            data.IBSNreps = 0; % Deterministic problems            
        else
            data.IBSNreps = 100; % Reps used for IBS estimator
        end
        
        % Read marginals from file
        %marginals = load([problem_name '_marginals.mat']);
        %y.Post.MarginalBounds = marginals.MarginalBounds{nid};
        %y.Post.MarginalPdf = marginals.MarginalPdf{nid};
        
        % Save data and coordinate transformation struct
        data.trinfo = trinfo;
        y.Data = data;
        y.DeterministicFlag = (data.IBSNreps == 0);
                
    end
    
else
    
    % Iteration call -- evaluate objective function
    
    % Transform unconstrained variables to original space
    x_orig = warpvars(x,'i',infprob.Data.trinfo);
    dy = warpvars(x,'logpdf',infprob.Data.trinfo);   % Jacobian correction
    
    % Compute log likelihood of data and possibly std of log likelihood
    if infprob.DeterministicFlag
        LL = akrami2018_llfun(x_orig,infprob.Data);
        y_std = 0;
    else
        ibs_opts = struct('Nreps',infprob.Data.IBSNreps,...
            'ReturnPositive',true,'ReturnStd',true);
        [LL,y_std] = ibslike(@akrami2018_gendata,x_orig,infprob.Data.y,[],ibs_opts,infprob.Data);
    end
    y = LL + dy;
    
end

end

%--------------------------------------------------------------------------
function y = logpost(x,infprob)    
    y = infbench_akrami2018b(x,infprob);
    lnp = infbench_lnprior(x,infprob);
    y = y + lnp;
end

function ll = akrami2018_llfun(theta,data)

% Model parameters
w = theta(1:end);
lambda = 0.01;

p_left = 1./(1+exp(sum(bsxfun(@times,w,data.X),2)));
p_lresp = (data.y == 1).*p_left + (data.y == 2).*(1-p_left);
p_lresp = (1-lambda)*p_lresp + lambda/2;

ll = sum(log(p_lresp));

end



function R = akrami2018_gendata(theta,idx,data)

% Model parameters
w = theta(1:end);
lambda = 0.01;

p_left = 1./(1+exp(sum(bsxfun(@times,w,data.X),2)));
p_lresp = (1-lambda)*p_left + lambda/2;

pl_vec = p_lresp(idx);

R = double(rand(size(pl_vec)) < pl_vec);
R(R == 0) = 2;

end


function akrami2018_test(infprob)

theta = randn(1,9);
LL = akrami2018_llfun(theta,infprob.Data);
ibs_opts = struct('Nreps',infprob.Data.IBSNreps,...
    'ReturnPositive',true,'ReturnStd',true);
[LL_ibs,y_std] = ibslike(@akrami2018_gendata,theta,infprob.Data.y,[],ibs_opts,infprob.Data);

[LL,LL_ibs,y_std]
[LL - LL_ibs, y_std]

end
