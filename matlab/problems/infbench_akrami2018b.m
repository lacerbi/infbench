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
                % R_max = 1.001. Ntot = 100000. Neff_min = 94032.4. Total funccount = 10209971.
                Mean_mcmc = [0.826097271561087 -1.18136471604667 -0.0189569331821486 -1.22888526306926 1.47028139531814 0.197472431871508 -0.14116344785186 1.02597956363325 -0.830105558788376];
                Cov_mcmc = [0.000935301680140307 -0.000735509765742965 -3.3482210729806e-05 -0.000204883120555423 0.000261022702666796 -3.94057589886134e-05 4.84451793519011e-05 0.000127998371046209 -0.000525167408161282;-0.000735509765742965 0.00122635381851453 3.12079532007797e-05 9.50297438257228e-06 -7.47126528160998e-05 2.75790512646124e-05 -3.95284030743019e-05 1.40903166738186e-05 0.00040367665177839;-3.3482210729806e-05 3.12079532007797e-05 0.000491220540566427 0.00128476198405997 -0.00127084688755487 -2.42429446088502e-05 3.83611559449289e-05 -0.000982106810434918 -0.000352103988836941;-0.000204883120555423 9.50297438257228e-06 0.00128476198405997 1.3778469605393 -1.37708207903235 -0.000283564866376109 0.00064538437983087 -1.03249681398292 -0.0126033684642543;0.000261022702666796 -7.47126528160998e-05 -0.00127084688755487 -1.37708207903235 1.37706406754047 0.00034602618097123 -0.000634606726068061 1.03216437692207 0.0112107226743759;-3.94057589886134e-05 2.75790512646124e-05 -2.42429446088502e-05 -0.000283564866376109 0.00034602618097123 0.0010121654633564 -0.000692525646422178 0.000188593203055965 -0.000976065572706183;4.84451793519011e-05 -3.95284030743019e-05 3.83611559449289e-05 0.00064538437983087 -0.000634606726068061 -0.000692525646422178 0.00111293837705694 -0.000428861407809761 -0.000310087609361833;0.000127998371046209 1.40903166738186e-05 -0.000982106810434918 -1.03249681398292 1.03216437692207 0.000188593203055965 -0.000428861407809761 0.774276683565385 0.00872618378972165;-0.000525167408161282 0.00040367665177839 -0.000352103988836941 -0.0126033684642543 0.0112107226743759 -0.000976065572706183 -0.000310087609361833 0.00872618378972165 0.0242403838424609];
                lnZ_mcmc = -6233.89696831064;
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
            data.IBSNreps = 200; % Reps used for IBS estimator
        end
        
        % Read marginals from file
        marginals = load([problem_name '_marginals.mat']);
        y.Post.MarginalBounds = marginals.MarginalBounds{nid};
        y.Post.MarginalPdf = marginals.MarginalPdf{nid};
        
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
