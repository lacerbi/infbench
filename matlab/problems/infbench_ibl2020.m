function [y,y_std] = infbench_ibl2020(x,infprob,mcmc_params)
%INFBENCH_IBL2020 Inference benchmark log pdf -- exponential averaging model for IBL data, IBL (2020).

if nargin < 3; mcmc_params = []; end

problem_name = 'ibl2020';
infbench_fun = str2func(['infbench_' problem_name]);


if isempty(x)
    if isempty(infprob) % Generate this document        
        fprintf('\n');

        % Add sampler directory to MATLAB path
        pathstr = fileparts(mfilename('fullpath'));
        addpath([pathstr,filesep(),'..',filesep(),'infalgo',filesep(),'parallel-GP-SL',filesep(),'utils',filesep(),'mcmcstat-master']);
                
        for n = 1:2            
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
                Nopts = 10;
                
                opts = struct('Display','iter','MaxFunEvals',2e3);

                for iOpt = 1:Nopts                
                    x0 = rand(1,D).*(PUB - PLB) + PLB;
                    fun = @(x) -infbench_fun(x,infprob);                    
                    [xnew(iOpt,:),fvalnew(iOpt)] = bads(fun,x0,LB,UB,PLB,PUB,[],opts);
                end
                
                [fvalnew,idx_best] = min(fvalnew);
                xnew = xnew(idx_best,:);
                
                fvalnew = -fvalnew;
                xmin = warpvars(xnew,'inv',trinfo);
                fval = fvalnew + warpvars(xnew,'logp',trinfo);

                x0 = xmin;
                x0 = warpvars(x0,'d',trinfo);   % Convert to unconstrained coordinates            
                fun = @(x) -logpost(x,infprob);
                [xnew,fvalnew] = bads(fun,x0,LB,UB,PLB,PUB,[],opts);

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
                sampleopts.Display = 'off';
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
        pathstr = fileparts(mfilename('fullpath'));
        addpath([pathstr,filesep(),problem_name]);

        model_name = 'exponential_contrastnoise_lapse';
        params = ibl2020_params_new(model_name);        
        
        % Define parameter upper/lower bounds
        bounds = ibl2020_setup_params([],params);
        lb = bounds.LB;
        ub = bounds.UB;
        plb = bounds.PLB;
        pub = bounds.PUB;
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
			case 1
				name = 'S1';    % CSHL_003 2019/08/16
				xmin = [-1.08594435100547 -1.80385774589036 1.46108739412787e-07 0.114367382666929 2.6352424018522e-05 1.21576715659691 2.53714885163401 2.24608793118832];
				fval = -199.424822353514;
				xmin_post = [-1.10461126363538 -1.81978994979406 9.52236470230794e-10 0.107049054730639 8.76854178238906e-06 1.26053476719053 2.54069789660959 2.25611330563372];
				fval_post = -209.721503262913;
                % R_max = 1.006. Ntot = 100000. Neff_min = 89595.7. Total funccount = 8396758.
                Mean_mcmc = [-1.36599353769662 -1.88960972718479 0.0541149000000237 0.110712038723506 0.0359808007844354 1.6748398057123 3.43106825333691 3.30107447162769];
                Cov_mcmc = [0.0502558292505181 0.0216278640554018 -0.00254778518307619 -7.75848335730976e-05 -0.0014119610394125 -0.0508255600849554 -0.0585779610097104 -0.0956270554699201;0.0216278640554018 0.10742117182685 -0.00237775520878403 -0.00168170904989171 -0.000877626601615884 -0.0809046687360803 -0.195924043610106 -0.134970199535634;-0.00254778518307619 -0.00237775520878403 0.00295271996583016 -0.000117087352124262 6.9754685366636e-06 0.00049432562494269 0.00119688396822461 -0.000858448310023431;-7.75848335730976e-05 -0.00168170904989171 -0.000117087352124262 0.00461886362504913 -0.000579389274078338 -0.000833219854770625 -0.0166362120976628 -0.00868947613304207;-0.0014119610394125 -0.000877626601615884 6.9754685366636e-06 -0.000579389274078338 0.00105031720397544 0.000576735958176922 -0.00283347115137788 -0.00221829576363554;-0.0508255600849554 -0.0809046687360803 0.00049432562494269 -0.000833219854770625 0.000576735958176922 0.646187570396557 0.735503726514117 0.75944822900212;-0.0585779610097104 -0.195924043610106 0.00119688396822461 -0.0166362120976628 -0.00283347115137788 0.735503726514117 2.89967920153727 2.55845926948532;-0.0956270554699201 -0.134970199535634 -0.000858448310023431 -0.00868947613304207 -0.00221829576363554 0.75944822900212 2.55845926948532 2.4220775128761];
                lnZ_mcmc = -219.580300204325;
			case 2
				name = 'S2';    % ZM_1091 2019/02/25
				xmin = [-1.30043055516692 -1.86697920048234 0.233118949406107 0.0466037712972337 1.45158595546424e-07 1.50053608647163 2.86863405395361 2.10350834309083];
				fval = -248.406336327939;
				xmin_post = [-1.30082224529663 -1.86699784751996 0.233914196902953 0.0468648614629464 3.03736175298608e-06 1.50343459329679 2.86987690154064 2.10571297693915];
				fval_post = -258.724751747531;
                % R_max = 1.007. Ntot = 100000. Neff_min = 95235.2. Total funccount = 9309791.
                Mean_mcmc = [-1.64494238124094 -1.80529022668457 0.174091297018239 0.0839583753937784 0.0492297942108852 1.63600687111949 2.28673412174734 1.70119539029634];
                Cov_mcmc = [0.243954601036117 0.000128299827557612 -0.0218439722634807 -0.010027438391431 -0.00289925061144382 -0.047892522650349 0.205166449859634 0.188854221290381;0.000128299827557612 0.0618399792700389 -0.00502610126236004 -0.000854430864828994 -0.000581248431864264 -0.0215563509890273 -0.0578284523249629 -0.00979280498911773;-0.0218439722634807 -0.00502610126236004 0.0140087364843992 0.000403986614405763 -0.0012655967898559 -0.000277574809870801 -0.0108723626099896 -0.0210085294605871;-0.010027438391431 -0.000854430864828994 0.000403986614405763 0.00432528003094709 -0.000498570937509632 -8.27020828035344e-05 -0.0105330980362888 -0.00317837867680053;-0.00289925061144382 -0.000581248431864264 -0.0012655967898559 -0.000498570937509632 0.00134695083183603 0.00146889174718617 -0.00597582782433392 -0.00660136810177476;-0.047892522650349 -0.0215563509890273 -0.000277574809870801 -8.27020828035344e-05 0.00146889174718617 0.111247582949372 0.0273218615965322 -0.00191929237464368;0.205166449859634 -0.0578284523249629 -0.0108723626099896 -0.0105330980362888 -0.00597582782433392 0.0273218615965322 0.453364885077011 0.341264140123874;0.188854221290381 -0.00979280498911773 -0.0210085294605871 -0.00317837867680053 -0.00660136810177476 -0.00191929237464368 0.341264140123874 0.329365093980618];
                lnZ_mcmc = -268.639328954035;                
        end

        temp = load('ibl2020_data.mat');        
        data = temp.data{n};
        data.params = params;
        
        %data.R = [data.choice,data.totfixdurbin];
        %data.Ng = 301;  % Grid size
                
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
                
        data.IBSNreps = 20; % Reps used for IBS estimator
        
        % Read marginals from file
        marginals = load([problem_name '_marginals.mat']);
        y.Post.MarginalBounds = marginals.MarginalBounds{n};
        y.Post.MarginalPdf = marginals.MarginalPdf{n};
        
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
        LL = -ibl2020_nllfun(x_orig,infprob.Data.params,infprob.Data);
        y_std = 0;
    else
        ibs_opts = struct('Nreps',infprob.Data.IBSNreps,...
            'ReturnPositive',true,'ReturnStd',true);
        [LL,y_std] = ibslike(@krajbich2010_gendata,x_orig,infprob.Data.R,[],ibs_opts,infprob.Data);
    end
    y = LL + dy;
    
end

end

%--------------------------------------------------------------------------
function y = logpost(x,infprob)    
    y = infbench_ibl2020(x,infprob);
    lnp = infbench_lnprior(x,infprob);
    y = y + lnp;
end


