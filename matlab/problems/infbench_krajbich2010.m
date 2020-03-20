function [y,y_std] = infbench_krajbich2010(x,infprob,mcmc_params)
%INFBENCH_KRAJBICH2010 Inference benchmark log pdf -- attentional drift-diffusion model from Krajbich et al. (2010).

if nargin < 3; mcmc_params = []; end

problem_name = 'krajbich2010';
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
                Nopts = 5;
                
                opts = struct('Display','iter','MaxFunEvals',1e3);

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
                
                sampleopts.Thin = 1;
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
                
                % Re-evaluate final lls with higher precision
                fprintf('Re-evaluate log posteriors...\n');
                infprob.Ng = 501;
                lls = zeros(size(Xs,1),1);
                for i = 1:size(Xs,1)
                    lls(i) = logpost(Xs(i,:),infprob);
                    if mod(i,ceil(size(Xs,1)/10))==0; fprintf('%d/%d..',i,size(Xs,1)); end
                end
                fprintf('\nDone.\n');
                
                
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
                
        D = 4;        
        xmin = NaN(1,D);       fval = Inf;
        xmin_post = NaN(1,D);  fval_post = Inf;
        Mean_mcmc = NaN(1,D);  Cov_mcmc = NaN(D,D);    lnZ_mcmc = NaN;
                
        % Define parameter upper/lower bounds
        lb = [0.1,0,0,0.01];
        ub = [2,5,1,0.2];
        plb = [0.2,0.1,0.1,0.03];
        pub = [1,2,0.9,0.1];
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
                subjid = 13;
				xmin = [0.47659450173378 0.103994668275118 0.789939886331558 0.0100000001490116];
				fval = -495.323269935645;
				xmin_post = [0.476600289344788 0.103943380713463 0.790082359313965 0.0100000011920929];
				fval_post = -495.913835514614;
                % R_max = 1.010. Ntot = 100000. Neff_min = 19386.9. Total funccount = 410547.
                Mean_mcmc = [0.476566274201714 0.120251079262292 0.778912996540066 0.0809507440748478];
                Cov_mcmc = [0.000718594409091987 -5.59208743356943e-05 -9.41064361782147e-05 -5.06057908934977e-05;-5.59208743356943e-05 0.000433943538372732 -0.000422462643621677 0.000420874211259792;-9.41064361782147e-05 -0.000422462643621677 0.0135182146782027 0.000200355157564871;-5.06057908934977e-05 0.000420874211259792 0.000200355157564871 0.00222846889918172];
                lnZ_mcmc = -505.192475816479;                
			case 2
				name = 'S2';
                subjid = 16;
				xmin = [0.817608234286308 0.22548133940436 0.510444793105125 0.0155491371965036];
				fval = -369.441276164266;
				xmin_post = [0.817557442188263 0.225220758281648 0.510107111930847 0.0152974556013942];
				fval_post = -370.031856505713;
                % R_max = 1.010. Ntot = 100000. Neff_min = 24894.5. Total funccount = 451165.
                Mean_mcmc = [0.824023443246209 0.249580104836662 0.509408703416727 0.0727217640175234];
                Cov_mcmc = [0.00269561608998282 -9.89998716103414e-05 -0.00021165083906951 -0.00019250898136798;-9.89998716103414e-05 0.00143263537896748 0.00177998984444992 0.00076004377978344;-0.00021165083906951 0.00177998984444992 0.0114928921899632 0.000386780644677561;-0.00019250898136798 0.00076004377978344 0.000386780644677561 0.00183894691827946];
                lnZ_mcmc = -378.471263223713;
        end

        data = get_krajbich2010_data(subjid);
        data.R = [data.choice,data.totfixdurbin];
        data.Ng = 301;  % Grid size
                
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
        LL = krajbich2010_loglike(x_orig,infprob.Data);
        y_std = 0;
    else
        ibs_opts = struct('Nreps',infprob.Data.IBSNreps,...
            'ReturnPositive',true,'ReturnStd',true);
        [LL,y_std] = ibslike(@krajbich2010_gendata,x_orig,infprob.Data.R,[],ibs_opts,infprob.Data);
    end
    y = LL - dy;
    
end

end

%--------------------------------------------------------------------------
function y = logpost(x,infprob)    
    y = infbench_krajbich2010(x,infprob);
    lnp = infbench_lnprior(x,infprob);
    y = y + lnp;
end


