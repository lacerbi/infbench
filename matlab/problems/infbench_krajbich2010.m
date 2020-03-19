function [y,y_std] = infbench_krajbich2010(x,infprob,mcmc_params)
%INFBENCH_KRAJBICH2010 Inference benchmark log pdf -- attentional drift-diffusion model from Krajbich et al. (2010).

if nargin < 3; mcmc_params = []; end

if isempty(x)
    if isempty(infprob) % Generate this document        
        fprintf('\n');

        % Add sampler directory to MATLAB path
        pathstr = fileparts(mfilename('fullpath'));
        addpath([pathstr,filesep(),'..',filesep(),'infalgo',filesep(),'parallel-GP-SL',filesep(),'utils',filesep(),'mcmcstat-master']);
                
        for n = 1:2            
            name = ['S' num2str(n)];
            
            infprob = infbench_krajbich2010([],n);
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
                    fun = @(x) -infbench_krajbich2010(x,infprob);                    
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
                
                
                filename = ['krajbich2010_mcmc_n' num2str(n) '_id' num2str(id) '.mat'];
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
        addpath([pathstr,filesep(),'krajbich2010']);
                
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
                % R_max = 1.001. Ntot = 100000. Neff_min = 54111.3. Total funccount = 11700125.
                Mean_mcmc = [4.43184948740104 4.03007356523433 9.44642934826877 5.37364365978457 0.0342797767874796 10.1956228010525];
                Cov_mcmc = [4.29368228671573 3.23503670174619 1.72234398523099 -3.18635536176589 -0.00103228727914319 0.0750531012021658;3.23503670174619 4.24229729932267 1.72434998521052 -3.2566406240311 -0.00116413766169158 0.0465711847985821;1.72234398523099 1.72434998521052 2.34875083144942 -1.65395411694875 -0.00202064165470951 0.128377443424245;-3.18635536176589 -3.2566406240311 -1.65395411694875 3.58517502135081 -0.00272662076011058 0.0885590382686063;-0.00103228727914319 -0.00116413766169158 -0.00202064165470951 -0.00272662076011058 0.000310341370240065 -0.0011075429545227;0.0750531012021658 0.0465711847985821 0.128377443424245 0.0885590382686063 -0.0011075429545227 0.194941301472202];
                lnZ_mcmc = -504.047571906284;
			case 2
				name = 'S2';
                subjid = 16;
				xmin = [0.817608234286308 0.22548133940436 0.510444793105125 0.0155491371965036];
				fval = -369.441276164266;
				xmin_post = [0.817557442188263 0.225220758281648 0.510107111930847 0.0152974556013942];
				fval_post = -370.031856505713;                
                % R_max = 1.000. Ntot = 100000. Neff_min = 96835.7. Total funccount = 10708766.
                Mean_mcmc = [4.89282924167333 14.2064928026467 23.3964204349182 5.30162898591846 0.120236492440139 21.0238440637653];
                Cov_mcmc = [7.51391079682306 2.98799024398108 2.73964978057133 -2.85795575335986 -0.0555621351236203 1.39742004275659;2.98799024398108 4.65203375996802 1.94008801317162 -1.00299752638156 -0.0458445852902849 1.18037552150395;2.73964978057133 1.94008801317162 8.91409956549985 0.440125319579184 -0.0516383101162724 1.96085992736082;-2.85795575335986 -1.00299752638156 0.440125319579184 7.98982474137416 -0.0658279353535865 1.64768589362127;-0.0555621351236203 -0.0458445852902849 -0.0516383101162724 -0.0658279353535865 0.00279225363491406 -0.0412531711408089;1.39742004275659 1.18037552150395 1.96085992736082 1.64768589362127 -0.0412531711408089 1.8312021774307];
                lnZ_mcmc = -446.402055742672;			
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
                
        data.IBSNreps = 200; % Reps used for IBS estimator
        
        % Read marginals from file
        %marginals = load('acerbidokka2018_marginals.mat');
        %y.Post.MarginalBounds = marginals.MarginalBounds{n};
        %y.Post.MarginalPdf = marginals.MarginalPdf{n};
        
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



%--------------------------------------------------------------------------
function [ll,ll_std] = ll_stochastic(theta,data)

sigma_vis = theta(1:3);
sigma_vest = theta(4);
lambda = theta(5);
kappa = theta(6);

ll = 0; ll_var = 0;

for iNoise = 1:3    
    X = data.X{iNoise}(:,3:4);
    R = data.X{iNoise}(:,5);
        
    theta = [sigma_vis(iNoise),sigma_vest,lambda,kappa];
    
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

sigma_vis = theta(1);
sigma_vest = theta(2);
lambda = theta(3);
kappa = theta(4);

Nt = size(X,1);

x_vis = randn(Nt,1)*sigma_vis + X(:,1);
x_vest = randn(Nt,1)*sigma_vest + X(:,2);

R = 2 - (abs(x_vis - x_vest) < kappa);

lapse_idx = rand(Nt,1) < lambda;
R(lapse_idx) = (rand(sum(lapse_idx),1) < 0.5) + 1;


end




