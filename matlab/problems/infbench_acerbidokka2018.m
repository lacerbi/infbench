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
                fval = fvalnew - warpvars(xnew,'logp',trinfo);

                x0 = xmin;
                x0 = warpvars(x0,'d',trinfo);   % Convert to unconstrained coordinates            
                fun = @(x) -logpost(x,infprob);
                [xnew,fvalnew] = fmincon(fun,x0,[],[],[],[],LB,UB,[],opts);

                fvalnew = -fvalnew;
                xmin_post = xmin;
                xmin_post = warpvars(xnew,'inv',trinfo);
                fval_post = fvalnew - warpvars(xnew,'logp',trinfo);

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
        lb = [0.5*ones(1,4),5e-3,0.25];
        ub = [80*ones(1,4),0.5,180];
        plb = [1*ones(1,4),0.01,1];
        pub = [40*ones(1,4),0.2,45];
        noise = [];
        
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
				name = 'S1';
				xmin = [6.20058541794124 5.88061641083698 10.1912429256554 4.04744822755874 0.0252735866538809 10.1385937945795];
				fval = -483.513343605127;
				xmin_post = [6.20058943029349 5.88061526563754 10.191246982267 4.04744742248115 0.0252736278818014 10.1385941974292];
				fval_post = -505.504741171978;
                % R_max = 1.001. Ntot = 100000. Neff_min = 54111.3. Total funccount = 11700125.
                Mean_mcmc = [4.43184948740104 4.03007356523433 9.44642934826877 5.37364365978457 0.0342797767874796 10.1956228010525];
                Cov_mcmc = [4.29368228671573 3.23503670174619 1.72234398523099 -3.18635536176589 -0.00103228727914319 0.0750531012021658;3.23503670174619 4.24229729932267 1.72434998521052 -3.2566406240311 -0.00116413766169158 0.0465711847985821;1.72234398523099 1.72434998521052 2.34875083144942 -1.65395411694875 -0.00202064165470951 0.128377443424245;-3.18635536176589 -3.2566406240311 -1.65395411694875 3.58517502135081 -0.00272662076011058 0.0885590382686063;-0.00103228727914319 -0.00116413766169158 -0.00202064165470951 -0.00272662076011058 0.000310341370240065 -0.0011075429545227;0.0750531012021658 0.0465711847985821 0.128377443424245 0.0885590382686063 -0.0011075429545227 0.194941301472202];
                lnZ_mcmc = -504.047571906284;
			case {2,102}
                nid = 2;
				name = 'S2';
				xmin = [4.34115830610531 13.6634127742158 21.9399171931226 4.17896878534463 0.143318687089021 19.9460973718226];
				fval = -429.775221706335;
				xmin_post = [4.3411532785444 13.6634151966651 21.9399235800863 4.17897270790704 0.143318663575363 19.9460973695369];
				fval_post = -451.766619273163;
                % R_max = 1.000. Ntot = 100000. Neff_min = 96835.7. Total funccount = 10708766.
                Mean_mcmc = [4.89282924167333 14.2064928026467 23.3964204349182 5.30162898591846 0.120236492440139 21.0238440637653];
                Cov_mcmc = [7.51391079682306 2.98799024398108 2.73964978057133 -2.85795575335986 -0.0555621351236203 1.39742004275659;2.98799024398108 4.65203375996802 1.94008801317162 -1.00299752638156 -0.0458445852902849 1.18037552150395;2.73964978057133 1.94008801317162 8.91409956549985 0.440125319579184 -0.0516383101162724 1.96085992736082;-2.85795575335986 -1.00299752638156 0.440125319579184 7.98982474137416 -0.0658279353535865 1.64768589362127;-0.0555621351236203 -0.0458445852902849 -0.0516383101162724 -0.0658279353535865 0.00279225363491406 -0.0412531711408089;1.39742004275659 1.18037552150395 1.96085992736082 1.64768589362127 -0.0412531711408089 1.8312021774307];
                lnZ_mcmc = -446.402055742672;			
            case {3,103}
                nid = 3;
				name = 'S3';
				xmin = [6.0285765440931 8.88084905314686 15.1350755736381 5.67520284353148 0.0867812417746608 17.9785316624426];
				fval = -583.757699194513;
				xmin_post = [6.0285793465732 8.88084921857473 15.1350785813087 5.67520375615201 0.0867811871268778 17.9785333080787];
				fval_post = -605.749096761348; 
                % R_max = 1.001. Ntot = 100000. Neff_min = 53804.7. Total funccount = 11321065.
                Mean_mcmc = [4.85120716206632 7.80217821705901 14.8144059031867 6.17441156930417 0.0939129918835423 18.1853036968548];
                Cov_mcmc = [6.77451126162769 4.36376451548747 2.40603225645798 -4.8793729506778 -0.0194611514544288 0.276629894500092;4.36376451548747 5.87076638649965 2.19701081134387 -4.03531156142191 -0.027544179287174 0.16082889768144;2.40603225645798 2.19701081134387 3.0501799401009 -1.72955121867392 -0.0266659581142849 0.387311053992236;-4.8793729506778 -4.03531156142191 -1.72955121867392 6.64866004871633 -0.030455803648185 0.480297016357256;-0.0194611514544288 -0.027544179287174 -0.0266659581142849 -0.030455803648185 0.0019910662050552 -0.0118101410151322;0.276629894500092 0.16082889768144 0.387311053992236 0.480297016357256 -0.0118101410151322 0.595241858536597];
                lnZ_mcmc = -601.865750332466;
        end

        temp = load('acerbidokka2018_data.mat');
        data.X = temp.unity_data{nid};  % Get subject's dataset
                
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
        marginals = load('acerbidokka2018_marginals.mat');
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
        LL = ll_deterministic(x_orig,infprob.Data);
        y_std = 0;
    else
        [LL,y_std] = ll_stochastic(x_orig,infprob.Data);
    end
    y = LL + dy;
    
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

sigma_vis = theta(1:3);
sigma_vest = theta(4);
lambda = theta(5);
kappa = theta(6);

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




