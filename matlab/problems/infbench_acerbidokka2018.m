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
                % R_max = 1.029. Ntot = 100000. Neff_min = 2378.9. Total funccount = 12814064.
                Mean_mcmc = [0.64709236994056 0.514658951518512 2.12154489538113 1.83348811941973 0.0354736708543006 2.31391349100702];
                Cov_mcmc = [0.56713183738407 0.242392190065536 0.045927450937598 -0.142445929217476 3.58557926957e-05 0.00108220655383944;0.242392190065536 0.539744526630762 0.044368971935126 -0.147091250980619 -2.34637144437161e-05 0.000337475183706663;0.045927450937598 0.044368971935126 0.0369346638298177 -0.0267583552612944 -0.000199827588828711 0.00121988609132642;-0.142445929217476 -0.147091250980619 -0.0267583552612944 0.133333928653792 -0.000509050097385322 0.00162501507235783;3.58557926957e-05 -2.34637144437161e-05 -0.000199827588828711 -0.000509050097385322 0.000326857955061985 -0.000118240894872711;0.00108220655383944 0.000337475183706663 0.00121988609132642 0.00162501507235783 -0.000118240894872711 0.00186067501230593];
            	lnZ_mcmc = -497.801320392954;    
			case 2
				name = 'S2';
				xmin = [1.38726609391213 2.60712631868221 3.08537558206964 1.5047752731606 0.143318676649389 2.99303342142661];
				fval = -429.775221706332;
				xmin_post = [1.38726676853333 2.60712639983091 3.08537554925837 1.50477233430562 0.143318797696291 2.99303327788359];
				fval_post = -437.453388225819;
                % R_max = 1.004. Ntot = 100000. Neff_min = 53660.6. Total funccount = 13661933.
                Mean_mcmc = [0.931873292033045 2.58648306039213 3.0997837475777 1.10970984492108 0.148003555532701 3.01027238055142];
                Cov_mcmc = [0.697224319433399 0.0425421798412317 0.019237179446457 -0.27345692874211 -0.00764326741626143 0.00713132959421961;0.0425421798412317 0.0285939813950335 0.00437634659151812 -0.0293384165618791 -0.00266321084791802 0.00216929505870157;0.019237179446457 0.00437634659151812 0.0153065351763119 -0.00345968769245005 -0.00148485979355223 0.00263119590710567;-0.27345692874211 -0.0293384165618791 -0.00345968769245005 0.712345288929726 -0.0119224169158888 0.011577946568483;-0.00764326741626143 -0.00266321084791802 -0.00148485979355223 -0.0119224169158888 0.00223491428849746 -0.0011350759484854;0.00713132959421961 0.00216929505870157 0.00263119590710567 0.011577946568483 -0.0011350759484854 0.00283418783181604];
                lnZ_mcmc = -443.292795250245;
			case 3
				name = 'S3';
				xmin = [1.19318104591285 1.98890187717717 2.6580800074951 2.02743456496656 0.0867816136434034 2.88917816646634];
				fval = -583.757699194538;
				xmin_post = [1.19394210042755 1.98905711125483 2.65812076676506 2.02729160630505 0.086781453413982 2.88917822194199];
				fval_post = -591.435865713995;                
                % R_max = 1.013. Ntot = 100000. Neff_min = 6145.3. Total funccount = 15198976.
                Mean_mcmc = [0.746194728589647 1.55007922838981 2.62598381048443 1.78075823279715 0.110284296973345 2.8908162572388];
                Cov_mcmc = [0.748559375960408 0.300296257197953 0.0443650083011604 -0.345461264682485 -0.00198630944568612 -0.000427887250993565;0.300296257197953 0.584580972844009 0.0384792232008505 -0.237314793019631 -0.00451627661372513 -0.00218169251768176;0.0443650083011604 0.0384792232008505 0.0165559465525891 -0.0310890367122426 -0.00175205560771741 0.00123511770927753;-0.345461264682485 -0.237314793019631 -0.0310890367122426 0.451674033813764 -0.0075858421994749 0.00489091229884737;-0.00198630944568612 -0.00451627661372513 -0.00175205560771741 -0.0075858421994749 0.00228594914890072 -0.000618052889406043;-0.000427887250993565 -0.00218169251768176 0.00123511770927753 0.00489091229884737 -0.000618052889406043 0.00180859559758757];
                lnZ_mcmc = -597.623155164008;
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




