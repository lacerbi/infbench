function [y,y_std] = infbench_goris2015b(x,infprob,mcmc_params)
%INFBENCH_GORIS2015 Inference benchmark log pdf -- neuronal model from Goris et al. (2015).

if nargin < 3; mcmc_params = []; end

problem_name = 'goris2015b';
infbench_fun = str2func(['infbench_' problem_name]);

if isempty(x)
    if isempty(infprob) % Generate this document        
        fprintf('\n');

%        for n = 7:12
        for n = 7:8
            switch n
                case 7;             name = 'm620r12';
                    xmin = [138.8505216444501 2.38317564009549 0.682321320237262 1.1613095607596 1 0.231748337632257 -0.272638945416596 3.10117864852662 72.8822298534178 0.00789002312857097 0.101380347749067 0.693895739234024];
                    fval = -2594.08310420223;
                case 8;             name = 'm620r35';
                    xmin = [227.530092052922 3.00555244729356 2.44308608399358 0.867188243443111 0.886173951591468 -0.039623616648953 -0.463062374186733 1.03689776623743 642.425443794774 0.001 8.65532672683422 0.0851310168322118];
                    fval = -6349.08642324277;
                case 9;             name = 'm624l54';
                    xmin = [219.1171252067330 1.39227769925336 3.03322024598022 0.617344863953975 0.355100712841601 0.067816661803776 -0.12988300841421 1.22284211173093 345.53039808358 0.001 4.29514870404523 0.308744430060413];
                    fval = -6497.67358979187;
                case 10;            name = 'm625r58';    
                    xmin = [305.847145521642 4.7748407736943 1.65588982473794 3.5 0.0559001728606333 0.0448467266543214 -0.192486743740829 4.81143259854959 349.659335938905 0.0700397166151502 0.213318145774832 1.05373239451993];
                    fval = -2118.86479506553;                
                case 11;            name = 'm625r62';    
                    xmin = [124.892643872064 3.69231469595064 2.86789754058724 0.364058707885479 0.00572852306727856 -0.529446086627358 0.410971814338923 4.63088807020964 59810475.5470703 0.001 0.001 0.173543268428183];
                    fval = -301.732736702732;
                case 12;            name = 'm627r58';    
                    xmin = [127.0183192213690 3.48851184430314 3.5 0.897489804982497 0.127319452033804 0.0551773381806372 -0.015460244000334 1.6645293750073 1068.80040238826 0.001 3.84129161213518 0.329767743116606];
                    fval = -5929.98517589428;                
            end                
            
            infprob = infbench_fun([],n);
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
            
                % First, check optimum
                opts = struct('Display','iter','MaxFunEvals',1e3);

                x0 = xmin(infprob.idxParams);
                x0 = warpvars(x0,'d',trinfo);
                fun = @(x) -infbench_fun(x,infprob);
                [xnew,fvalnew] = bads(fun,x0,LB,UB,PLB,PUB,[],opts);

                fvalnew = -fvalnew;
                xmin(infprob.idxParams) = warpvars(xnew,'inv',trinfo);
                fval = fvalnew - warpvars(xnew,'logp',trinfo);

                x0 = xmin(infprob.idxParams);
                x0 = warpvars(x0,'d',trinfo); 
                fun = @(x) -logpost(x,infprob);
                [xnew,fvalnew] = bads(fun,x0,LB,UB,PLB,PUB,[],opts);

                fvalnew = -fvalnew;
                xmin_post = xmin;
                xmin_post(infprob.idxParams) = warpvars(xnew,'inv',trinfo);
                fval_post = fvalnew - warpvars(xnew,'logp',trinfo);

                fprintf('\t\t\tcase %d\n',n);
                fprintf('\t\t\t\tname = ''%s'';\n\t\t\t\txmin = %s;\n\t\t\t\tfval = %s;\n',name,mat2str(xmin),mat2str(fval));
                fprintf('\t\t\t\txmin_post = %s;\n\t\t\t\tfval_post = %s;\n',mat2str(xmin_post),mat2str(fval_post));
                
            elseif id > 0 && n == mcmc_params(2)
                
                rng(id);
                widths = 0.5*(PUB - PLB);
                logpfun = @(x) logpost(x,infprob);
                
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
                
                [Xs,lls,exitflag,output] = eissample_lite(logpfun,x0,Ns,W,widths,LB,UB,sampleopts);
                
                filename = [problem_name '_mcmc_n' num2str(n) '_id' num2str(id) '.mat'];
                save(filename,'Xs','lls','exitflag','output');                
            end            
        end        
        
    else
        % Initialization call -- define problem and set up data
        n = infprob(1);
        
        % Are we transforming the entire problem to unconstrained space?
        transform_to_unconstrained_coordinates = false;        
        
        % The parameters and their bounds
        % 01 = preferred direction of motion (degrees), unbounded (periodic [0,360]), logical to use most effective stimulus value for family 1, high contrast as starting point
        % 02 = preferred spatial frequency (cycles per degree), values between [.05 15], logical to use most effective stimulus frequency as starting point
        % 03 = aspect ratio 2-D Gaussian, values between [.1 3.5], 1 is reasonable starting point
        % 04 = derivative order in space, values between [.1 3.5], 1 is reasonable starting point
        % 05 = directional selectivity, values between [0 1], 0.5 is reasonable starting point
        % 06 = gain inhibitory channel, values between [-1 1], but majority of cells between [-.2 .2], -0.1 is reasonable starting point
        % 07 = normalization constant, log10 basis, values between [-1 1], 0 is reasonable starting point 
        % 08 = response exponent, values between [1 6.5], 3 is reasonable starting point
        % 09 = response scalar, values between [1e-3 1e9], 4e3 is reasonable starting point (depending on choice of other starting points)
        % 10 = early additive noise, values between [1e-3 1e1], 0.1 is reasonable starting point
        % 11 = late additive noise, values between [1e-3 1e1], 0.1 is reasonable starting point
        % 12 = variance of response gain, values between [1e-3 1e1], 0.1 is reasonable starting point    

        % The most interesting parameters are 01–04, 06, 08, and 12
        
        % Datasets 1-6 are fake neurons
        % (used in parameter recovery analysis of Goris et al. 2015)    
        % Datasets 7-12 are real neurons (cells 9, 16, 45, 64, 65, and 78)
        
        xmin = NaN(1,12);       fval = Inf;
        xmin_post = NaN(1,12);  fval_post = Inf;
        Mean_laplace = NaN(1,7);    Cov_laplace = NaN(7,7); lnZ_laplace = NaN;
        Mean_mcmc = NaN(1,7);       Cov_mcmc = NaN(7,7);    lnZ_mcmc = NaN;
        
        switch n
            case {1,2,3,4,5,6}; name = ['fake0', int2str(n)];
            case 7                
                name = 'm620r12';
                xmin = [138.847253657877 2.38305768733247 0.682799619682061 1.1610646085556 1 0.231798049807549 -0.272638945416596 3.10144125856459 72.8822298534178 0.00789002312857097 0.101380347749067 0.69365902859132];
                fval = -2594.08312619286;
                xmin_post = [138.854185827076 2.383919716066 0.683555957652529 1.16017086317862 1 0.231799114868045 -0.272638945416596 3.10158656165004 72.8822298534178 0.00789002312857097 0.101380347749067 0.693926368767479];
                fval_post = -2609.82188588312;                                
                % R_max = 1.003. Ntot = 800000. Neff_min = 743680.7. Total funccount = 5172327.
                Mean_mcmc = [138.853267099868 2.42914044129109 0.615860859044795 1.22772040410999 0.232269065743085 3.10225296737271 0.710101903921601];
                Cov_mcmc = [1.67664293842744 0.000301496866863507 0.0108771297737337 -0.0105855160026286 0.000479461218380233 0.0021638085196806 -0.000456216056413771;0.000301496866863507 0.0322886348895838 0.00448910809268653 -0.00215753294025638 -5.68709372523392e-05 -0.000340438784217715 1.74551433236102e-05;0.0108771297737337 0.00448910809268653 0.0360050119739727 -0.0336130593285287 0.000273096312657433 0.00211200191145928 1.45437878520503e-05;-0.0105855160026286 -0.00215753294025638 -0.0336130593285287 0.0465405433220852 0.000155756318982554 -0.00308123432795199 -0.000179282353614662;0.000479461218380233 -5.68709372523392e-05 0.000273096312657433 0.000155756318982554 0.000225038258043983 0.00143399420907252 -1.48174894560478e-05;0.0021638085196806 -0.000340438784217715 0.00211200191145928 -0.00308123432795199 0.00143399420907252 0.0120291551779573 -0.00010795629972068;-0.000456216056413771 1.74551433236102e-05 1.45437878520503e-05 -0.000179282353614662 -1.48174894560478e-05 -0.00010795629972068 0.00385046783125344];
                lnZ_mcmc = -2618.97716613641;                
			case 8
                name = 'm620r35';
                xmin = [227.529653459492 3.00549354866858 2.44263708293926 0.867454206009722 0.886173951591468 -0.0396211358007847 -0.463062374186733 1.03691408842176 642.425443794774 0.001 8.65532672683422 0.085140255600458];
                fval = -6349.08645023851;
                xmin_post = [227.527988101867 3.00613478791754 2.44277087582304 0.867397258023728 0.886173951591468 -0.0395996648621448 -0.463062374186733 1.03699629821814 642.425443794774 0.001 8.65532672683422 0.0851209077401497];
                fval_post = -6364.82528308428;                    
                % R_max = 1.002. Ntot = 800000. Neff_min = 765771.8. Total funccount = 5212417.
                Mean_mcmc = [227.529689492625 2.97640807282226 2.43698656156055 0.872024378952079 -0.0394764566247722 1.03946854684131 0.0861299347560101];
                Cov_mcmc = [0.147717919544215 0.000510547054285867 -0.00048915268187323 0.000328648064162863 -1.15107407186054e-05 -9.38210642639461e-05 7.43598146479497e-06;0.000510547054285867 0.00776558741696669 -0.00114604344921399 0.000588583641993853 1.32303545253531e-05 -0.000280452441753885 -9.1355454464638e-06;-0.00048915268187323 -0.00114604344921399 0.00442811612829913 -0.00232863835101328 1.25461252355845e-06 -8.09452242583214e-05 1.91815801602183e-07;0.000328648064162863 0.000588583641993853 -0.00232863835101328 0.00163220188573308 1.98864578579399e-05 -6.00423939653006e-05 8.88387943966209e-07;-1.15107407186054e-05 1.32303545253531e-05 1.25461252355845e-06 1.98864578579399e-05 5.99313413954274e-06 1.41749025960364e-05 2.36407095638823e-07;-9.38210642639461e-05 -0.000280452441753885 -8.09452242583214e-05 -6.00423939653006e-05 1.41749025960364e-05 0.000181443737055837 1.82968682602167e-06;7.43598146479497e-06 -9.1355454464638e-06 1.91815801602183e-07 8.88387943966209e-07 2.36407095638823e-07 1.82968682602167e-06 2.7019565519849e-05];
                lnZ_mcmc = -6384.63215679431;
            case 9;             name = 'm624l54';
                xmin = [219.1171252067330 1.39227769925336 3.03322024598022 0.617344863953975 0.355100712841601 0.067816661803776 -0.12988300841421 1.22284211173093 345.53039808358 0.001 4.29514870404523 0.308744430060413];
                fval = -6497.67358979187;
            case 10;            name = 'm625r58';    
                xmin = [305.847145521642 4.7748407736943 1.65588982473794 3.5 0.0559001728606333 0.0448467266543214 -0.192486743740829 4.81143259854959 349.659335938905 0.0700397166151502 0.213318145774832 1.05373239451993];
                fval = -2118.86479506553;                
			case 11
				name = 'm625r62';
				xmin = [124.88710608225 3.68586576122773 2.96330657248823 0.328529880587874 0.00572852306727856 -0.547736645287558 0.410971814338923 4.56809152611831 59810475.5470703 0.001 0.001 0.129162534142113];
				fval = -302.033482141767;
				xmin_post = [124.865607867591 3.68556099680715 2.92014194714368 0.365291069910745 0.00572852306727856 -0.511617271295265 0.410971814338923 4.69464968545128 59810475.5470703 0.001 0.001 0.129337540160696];
				fval_post = -314.509056865964;
            case 12;            name = 'm627r58';    
                xmin = [127.0183192213690 3.48851184430314 3.5 0.897489804982497 0.127319452033804 0.0551773381806372 -0.015460244000334 1.6645293750073 1068.80040238826 0.001 3.84129161213518 0.329767743116606];
                fval = -5929.98517589428;                
        end

        % Add data directory to MATLAB path
        pathstr = fileparts(mfilename('fullpath'));
        addpath([pathstr,filesep(),'goris2015']);
        
        % Load the spatial frequency data
        temp = load([name, '_sf.mat']);
        data.sSF = temp.S;

        % Load the orientation mixture data
        temp = load([name, '_tp.mat']);
        data.sTP = temp.S;    
                
        % nvec = 1:12;
        lb = [0,.05,.1,.1, 0,-1,-1,1, 1e-3,1e-3,1e-3,1e-3];  % First variable is periodic
        ub = [360,15,3.5,3.5, 1,1,1,6.5, 1e9,1e1,1e1,1e1];
        plb =   [90,.5,.3,.3, 0,-0.3,-1,1.01, 10,1e-3,1e-2,1.5e-2];
        pub =   [270,10,3.2,3.2, 1,0.3,1,5, 1e5,1e-1,1,1];
        noise = [];
        
        % LL = -TPGiveBof(xmin, data.sTP, data.sSF);
        % [LL,fval]
        
        idx_params = [1:4,6,8,12];  % Only consider a subset of params
        D = numel(idx_params);
        
        if transform_to_unconstrained_coordinates
            trinfo = warpvars(D,lb(idx_params),ub(idx_params),plb(idx_params),pub(idx_params));     % Transform to unconstrained space
            trinfo.mu = zeros(1,D);     % Necessary for retro-compatibility
            trinfo.delta = ones(1,D);
        else
            trinfo = [];
        end
        y.xBaseFull = xmin;
        xmin = warpvars(xmin(idx_params),'d',trinfo);
        fval = fval + warpvars(xmin,'logp',trinfo);
        
        Mean = zeros(1,D);
        Cov = eye(D);
        Mode = xmin;
                
        y.D = D;
        y.LB = warpvars(lb(idx_params),'d',trinfo);
        y.UB = warpvars(ub(idx_params),'d',trinfo);
        y.PLB = warpvars(plb(idx_params),'d',trinfo);
        y.PUB = warpvars(pub(idx_params),'d',trinfo);
                
        y.lnZ = 0;              % Log normalization factor
        y.Mean = Mean;          % Distribution moments
        y.Cov = Cov;
        y.Mode = Mode;          % Mode of the pdf
        y.ModeFval = fval;
        y.idxParams = idx_params;
        
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
%        marginals = load([problem_name '_marginals.mat']);
%        y.Post.MarginalBounds = marginals.MarginalBounds{n};
%        y.Post.MarginalPdf = marginals.MarginalPdf{n};
        
        % Save data and coordinate transformation struct
        data.trinfo = trinfo;
        y.Data = data;        
    end
    
else
    
    % Iteration call -- evaluate objective function
    
    % Transform unconstrained variables to original space
    x_orig = warpvars(x,'i',infprob.Data.trinfo);
    dy = warpvars(x,'logpdf',infprob.Data.trinfo);   % Jacobian correction
    
    if all(isfinite(infprob.Post.xBaseFull))
        xfull = infprob.Post.xBaseFull;
    else
        xfull = infprob.xBaseFull;        
    end
    xfull(infprob.idxParams) = x_orig;
    
    % Compute log likelihood of data (fcn returns nLL)
    LL = -TPGiveBof(xfull, infprob.Data.sTP, infprob.Data.sSF);    
    y = LL + dy;
    y_std = 0;
    
end

end

%--------------------------------------------------------------------------
function y = logpost(x,infprob)
    y = infbench_goris2015(x,infprob);  % Get log-likelihood
    lnp = infbench_lnprior(x,infprob);
    y = y + lnp;
end