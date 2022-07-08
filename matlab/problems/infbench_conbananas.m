function [y,y_std] = infbench_conbananas(x,infprob,mcmc_params)
%INFBENCH_CONBANANAS Inference benchmark log pdf -- constrained product of Rosenbrock functions.

if nargin < 3; mcmc_params = []; end

lb = 0;
ub = 10;

if isempty(x)
    if isempty(infprob) % Generate this document        
                
        % Compute normalization constant (note that this also depends on the bounds)
        Z = integral2(@(x1,x2) reshape(exp(ll_rosenbrock([x1(:),x2(:)])),size(x1)),lb,ub,lb,ub);
        
        % Compute moments
        m1 = integral2(@(x1,x2) reshape(x1(:).*exp(ll_rosenbrock([x1(:),x2(:)])),size(x1)),lb,ub,lb,ub)/Z;
        m2 = integral2(@(x1,x2) reshape(x2(:).*exp(ll_rosenbrock([x1(:),x2(:)])),size(x1)),lb,ub,lb,ub)/Z;
        c1 = integral2(@(x1,x2) reshape((x1(:)-m1).^2.*exp(ll_rosenbrock([x1(:),x2(:)])),size(x1)),lb,ub,lb,ub)/Z;
        c2 = integral2(@(x1,x2) reshape((x2(:)-m2).^2.*exp(ll_rosenbrock([x1(:),x2(:)])),size(x1)),lb,ub,lb,ub)/Z;
        c12 = integral2(@(x1,x2) reshape((x1(:)-m1).*(x2(:)-m2).*exp(ll_rosenbrock([x1(:),x2(:)])),size(x1)),lb,ub,lb,ub)/Z;        
        Mu2 = [m1,m2]
        Sigma2 = [c1 c12; c12 c2]

        Z = Z / (ub - lb)^2;    % Renormalize for uniform prior in 2D
        
        % Compute mode
        x0 = ones(1,2);
        x_mode2 = fmincon(@(x) -ll_rosenbrock(x),x0,[],[],[],[],lb*[1,1],ub*[1,1])
        
        fprintf('switch D\n');        
        for D = 2:2:20
            Mu = zeros(1,D);
            Sigma = zeros(D,D);
            for d = 1:D/2
                Mu(d*2-1:d*2) = Mu2;
                Sigma(d*2-1:d*2,d*2-1:d*2) = Sigma2;                
            end
            Mode = repmat(x_mode2,[1,D/2]);
            fprintf('\tcase %d\n',D);
            fprintf('\t\tMean = %s;\n',mat2str(Mu));
            fprintf('\t\tCov = %s;\n',mat2str(Sigma));
            fprintf('\t\tlnZ = %s;\n',mat2str(log(Z)*(D/2)));
            fprintf('\t\tMode = %s;\n',mat2str(Mode));
        end
        fprintf('\totherwise\n\t\terror(''%s:TooHighDim'',''Benchmark function supports up to D=20 dimensions.'');\n',mfilename);
        fprintf('end\n');
    else    
        % Initialization call -- define problem and set up data
        D = infprob(1);

        switch D
            case 2
                Mean = [1.1026780672163 1.68840019758032];
                Cov = [0.383121902786579 0.755785036140133;0.755785036140133 2.10025209165182];
                lnZ = -7.75728877736958;
                Mode = [0.147411032194048 0.0205881504343637];
            case 4
                Mean = [1.1026780672163 1.68840019758032 1.1026780672163 1.68840019758032];
                Cov = [0.383121902786579 0.755785036140133 0 0;0.755785036140133 2.10025209165182 0 0;0 0 0.383121902786579 0.755785036140133;0 0 0.755785036140133 2.10025209165182];
                lnZ = -15.5145775547392;
                Mode = [0.147411032194048 0.0205881504343637 0.147411032194048 0.0205881504343637];
            case 6
                Mean = [1.1026780672163 1.68840019758032 1.1026780672163 1.68840019758032 1.1026780672163 1.68840019758032];
                Cov = [0.383121902786579 0.755785036140133 0 0 0 0;0.755785036140133 2.10025209165182 0 0 0 0;0 0 0.383121902786579 0.755785036140133 0 0;0 0 0.755785036140133 2.10025209165182 0 0;0 0 0 0 0.383121902786579 0.755785036140133;0 0 0 0 0.755785036140133 2.10025209165182];
                lnZ = -23.2718663321087;
                Mode = [0.147411032194048 0.0205881504343637 0.147411032194048 0.0205881504343637 0.147411032194048 0.0205881504343637];
            case 8
                Mean = [1.1026780672163 1.68840019758032 1.1026780672163 1.68840019758032 1.1026780672163 1.68840019758032 1.1026780672163 1.68840019758032];
                Cov = [0.383121902786579 0.755785036140133 0 0 0 0 0 0;0.755785036140133 2.10025209165182 0 0 0 0 0 0;0 0 0.383121902786579 0.755785036140133 0 0 0 0;0 0 0.755785036140133 2.10025209165182 0 0 0 0;0 0 0 0 0.383121902786579 0.755785036140133 0 0;0 0 0 0 0.755785036140133 2.10025209165182 0 0;0 0 0 0 0 0 0.383121902786579 0.755785036140133;0 0 0 0 0 0 0.755785036140133 2.10025209165182];
                lnZ = -31.0291551094783;
                Mode = [0.147411032194048 0.0205881504343637 0.147411032194048 0.0205881504343637 0.147411032194048 0.0205881504343637 0.147411032194048 0.0205881504343637];
            case 10
                Mean = [1.1026780672163 1.68840019758032 1.1026780672163 1.68840019758032 1.1026780672163 1.68840019758032 1.1026780672163 1.68840019758032 1.1026780672163 1.68840019758032];
                Cov = [0.383121902786579 0.755785036140133 0 0 0 0 0 0 0 0;0.755785036140133 2.10025209165182 0 0 0 0 0 0 0 0;0 0 0.383121902786579 0.755785036140133 0 0 0 0 0 0;0 0 0.755785036140133 2.10025209165182 0 0 0 0 0 0;0 0 0 0 0.383121902786579 0.755785036140133 0 0 0 0;0 0 0 0 0.755785036140133 2.10025209165182 0 0 0 0;0 0 0 0 0 0 0.383121902786579 0.755785036140133 0 0;0 0 0 0 0 0 0.755785036140133 2.10025209165182 0 0;0 0 0 0 0 0 0 0 0.383121902786579 0.755785036140133;0 0 0 0 0 0 0 0 0.755785036140133 2.10025209165182];
                lnZ = -38.7864438868479;
                Mode = [0.147411032194048 0.0205881504343637 0.147411032194048 0.0205881504343637 0.147411032194048 0.0205881504343637 0.147411032194048 0.0205881504343637 0.147411032194048 0.0205881504343637];
            case 12
                Mean = [1.1026780672163 1.68840019758032 1.1026780672163 1.68840019758032 1.1026780672163 1.68840019758032 1.1026780672163 1.68840019758032 1.1026780672163 1.68840019758032 1.1026780672163 1.68840019758032];
                Cov = [0.383121902786579 0.755785036140133 0 0 0 0 0 0 0 0 0 0;0.755785036140133 2.10025209165182 0 0 0 0 0 0 0 0 0 0;0 0 0.383121902786579 0.755785036140133 0 0 0 0 0 0 0 0;0 0 0.755785036140133 2.10025209165182 0 0 0 0 0 0 0 0;0 0 0 0 0.383121902786579 0.755785036140133 0 0 0 0 0 0;0 0 0 0 0.755785036140133 2.10025209165182 0 0 0 0 0 0;0 0 0 0 0 0 0.383121902786579 0.755785036140133 0 0 0 0;0 0 0 0 0 0 0.755785036140133 2.10025209165182 0 0 0 0;0 0 0 0 0 0 0 0 0.383121902786579 0.755785036140133 0 0;0 0 0 0 0 0 0 0 0.755785036140133 2.10025209165182 0 0;0 0 0 0 0 0 0 0 0 0 0.383121902786579 0.755785036140133;0 0 0 0 0 0 0 0 0 0 0.755785036140133 2.10025209165182];
                lnZ = -46.5437326642175;
                Mode = [0.147411032194048 0.0205881504343637 0.147411032194048 0.0205881504343637 0.147411032194048 0.0205881504343637 0.147411032194048 0.0205881504343637 0.147411032194048 0.0205881504343637 0.147411032194048 0.0205881504343637];
            case 14
                Mean = [1.1026780672163 1.68840019758032 1.1026780672163 1.68840019758032 1.1026780672163 1.68840019758032 1.1026780672163 1.68840019758032 1.1026780672163 1.68840019758032 1.1026780672163 1.68840019758032 1.1026780672163 1.68840019758032];
                Cov = [0.383121902786579 0.755785036140133 0 0 0 0 0 0 0 0 0 0 0 0;0.755785036140133 2.10025209165182 0 0 0 0 0 0 0 0 0 0 0 0;0 0 0.383121902786579 0.755785036140133 0 0 0 0 0 0 0 0 0 0;0 0 0.755785036140133 2.10025209165182 0 0 0 0 0 0 0 0 0 0;0 0 0 0 0.383121902786579 0.755785036140133 0 0 0 0 0 0 0 0;0 0 0 0 0.755785036140133 2.10025209165182 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0.383121902786579 0.755785036140133 0 0 0 0 0 0;0 0 0 0 0 0 0.755785036140133 2.10025209165182 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0.383121902786579 0.755785036140133 0 0 0 0;0 0 0 0 0 0 0 0 0.755785036140133 2.10025209165182 0 0 0 0;0 0 0 0 0 0 0 0 0 0 0.383121902786579 0.755785036140133 0 0;0 0 0 0 0 0 0 0 0 0 0.755785036140133 2.10025209165182 0 0;0 0 0 0 0 0 0 0 0 0 0 0 0.383121902786579 0.755785036140133;0 0 0 0 0 0 0 0 0 0 0 0 0.755785036140133 2.10025209165182];
                lnZ = -54.301021441587;
                Mode = [0.147411032194048 0.0205881504343637 0.147411032194048 0.0205881504343637 0.147411032194048 0.0205881504343637 0.147411032194048 0.0205881504343637 0.147411032194048 0.0205881504343637 0.147411032194048 0.0205881504343637 0.147411032194048 0.0205881504343637];
            case 16
                Mean = [1.1026780672163 1.68840019758032 1.1026780672163 1.68840019758032 1.1026780672163 1.68840019758032 1.1026780672163 1.68840019758032 1.1026780672163 1.68840019758032 1.1026780672163 1.68840019758032 1.1026780672163 1.68840019758032 1.1026780672163 1.68840019758032];
                Cov = [0.383121902786579 0.755785036140133 0 0 0 0 0 0 0 0 0 0 0 0 0 0;0.755785036140133 2.10025209165182 0 0 0 0 0 0 0 0 0 0 0 0 0 0;0 0 0.383121902786579 0.755785036140133 0 0 0 0 0 0 0 0 0 0 0 0;0 0 0.755785036140133 2.10025209165182 0 0 0 0 0 0 0 0 0 0 0 0;0 0 0 0 0.383121902786579 0.755785036140133 0 0 0 0 0 0 0 0 0 0;0 0 0 0 0.755785036140133 2.10025209165182 0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0.383121902786579 0.755785036140133 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0.755785036140133 2.10025209165182 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0.383121902786579 0.755785036140133 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0.755785036140133 2.10025209165182 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0 0.383121902786579 0.755785036140133 0 0 0 0;0 0 0 0 0 0 0 0 0 0 0.755785036140133 2.10025209165182 0 0 0 0;0 0 0 0 0 0 0 0 0 0 0 0 0.383121902786579 0.755785036140133 0 0;0 0 0 0 0 0 0 0 0 0 0 0 0.755785036140133 2.10025209165182 0 0;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.383121902786579 0.755785036140133;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.755785036140133 2.10025209165182];
                lnZ = -62.0583102189566;
                Mode = [0.147411032194048 0.0205881504343637 0.147411032194048 0.0205881504343637 0.147411032194048 0.0205881504343637 0.147411032194048 0.0205881504343637 0.147411032194048 0.0205881504343637 0.147411032194048 0.0205881504343637 0.147411032194048 0.0205881504343637 0.147411032194048 0.0205881504343637];
            case 18
                Mean = [1.1026780672163 1.68840019758032 1.1026780672163 1.68840019758032 1.1026780672163 1.68840019758032 1.1026780672163 1.68840019758032 1.1026780672163 1.68840019758032 1.1026780672163 1.68840019758032 1.1026780672163 1.68840019758032 1.1026780672163 1.68840019758032 1.1026780672163 1.68840019758032];
                Cov = [0.383121902786579 0.755785036140133 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;0.755785036140133 2.10025209165182 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;0 0 0.383121902786579 0.755785036140133 0 0 0 0 0 0 0 0 0 0 0 0 0 0;0 0 0.755785036140133 2.10025209165182 0 0 0 0 0 0 0 0 0 0 0 0 0 0;0 0 0 0 0.383121902786579 0.755785036140133 0 0 0 0 0 0 0 0 0 0 0 0;0 0 0 0 0.755785036140133 2.10025209165182 0 0 0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0.383121902786579 0.755785036140133 0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0.755785036140133 2.10025209165182 0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0.383121902786579 0.755785036140133 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0.755785036140133 2.10025209165182 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0 0.383121902786579 0.755785036140133 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0 0.755785036140133 2.10025209165182 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0 0 0 0.383121902786579 0.755785036140133 0 0 0 0;0 0 0 0 0 0 0 0 0 0 0 0 0.755785036140133 2.10025209165182 0 0 0 0;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.383121902786579 0.755785036140133 0 0;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.755785036140133 2.10025209165182 0 0;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.383121902786579 0.755785036140133;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.755785036140133 2.10025209165182];
                lnZ = -69.8155989963262;
                Mode = [0.147411032194048 0.0205881504343637 0.147411032194048 0.0205881504343637 0.147411032194048 0.0205881504343637 0.147411032194048 0.0205881504343637 0.147411032194048 0.0205881504343637 0.147411032194048 0.0205881504343637 0.147411032194048 0.0205881504343637 0.147411032194048 0.0205881504343637 0.147411032194048 0.0205881504343637];
            case 20
                Mean = [1.1026780672163 1.68840019758032 1.1026780672163 1.68840019758032 1.1026780672163 1.68840019758032 1.1026780672163 1.68840019758032 1.1026780672163 1.68840019758032 1.1026780672163 1.68840019758032 1.1026780672163 1.68840019758032 1.1026780672163 1.68840019758032 1.1026780672163 1.68840019758032 1.1026780672163 1.68840019758032];
                Cov = [0.383121902786579 0.755785036140133 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;0.755785036140133 2.10025209165182 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;0 0 0.383121902786579 0.755785036140133 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;0 0 0.755785036140133 2.10025209165182 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;0 0 0 0 0.383121902786579 0.755785036140133 0 0 0 0 0 0 0 0 0 0 0 0 0 0;0 0 0 0 0.755785036140133 2.10025209165182 0 0 0 0 0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0.383121902786579 0.755785036140133 0 0 0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0.755785036140133 2.10025209165182 0 0 0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0.383121902786579 0.755785036140133 0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0.755785036140133 2.10025209165182 0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0 0.383121902786579 0.755785036140133 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0 0.755785036140133 2.10025209165182 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0 0 0 0.383121902786579 0.755785036140133 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0 0 0 0.755785036140133 2.10025209165182 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.383121902786579 0.755785036140133 0 0 0 0;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.755785036140133 2.10025209165182 0 0 0 0;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.383121902786579 0.755785036140133 0 0;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.755785036140133 2.10025209165182 0 0;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.383121902786579 0.755785036140133;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.755785036140133 2.10025209165182];
                lnZ = -77.5728877736958;
                Mode = [0.147411032194048 0.0205881504343637 0.147411032194048 0.0205881504343637 0.147411032194048 0.0205881504343637 0.147411032194048 0.0205881504343637 0.147411032194048 0.0205881504343637 0.147411032194048 0.0205881504343637 0.147411032194048 0.0205881504343637 0.147411032194048 0.0205881504343637 0.147411032194048 0.0205881504343637 0.147411032194048 0.0205881504343637];
            otherwise
                if D < 20
                    error('infbench_conbananas:OddDim','Benchmark function supports only even D.');
                else
                    error('infbench_conbananas:TooHighDim','Benchmark function supports up to D=20 dimensions.');
                end
        end
                
        % Are we transforming the entire problem to unconstrained space?
        transform_to_unconstrained_coordinates = false;
        
        % Define parameter upper/lower bounds
        lb = lb*ones(1,D);
        ub = ub*ones(1,D);
        plb = 0.1*ones(1,D);
        pub = 3*ones(1,D);
        noise = [];

        if transform_to_unconstrained_coordinates
            trinfo = warpvars(D,lb,ub,plb,pub);     % Transform to unconstrained space
            trinfo.mu = zeros(1,D);     % Necessary for retro-compatibility
            trinfo.delta = ones(1,D);
        else
            trinfo = [];
        end
        
        xmin = Mode;
        fval = ll_rosenbrock(xmin);
        xmin_post = xmin;
        fval_post = fval;

        xmin = warpvars(xmin,'d',trinfo);
        fval = fval + warpvars(xmin,'logp',trinfo);

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

        y.Post.Mean = Mean;
        y.Post.Mode = xmin_post;          % Mode of the posterior
        y.Post.ModeFval = fval_post;        
        y.Post.lnZ = lnZ;
        y.Post.Cov = Cov;

        % Compute marginals
        nx = 2^13;
        xmesh = linspace(lb(1), ub(1), nx);
                
        xpdf = zeros(1,nx);
        ypdf = zeros(1,nx);
        for i = 1:nx
            xpdf(i) = qtrapz(exp(ll_rosenbrock([xmesh(i)*ones(nx,1),xmesh(:)])));
            ypdf(i) = qtrapz(exp(ll_rosenbrock([xmesh(:),xmesh(i)*ones(nx,1)])));
        end
        xpdf = xpdf / (qtrapz(xpdf)*(xmesh(2)-xmesh(1)));
        ypdf = ypdf / (qtrapz(ypdf)*(xmesh(2)-xmesh(1)));
        
        for d = 1:D/2
            y.Post.MarginalBounds(:,[2*d-1:2*d]) = [lb(1:2); ub(1:2)];
            y.Post.MarginalPdf(2*d-1,:) = xpdf;
            y.Post.MarginalPdf(2*d,:) = ypdf;
        end

        % Save data and coordinate transformation struct
        data.trinfo = trinfo;
        y.Data = data;
        y.DeterministicFlag = true;
    end
else
    
    % Iteration call -- evaluate objective function
    
    % Transform unconstrained variables to original space
    x_orig = warpvars(x,'i',infprob.Data.trinfo);
    dy = warpvars(x,'logpdf',infprob.Data.trinfo);   % Jacobian correction
    
    % Compute log likelihood of data and possibly std of log likelihood
    LL = ll_rosenbrock(x_orig,infprob.Data);
    y_std = 0;
    y = LL + dy;
    
end

end

%--------------------------------------------------------------------------
function ll = ll_rosenbrock(theta,data)

if nargin < 2; data = []; end

sigma2 = 9;
[N,D] = size(theta);

ll = zeros(N,1);

for d = 1:(D/2)
    x = [theta(:,d*2-1), theta(:,d*2)];    
    ll = ll - ((x(:,1) .^2 - x(:,2)).^2 + (x(:,1)-1).^2/100);    
end

ll = ll - 0.5*sum(theta.^2,2)/sigma2 - 0.5*D*log(2*pi*sigma2);    

end