function [th_opt, f_opt] = optimize_acq(f, th_grid, th_tr, acq_opt, rep_th)
% Wrapper of some algorithms for optimising a given acquisition function f. 
% NOTE: 'rep_th' tells if the optimisation is in batch case so that instead of optimising
%   over Theta, we optimise jointly over Theta^rep_th.
% NOTE: function handle f must be to a function to be minimized!
% TODO: improved fmincon with gradient computations. 

display_type = acq_opt.display_type;
if ~strcmp(display_type,'off')
    disp('Optimizing the acquisition function...');
end
th_range = th_grid.range;
lb = repmat(th_range(:,1)',1,rep_th);
ub = repmat(th_range(:,2)',1,rep_th);
d = th_grid.dim;

if strcmp(acq_opt.optim_alg,'grid') 
    %% Naive minimization in a grid, useful for dim <= 2 (or maybe 3)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if rep_th > 1
        error('Grid search not implemented for the "joint batch" case.');
    elseif d > 2
        error('Grid search is for dim <= 2.');
    end
    if ~strcmp(display_type,'off')
        disp(['Using grid search with ', num2str(size(th_grid,1)), ' points.']);
    end
    if d==1
        f_all = f(th_grid.theta');
        [f_opt,opt_ind] = min(f_all);
        th_opt = th_grid.theta(:,opt_ind);
    else
        f_all = f(th_grid.theta2d');
        [f_opt,opt_ind] = min(f_all);
        th_opt = th_grid.theta2d(:,opt_ind);
    end
    
%     % Idea: improve the result by grid search by optimization
%     if 0
%         alg = 'interior-point'; %'sqp'
%         fmin_opts = optimset('Algorithm',alg, 'TolFun',...
%             acq_opt.fmincon.tolf,'TolX',acq_opt.fmincon.tolx, 'Display',display_type);
%         th0 = x_opt;
%         try
%             [th,f_val] = fmincon(f,th0,[],[],[],[],lb,ub, [], fmin_opts);
%         catch
%             warning('Optimization with current point failed.');
%             %x0, f(x0)
%             f_val = inf;
%         end
%         % check if the optimized point is better than the best point by random search
%         if f_val < f_opt
%             f_opt = f_val;
%             th_opt = th;
%         end
%     end
    
elseif strcmp(acq_opt.optim_alg,'rs') 
    %% Random search
    %%%%%%%%%%%%%%%%
    
    inits = inits_for_opt(acq_opt.rs.nr_init, th_grid, rep_th); % select points
    f_all = f(inits); % evaluate
    [f_opt,opt_ind] = min(f_all);
    th_opt = inits(opt_ind,:);
    
elseif strcmp(acq_opt.optim_alg,'fmincon') 
    %% Multistart gradient based optimization
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    inits = inits_for_opt(acq_opt.fmincon.nr_inits, th_grid, rep_th);
    
    % interior point algorithm is the default
    alg = 'interior-point'; %'sqp'
    fmin_opts = optimset('Algorithm',alg, 'TolFun',...
        acq_opt.fmincon.tolf,'TolX',acq_opt.fmincon.tolx, 'Display',display_type);
    f_opt = inf;
    th_opt = NaN;
    
    % Choose (possibly) only some of the provided initial points to be used 
    % as initial points for the optimization procedure 
    % -> if the initial point has very low value, it might lead to local 
    % optimum only (or fail to optimize) and thus it might be useful to 
    % not waste computation time on such optimization attempts.
    %
    % Additional heuristics idea: use also current data points as initial points for 
    % multistart optimization (this might help to avoid the (unlikely) 
    % situation where all the initial points are in areas with very small 
    % acq. values and/or flat region and thus might lead to local optima)
    %inits = [inits; th_tr];
    
    f_all = f(inits);
    [~,opt_ind] = sort(f_all);
    opt_ind = opt_ind(1:min(length(opt_ind),acq_opt.fmincon.nr_inits_grad));
    inits = inits(opt_ind,:); 
    f_init = f_all(opt_ind);
    
    % optimise with each chosen initial point
    for rep = 1:size(inits,1)
        th0 = inits(rep,:);
        try
            [th,f_val] = fmincon(f,th0,[],[],[],[],lb,ub, [], fmin_opts);
        catch
            warning('Optimization with an initial point failed.');
            %th0, f_init(rep)
            %keyboard;
            continue;
        end
        
        % check if the current point is better than any point found this far
        if f_val < f_opt
            f_opt = f_val;
            th_opt = th;
            ind_opt = rep;
        end
    end
    if ~strcmp(display_type,'off')
        disp(['Best acq. value found with index ', num2str(ind_opt)]);
    end

elseif strcmpi(acq_opt.optim_alg,'direct') 
    %% DIRECT algorithm
    %%%%%%%%%%%%%%%%%%%
    
    if ~strcmp(display_type,'off')
        disp('Using DIRECT to optimize.');
    end
    problem.f = @(th) f(th');
    opts.maxevals = acq_opt.direct.maxevals; % max. number of function evals (default is 20)
    opts.maxits = acq_opt.direct.maxits; % max. number of iterations (default is 10)
	opts.maxdeep = acq_opt.direct.maxdeep; % max. number of rect. divisions (default is 100)
    
    % evaluate!
    [ret_minval,final_xatmin,history] = direct(problem, [lb(:),ub(:)], opts); 
    %nr_eval = history(end,2)
    
    th_opt = final_xatmin';
    f_opt = ret_minval;
    if ~strcmp(display_type,'off')
        final_res = history(end,:)
    end
    
elseif strcmpi(acq_opt.optim_alg,'cmaes') 
    %% CMAES
    %%%%%%%%
    
    if ~strcmp(display_type,'off')
        disp('Using CMAES to optimize.');
    end
    
    % using similar options as in AGP method, Neurocomputing 2018
    opts.LBounds = lb(:);
    opts.UBounds = ub(:);
    opts.DispModulo  = Inf;
    opts.DispFinal = 'off';
    opts.DispFinal = 'on';
    if strcmp(display_type,'off')
        opts.DispFinal = 'off';
    end
    opts.SaveVariables = 'off';
    opts.LogModulo = 0;
    opts.CMA.active = 1;
    
    init = inits_for_opt(10, th_grid, rep_th); % select 10 points from unif prior
    %init = init(:);
    stdevs = [];
    
    % optimize!
    f1 = @(th) f(th');
    [xmin,fmin,counteval,stopflag] = cmaes_new(f1, init', stdevs, opts);
    th_opt = xmin;
    f_opt = fmin;
    if ~strcmp(display_type,'off')
        init, counteval,stopflag
    end
    
else
    %% Other optimization algorithms could be added here
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    error('Incorrect optimization algorithm or the algorithm is not implemented.');
end
if ~strcmp(display_type,'off')
    disp('Optimization of the acq function completed.');
end
%% check if some reasonable solution was found
if ~isfinite(th_opt)
    error('Optimization of the acquisition function failed.');
end
end


function inits = inits_for_opt(nr_pts, th_grid, rep_th)
% Select initial points for optimisation from uniform distribution. Handles the case where
% a batch of points is selected and optimised jointly. 

d = th_grid.dim;
dopt = d*rep_th;
inits = zeros(nr_pts,dopt);
for i = 1:rep_th
    inits(:,(i-1)*d+1:(i*d)) = acquire_init_batch(nr_pts, th_grid); % nr_pts x th.dim
end
end





