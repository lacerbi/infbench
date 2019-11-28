function [] = demo_acquired_points()
% For plotting *Figures 5 and 11-13* in the paper. 
%
% Note: Instead of using the plotting functionality inside the main algorithm, we just run
% the algorithm and do the 'custom plottings' here

close all;
subplot = @(m,n,p) subtightplot (m, n, p, [0.03 0.03], [0.2 0.1], [0.1 0.1]); % better subplots

ONE_FIG = 0; % If 1 then all 9 plots to same figure, otherwise creates 2 figures
compute_only = 0;
load_from_file = 1; 
%b = 2;
b = 4;
root = ['../results/viz_demo_b',num2str(b),'/'];
gp_abc_seed = 42;

toy_model = 'banana2d';
toy_model = 'bimodal2d';
%toy_model = 'simple2d';
noise_stdev = 1;

acqs = {'MAXV','MAXV','EIV','EIV','MAXIQR','MAXIQR','IMIQR','IMIQR','IMIQR'};
nacq = length(acqs);
batch_methods = {'seq','greedy','seq','greedy','seq','greedy','seq','greedy','joint'};
m = 81;
%m = 13; % FOR QUICK TESTING
batch_sizes = [1,b,  1,b,  1,b,  1,b,b  ];
ib = (m-1)/b+1;
n_iters =     [m,ib, m,ib, m,ib, m,ib,ib];


%% main settings
nr_init = 10;
%%%nr_iter = 50;
%%%batch_size = 1;
graphics_on = -1; % plot only final result=1 / plot after each batch=2 / plot all and pause=3 / compute estimates but no plotting=-1 (for cluster)

%% gp settings
%%%gp_opt.noise_model = 1; % if 1, then use bootstrapped noise variance estimates in GP
gp_opt.meanf = 1; % 0 is zero mean GP, 1 enables const/lin/quadratic terms
gp_opt.hyp_upd_freq = 1;
gp_opt.display_type = 'off';

%% acq settings
%%%acq_opt.method = 'unif';
%%%acq_opt.batch_method = 'greedy';
acq_opt.optim_alg = 'fmincon';
acq_opt.rs.nr_init = 1000; % number of evals in random search
acq_opt.fmincon.nr_inits = 1000; % number of initial points for the multistart optimization
acq_opt.fmincon.nr_inits_grad = 10; % number of the best initial points that are actually used for multistart optimization
acq_opt.fmincon.tolf = 1e-5; % error tolerances for fmincon optimizer
acq_opt.fmincon.tolx = 1e-5;
acq_opt.direct.maxevals = 1000; % max. number of function evals  (default is 20)
acq_opt.direct.maxits = 100; % max. number of iterations  (default is 10)
acq_opt.direct.maxdeep = 100; % max. number of rect. divisions (default is 100)
acq_opt.exp.is_samples = 200; % how many samples from the importance distribution 
acq_opt.exp.nr_grid.dim1 = 100; % number of grid points for grid integration
acq_opt.exp.nr_grid.dim2 = 50;
acq_opt.display_type = 'off';

%% 2d so mcmc settings not needed here
mcmc_opt = [];

%% abc/sl lik estimator settings
lik_opt.method = 'exact'; % HERE WE HAVE EXACT EVALS
gp_opt.noise_model = 0;
lik_opt.sl.estimator = 'sl'; % 'sl', 'ubsl', 'ublogsl'
lik_opt.sl.N = [];

%% other misc settings
%%%other_opt.res_ind = []; % iterations when to compute TV/KL if cluster computation
%%%other_opt.viz_ind = 1:nr_iter; % iterations when to plot the figures
other_opt.viz_save = 0; % whether the save the plotted figures
other_opt.viz_demomode2d = 0; % for plotting 2d illustrations
other_opt.display_type = 'off';
%%%other_opt.output_folder = ['../results/',sim_model.name];
%%%other_opt.output_filename = [acq_opt.method,acq_opt.batch_method];

% get toy model settings
[th_grid,sim_model] = get_test_model(toy_model,noise_stdev);

results = cell(nacq,1);
for i = 1:nacq
    fn_i = [root,toy_model,'_simul_data_acq',num2str(i)];
    if load_from_file
        % load cached output files for faster redrawing
        load(fn_i); % get res_i
    else
        % simulate now and cache
        i
        t0=tic;
        
        % update settings with current values
        nr_iter = n_iters(i);
        batch_size = batch_sizes(i);
        
        acq_opt.method = acqs{i};
        acq_opt.batch_method = batch_methods{i};
        other_opt.res_ind = 1:n_iters(i);
        return_lvl = 2;
        
        rng(gp_abc_seed);
        res_i = run_algorithm(nr_init, nr_iter, batch_size, graphics_on, th_grid, sim_model, ...
            gp_opt, acq_opt, mcmc_opt, lik_opt, other_opt, return_lvl);
        save(fn_i,'res_i');
        
        toc(t0)
    end
    results{i} = res_i;
end

if compute_only
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% change the order of the methods for easier comparison
if 1
    ind = [1 3 5 2 4 6 7:9]; 
    results = results(ind);
    acqs = acqs(ind);
    batch_methods = batch_methods(ind);
    batch_sizes = batch_sizes(ind);
end

% IMIQR was called previously (confusingly) MIIQR/EIIQR, this changes old name to new one 
% so that also old results can be plotted using this code file
if 1
    for i = 1:nacq
        if strcmp(acqs{i},'MIIQR') || strcmp(acqs{i},'EIIQR')
            acqs{i} = 'IMIQR';
        end
    end
end

if ONE_FIG
    figure(1);
    set(gcf,'Position',[5 5 1100 1300]);
    nr_contour_lines = 25;
    nx = 3; ny = 3;
    for i = 1:nacq
        subplot(ny,nx,i);
        hold on;
        thx = th_grid.theta(1,:);
        thy = th_grid.theta(2,:);
        epost = results{i}.iter{end}.estim_post; % take the estimate of the final iteration
        epost = vec_to_grid_matrix(epost, th_grid);
        [epost,c] = normalise_pdf_in_grid(epost, th_grid);
        contour(thx, thy, epost, nr_contour_lines); % estim post
        
        % initial points
        th_tri = results{i}.theta_tr;
        plot(th_tri(1:nr_init,1),th_tri(1:nr_init,2),'xk','MarkerSize',8,'Linewidth',1.5);
        
        % acquired points
        if batch_sizes(i) > 1
            plot(th_tri((nr_init+1):(end-2*b),1),th_tri((nr_init+1):(end-2*b),2),'.k','MarkerSize',11,'Linewidth',1.5);
            % last two batches shown with different color as the other points:
            plot(th_tri((end-2*b+1):(end-b),1),th_tri((end-2*b+1):(end-b),2),'dr','MarkerSize',5,'MarkerFaceColor','r');
            plot(th_tri((end-b+1):end,1),th_tri((end-b+1):end,2),'sr','MarkerSize',5,'MarkerFaceColor','r');
        else
            plot(th_tri((nr_init+1):end,1),th_tri((nr_init+1):end,2),'.k','MarkerSize',11,'Linewidth',1.5);
        end
        hold off;
        box on;
        titl_i = [acqs{i},' (',batch_methods{i},')'];
        if 1
            % add TV to the title
            tv_i = sprintf('%.2f',results{i}.iter{end}.tv);
            titl_i = [titl_i,', TV: ',tv_i];
        end
        title(titl_i);
        if i >= 7; xlabel('\theta_1'); else; set(gca,'xtick',[]); end
        if any(i==[1 4 7]); ylabel('\theta_2'); else; set(gca,'ytick',[]); end
    end
    
    if 1
        fn = [root,toy_model,'acq_points_demo'];
        my_export_fig(fn,'-transparent','-pdf');
    end
    
else
    %% Plot two figures instead of one 3x3 figure...
    subplot = @(m,n,p) subtightplot (m, n, p, [0.04 0.03], [0.2 0.1], [0.1 0.1]); % better subplots
    nr_contour_lines = 25;    
    for F = 1:2
        figure(F);
        if F == 1
            nx = 3; ny = 2;
            set(gcf,'Position',[5 5 1100 690]);
            acqs2plot = 1:6;
        else
            nx = 3; ny = 1;
            set(gcf,'Position',[5 5 1100 320]);
            acqs2plot = 7:9;
        end

        nr = 1;
        for i = acqs2plot
            subplot(ny,nx,nr);
            hold on;
            thx = th_grid.theta(1,:);
            thy = th_grid.theta(2,:);
            epost = results{i}.iter{end}.estim_post; % take the estimate of the final iteration
            epost = vec_to_grid_matrix(epost, th_grid);
            [epost,c] = normalise_pdf_in_grid(epost, th_grid);
            contour(thx, thy, epost, nr_contour_lines); % estim post
            
            % initial points
            th_tri = results{i}.theta_tr;
            plot(th_tri(1:nr_init,1),th_tri(1:nr_init,2),'xk','MarkerSize',8,'Linewidth',1.5);
            
            % acquired points
            if batch_sizes(i) > 1
                plot(th_tri((nr_init+1):(end-2*b),1),th_tri((nr_init+1):(end-2*b),2),'.k','MarkerSize',11,'Linewidth',1.5);
                % last two batches shown with different color as the other points:
                plot(th_tri((end-2*b+1):(end-b),1),th_tri((end-2*b+1):(end-b),2),'dr','MarkerSize',5,'MarkerFaceColor','r');
                plot(th_tri((end-b+1):end,1),th_tri((end-b+1):end,2),'sr','MarkerSize',5,'MarkerFaceColor','r');
            else
                plot(th_tri((nr_init+1):end,1),th_tri((nr_init+1):end,2),'.k','MarkerSize',11,'Linewidth',1.5);
            end
            hold off;
            box on;
            titl_i = [acqs{i},' (',batch_methods{i},')'];
            if 1
                % add TV to the title
                tv_i = sprintf('%.2f',results{i}.iter{end}.tv);
                titl_i = [titl_i,', TV: ',tv_i];
            end
            title(titl_i);
            if F==1
                if i >= 4; xlabel('\theta_1'); else; set(gca,'xtick',[]); end
                if any(i==[1 4]); ylabel('\theta_2'); else; set(gca,'ytick',[]); end
            else
                xlabel('\theta_1');
                if i==7; ylabel('\theta_2'); else; set(gca,'ytick',[]); end
            end
            nr = nr + 1;
        end
        
        if 1
            fn = [root,toy_model,'acq_points_demo_',num2str(F)];
            my_export_fig(fn,'-transparent','-pdf');
        end
    end
end




