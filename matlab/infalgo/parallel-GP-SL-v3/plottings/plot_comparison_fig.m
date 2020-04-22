function [] = plot_comparison_fig(th_grid, log_post_tr, theta_tr, estim_post, sim_model,...
    nr_init, nr_batch, res, other_opt)
% General function for plotting the estimated posterior and compare it to the baseline 
% (if available). 

d = th_grid.dim;
% check labels for plotting
if ~isfield(sim_model,'theta_names')
    names = cell(1,d); 
    for i = 1:d
        names{i} = ['\theta_',num2str(i)];
    end
else
    names = sim_model.theta_names;
end

if d == 1
    thx = th_grid.theta;
    
    %% fig1:
    subplot(1,2,1); % log lik
    my_shadedplot(thx, estim_post.nloglik_lb, estim_post.nloglik_ub, ...
        0.8*[1,1,1], [1,1,1]); % new loglik uncertainty
    hold on;
    plot(thx, estim_post.eloglik,'-k'); % log lik
    plot(thx, estim_post.loglik_lb,'--k'); % log lik CI (latent)
    plot(thx, estim_post.loglik_ub,'--k');
    plot(theta_tr,log_post_tr,'*b'); % data points
    xlabel('\theta');
    title('log likelihood');
    xlim([th_grid.range(1),th_grid.range(2)]);
    hold off;
    set(gcf,'Position',[50 1000 1200 450]);

    %% fig2:
    subplot(1,2,2); % estimated posterior
    
    % NOTE: WE RESCALE THESE AGAIN SO THAT THE POINT ESTIMATE OF THE DENSITY INTEGRATES TO 1
    epost = estim_post.epost; 
    [epost,c] = normalise_pdf_in_grid(epost, th_grid);
    
    my_shadedplot(thx, estim_post.post_lb/c, estim_post.post_ub/c, ...
        0.8*[1,1,1], [1,1,1]); % post uncertainty
    hold on;
    
    plot(thx, epost,'-k'); % estimated post
    if isfield(sim_model,'true_post_pdf') && ~isempty(sim_model.true_post_pdf)
    	plot(thx, sim_model.true_post_pdf,'-r'); % true post
    else
        % plot only true parameter instead
        plot(sim_model.true_theta, 0, '+r','MarkerSize',14);
    end
    
    xlabel('\theta');
    title('posterior');
    xlim([th_grid.range(1),th_grid.range(2)]);
    hold off;
    add_general_title(res);
    
elseif d == 2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % fig1: estimated posterior, fig2: uncertainty of post, fig3: true posterior
    nr_contour_lines = 25;
    thx = th_grid.theta(1,:);
    thy = th_grid.theta(2,:);
    is_truep = isfield(sim_model,'true_post_pdf2d') && ~isempty(sim_model.true_post_pdf2d);
    
    %% demo mode: plots only the estimated posterior and true posterior
    if other_opt.viz_demomode2d
        % estimated posterior
        subplot(1,2,1);
        set(gcf,'Position',[50 1000 730 320]);
        epost = estim_post.epost;
        epost = vec_to_grid_matrix(epost, th_grid);
        [epost,c] = normalise_pdf_in_grid(epost, th_grid);
        hold on;
        contour(thx, thy, epost, nr_contour_lines); % estim post
        add_datapoints_2dplot(theta_tr,nr_init,nr_batch,sim_model);
        %xlabel(names{1}); ylabel(names{2});
        hold off; box on;
        title('estimated posterior');
        
        % true posterior
        if is_truep
            subplot(1,2,2);
            true_post = vec_to_grid_matrix(sim_model.true_post_pdf2d, th_grid);
            hold on;
            contour(thx, thy, true_post, nr_contour_lines); % true post
            add_datapoints_2dplot([],[],[],sim_model);
            %xlabel(names{1}); ylabel(names{2});
            hold off; box on;
            title('exact posterior');
        end
        add_general_title(res);
        return;
    end
    
    %% fig1: loglik
    subplot(2,3,1);
    % set(gcf,'Position',[50 1000 1400 550]);
    eloglik = vec_to_grid_matrix(estim_post.eloglik, th_grid);
    hold on;
    contour(thx, thy, eloglik, nr_contour_lines); % post variance
    colorbar;
    add_datapoints_2dplot(theta_tr,nr_init,nr_batch,sim_model);
    xlabel(names{1}); ylabel(names{2});
    hold off; box on;
    title('mean of log lik');
    
    %% fig2: var of loglik
    subplot(2,3,2);
    eloglik = vec_to_grid_matrix(estim_post.varloglik, th_grid);
    hold on;
    contour(thx, thy, eloglik, nr_contour_lines); % post variance
    colorbar;
    add_datapoints_2dplot(theta_tr,nr_init,nr_batch,sim_model);
    xlabel(names{1}); ylabel(names{2});
    hold off; box on;
    title('var of log lik');
    
%     %% fig4: uncertainty of post
%     subplot(1,nr_plots,3);
%     varpost = vec_to_grid_matrix(estim_post.varpost, th_grid);
%     hold on;
%     contour(thx, thy, varpost, nr_contour_lines); % post variance
%     add_datapoints_2dplot(theta_tr,nr_init,nr_batch,sim_model);
%     xlabel(names{1}); ylabel(names{2});
%     hold off; box on;
%     title('var of estimated posterior');
    
    %% fig3: estimated posterior
    % NOTE: WE RESCALE THESE AGAIN SO THAT THE POINT ESTIMATE OF THE DENSITY INTEGRATES TO 1
    subplot(2,3,3);
    epost = estim_post.epost;
    epost = vec_to_grid_matrix(epost, th_grid);
    [epost,c] = normalise_pdf_in_grid(epost, th_grid);
    hold on;
    contour(thx, thy, epost, nr_contour_lines); % estim post
    add_datapoints_2dplot(theta_tr,nr_init,nr_batch,sim_model);
    xlabel(names{1}); ylabel(names{2});
    hold off; box on;
    title('estimated posterior');
    
    
    %% fig6&7: marginals of posterior
    % compute also marginals from the point estimate and plot them
    epost_marg = marginals_2dpdf(th_grid, epost);
    for i = 1:2
        subplot(2,3,3+i);
        hold on;
        plot(thx*(i==1)+thy*(i==2), epost_marg(i,:), '-k');
        if ~isempty(sim_model.true_theta)
            yl = ylim;
            plot(sim_model.true_theta(i)*[1,1],[0,yl(end)],'-r');
            ylim([0,yl(end)]);
        end
        hold off; box on;
        xlabel(names{i});
        xlim([th_grid.range(i,1),th_grid.range(i,2)]);
        title('estimated marginal posterior');
    end
    
    
    %% fig5: true posterior (plot only if available)
    if is_truep
        subplot(2,3,6);
        true_post = vec_to_grid_matrix(sim_model.true_post_pdf2d, th_grid);
        hold on;
    	contour(thx, thy, true_post, nr_contour_lines); % true post
        add_datapoints_2dplot([],[],[],sim_model);
        xlabel(names{1}); ylabel(names{2});
        hold off; box on;
        title('exact posterior');
    end
    add_general_title(res);
    
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Plot slices
    set(gcf,'Position',[50 1000 1600 800]);
    pls = isfield(estim_post.slice,'eloglik');
    k = 1; % height of the plot
    if pls
        k = k + 3;
        for i = 1:d
            limsi = [th_grid.range(i,1),th_grid.range(i,2)];
            % first plot row
            subplot(k,d,i);
            plot(th_grid.theta(i,:), estim_post.slice.eloglik(i,:),'-k'); % slice of estimated loglik at true val
            xlim(limsi);
            xlabel(names{i});
            if i==1
                ylabel('mean loglik (slice)');
            end
            
            % second plot row
            subplot(k,d,d+i);
            plot(th_grid.theta(i,:), sqrt(max(estim_post.slice.varloglik(i,:),0)),'-k'); % slice of stdev of loglik at true val
            xlim(limsi);
            xlabel(names{i});
            if i==1
                ylabel('stdev loglik (slice)');
            end
            
            % third plot row
            subplot(k,d,2*d+i);
            plot(th_grid.theta(i,:), estim_post.slice.epost(i,:),'-k'); % slice of estimated post at true val
            xlim(limsi);
            xlabel(names{i});
            if i==1
                ylabel('estim post (slice)');
            end
        end
    end
    
    %% Plot posterior marginals (obtained from MCMC&KDE)
    for i = 1:d
        % fourth plot row
        subplot(k,d,k*d-d+i);
        thi = th_grid.theta(i,:);
        hold on;
        plot(thi, estim_post.epost(i,:),'-k'); % point estimate of marginal post
        xlim([th_grid.range(i,1),th_grid.range(i,2)]);
        if isfield(sim_model,'true_post_pdf') && ~isempty(sim_model.true_post_pdf)
            plot(thi, sim_model.true_post_pdf(i,:),'-r'); % true post
        else
            % true post unavailable, plot true parameter instead
            yl = ylim;
            plot(sim_model.true_theta(i)*[1,1],[0,yl(end)],'-r');
            ylim([0,yl(end)]);
        end
        add_datapoint_proj_ndplot(theta_tr, nr_init, nr_batch, i); % data point projections
        xlabel(names{i});
        hold off; box on;
        if i==1
            ylabel('estim post');
        end
    end
    add_general_title(res);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting help functions:

function [] = add_general_title(res,gen_titl)
% Adds a general level title to the plot which shows TV/KL if these values are computed.

if nargin < 2
    gen_titl = 1;
end
if isfield(res,'tv') && ~isempty(res.tv)
    titl = ['TV=',num2str(mean(res.tv)),', KL=',num2str(mean(res.kl)),', L2=',num2str(mean(res.l2))];
    if gen_titl
        % suptitle(titl);
    else
        title(titl);
    end
end
end


function [] = add_datapoint_proj_ndplot(theta_tr, nr_init, nr_batch, i)
% Adds the datapoints projected to coordinate i.

if ~isempty(theta_tr)
    %plot(theta_tr(:,i),zeros(size(theta_tr(:,i))),'*k'); % data point projections
    nr_pts = size(theta_tr,1);
    plot(theta_tr(1:nr_init,i),zeros(nr_init,1),'xk');
    if nr_pts > nr_init % at least one batch
        plot(theta_tr(nr_init+1:end-nr_batch,i),zeros(size(theta_tr(nr_init+1:end-nr_batch,1))),'*k'); % old batches
        plot(theta_tr(end-nr_batch+1:end,i),zeros(size(theta_tr(end-nr_batch+1:end,1))),'*r'); % current batch
    end
end
end


function [] = add_datapoints_2dplot(theta_tr, nr_init, nr_batch, sim_model)
% Adds the training data points/old batches of selected points/current selected batch of 
% points to the figure. Also adds labels.

if ~isempty(theta_tr)
    nr_pts = size(theta_tr,1);
    plot(theta_tr(1:nr_init,1),theta_tr(1:nr_init,2),'xk'); % init points
    if nr_pts > nr_init % at least one batch
        plot(theta_tr(nr_init+1:end-nr_batch,1),theta_tr(nr_init+1:end-nr_batch,2),'*k'); % old batches
        plot(theta_tr(end-nr_batch+1:end,1),theta_tr(end-nr_batch+1:end,2),'*r'); % current batch
    end
end
if nargin == 4 && ~isempty(sim_model.true_theta)
    plot(sim_model.true_theta(1),sim_model.true_theta(2),'+k','MarkerSize',16); % true param
end
end




