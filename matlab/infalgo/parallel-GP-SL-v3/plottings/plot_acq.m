function [] = plot_acq(th_grid, acq_grid, th_tr, y_tr, acq_opt)
% Plots the results of the acq optimization in sequential or greedy batch case. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Additional settings
plot_acq_exact = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d = th_grid.dim;
% check the input struct
opt_f = acq_grid.opt_f;
opt_th = acq_grid.opt_th;
acq = []; acq_exact = []; 
if isfield(acq_grid,'acq')
    acq = acq_grid.acq;
end
if isfield(acq_grid,'acq_exact')
    acq_exact = acq_grid.acq_exact;
end

plot_acq_exact = (plot_acq_exact && ~isempty(acq_exact));
nr_subplots = 1 + plot_acq_exact;
titl = [acq_opt.method, ', t = ', num2str(length(y_tr)+1)];

%""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""  
%% plot acq. function values in the grid for closer inspection
if d == 1
    max_a = max(acq(isfinite(acq)));
    min_a = min(acq(isfinite(acq)));
    if isempty(max_a) || isempty(min_a)
        max_a = 1; min_a = 0;
    end
    if max_a <= min_a
        max_a = min_a+1; % for avoiding rare issue
    end
    
    if plot_acq_exact
        max_var = max(acq_exact);
        min_var = min(acq_exact);
        if max_var <= min_var
            max_var = min_var+1; % for avoiding rare issue
        end
        
        subplot(1,nr_subplots,2);
        hold on;
        plot(th_grid.theta,acq_exact,'-k');
        if isfield(acq_grid,'cur_acq')
            plot([th_grid.range(1), th_grid.range(2)],acq_grid.cur_acq*[1,1],'--r'); % current uncertainty value
            max_var = max(max_var,acq_grid.cur_acq);
            min_var = min(min_var,acq_grid.cur_acq);
        end
        for i = 1:length(opt_th)
            plot(opt_th(i)*[1,1],[min_var,max_var],'b-'); % plot optimum (i.e. the latest eval point)
        end
        hold off;
        box on;
        xlim([th_grid.range(1), th_grid.range(2)]);
        ylim([min_var,max_var]);
        ylabel('Acquisition value');
        title('Exact acq (simulation)');
    end
    
    subplot(1,nr_subplots,1);
    hold on;
    plot(th_grid.theta,acq,'-k'); % plot acq-value
    if isfield(acq_grid,'cur_acq')
        plot([th_grid.range(1), th_grid.range(2)],acq_grid.cur_acq*[1,1],'--r'); % current uncertainty value
        max_a = max(max_a,acq_grid.cur_acq);
        min_a = min(min_a,acq_grid.cur_acq);
    end
    for i = 1:length(opt_th)
        plot(opt_th(i)*[1,1],[min_a,max_a],'b-'); % plot optimum (i.e. the latest eval point)
    end
    hold off;
    box on;
    % data points not plotted here...
    xlim([th_grid.range(1), th_grid.range(2)]);
    ylim([min_a,max_a]);
    set(gcf,'Position',[50 50 nr_subplots*440 300]);
    ylabel('Acquisition value');
    title(titl,'Interpreter','none');

%""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""      
elseif d == 2
    
    if plot_acq_exact
        acq_exact = vec_to_grid_matrix(acq_exact, th_grid);
    end
    
    thx = th_grid.theta(1,:);
    thy = th_grid.theta(2,:);
    nr_cont_lines = 50;
    if plot_acq_exact
        subplot(1,nr_subplots,2);
        hold on;
        contour(thx, thy, acq_exact, nr_cont_lines); % var of the post curve
        plot(opt_th(:,1),opt_th(:,2),'b*'); % plot optimum (i.e. the latest eval point)
        plot(th_tr(:,1),th_tr(:,2),'k.','MarkerSize',10); % plot also datapoints
        colorbar;
        hold off;
        box on;
        title('Exact acq (simulation)');
    end
    
    subplot(1,nr_subplots,1);
    hold on;
    contour(thx, thy, acq, nr_cont_lines);
    plot(opt_th(:,1),opt_th(:,2),'b*'); % plot optimum(s) i.e. the latest eval point(s)
    plot(th_tr(:,1),th_tr(:,2),'k.','MarkerSize',10); % plot also datapoints
    hold off;
    colorbar;
    xlim([th_grid.range(1,1), th_grid.range(1,2)]);
    ylim([th_grid.range(2,1), th_grid.range(2,2)]);
    box on;
    set(gcf,'Position',[50 50 nr_subplots*440 290]);
    title(titl,'Interpreter','none');
    
%""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""    
else % dim > 2

    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
    batch_size = size(opt_th,1);
    for i = 1:d
        for j = 1:batch_size
            loc = [0.05 + 0.2*(j-1), (d-i+1)/d/1.15];
            str = num2str(opt_th(j,i));
            text(loc(1),loc(2),str);
        end
    end
    box on;
    title('Acquired location(s)');
    set(gcf,'Position',[50, 50, 100 + 100*batch_size, 250]);
    
%     disp(' ');
%     if min(size(opt_th)) == 1
%         disp(['The latest acquired point = ', num2str(opt_th(:)')]);
%     else
%         disp('The latest batch of acquired points:');
%         opt_th
%     end
%     disp(' ');
end
end




