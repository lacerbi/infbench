function [] = plot_acq_batch(th_grid, acq_grid, th_tr, y_tr, acq_opt)
% Plots the results of the acq optimization in the (non-greedy) batch case. 
% In 1d case, plots a 2d figure of the acq surface (should be symmetric).

% check the input struct
opt_f = acq_grid.opt_f;
opt_th = acq_grid.opt_th;
acq = []; 
if isfield(acq_grid,'acq')
    acq = acq_grid.acq;
end

titl = [acq_opt.method, ', t = ', num2str(length(y_tr)+1)];

thx = th_grid.theta(1,:);
thy = thx;
nr_cont_lines = 50;
ra = [th_grid.range(1,1), th_grid.range(1,2)]; % limits

hold on;
contour(thx, thy, acq, nr_cont_lines);
plot(opt_th(1),opt_th(2),'b*'); % plot optimum(s) i.e. the latest eval point(s)
if 1
    % Plot also datapoints
    plot(ra(1)*ones(length(th_tr)),th_tr(:)','k*');
    plot(th_tr(:)',ra(1)*ones(length(th_tr)),'k*');
end
plot(ra,ra,'k--') % diagonal plot
hold off;
colorbar;
xlim(ra); ylim(ra);
box on;
set(gcf,'Position',[50 50 480 320]);
title(titl,'Interpreter','none');

end



