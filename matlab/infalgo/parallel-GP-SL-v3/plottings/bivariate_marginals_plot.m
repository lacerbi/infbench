function bivariate_marginals_plot(th_grid,sim_model,samples,post)
% Plots all bivariate marginal densities nicely.

FULLRANGE = 1; % if 1, plots full range of the parameters in the figures
TRUE_PARAM = 1; % if 1, plots the true parameter to the figures
MAX_S = 1000; % maximum number of samples to plot

MAX_S = min(size(samples,1),MAX_S);
inds = floor(linspace(1,size(samples,1),MAX_S));
samples1 = samples(inds,:);
d = sim_model.dim;
for i = 1:d
    for j = i:d
        if nargin < 4 && i == j
            continue;
        end
        if nargin < 4
            subplot(d-1,d-1,(i-1)*(d-1)+j-1);
        else
            subplot(d,d,(i-1)*d+j);
        end
        if j == i % diagonal -> draw marginal density
            plot(th_grid.theta(i,:),post(:,i),'-k');
            if TRUE_PARAM && isfield(sim_model,'true_theta')
                hold on;
                plot(sim_model.true_theta(i),0,'dr','MarkerSize',6,'MarkerFaceColor','r');
                hold off;
            end
            if isfield(sim_model,'theta_names')
                xlabel(sim_model.theta_names{i}); 
            end
            if FULLRANGE
                xlim([th_grid.range(i,1),th_grid.range(i,2)]);
            end
            set(gca,'ytick',[]);
        else % upper diagonal
            plot(samples1(:,j),samples1(:,i),'.k');
            if TRUE_PARAM && isfield(sim_model,'true_theta')
                hold on;
                plot(sim_model.true_theta(j),sim_model.true_theta(i),'xr',...
                    'MarkerSize',8,'MarkerFaceColor','r','Linewidth',2);
                hold off;
            end
            if isfield(sim_model,'theta_names')
                xlabel(sim_model.theta_names{j}); ylabel(sim_model.theta_names{i}); 
            end
            if FULLRANGE
                xlim([th_grid.range(j,1),th_grid.range(j,2)]);
                ylim([th_grid.range(i,1),th_grid.range(i,2)]);
            end
        end
        box on;
    end
end
end


