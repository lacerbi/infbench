function [] = demo_viz_synthetic2d_densities()
% Plots *Figure 3* for the paper. 
%
% Visualize the three synthetic 2d test problems.

close all;

dens = {'simple2d','banana2d','bimodal2d'};
n = length(dens);

for i = 1:n
    [th_grid,sim_model] = get_test_model(dens{i},2,[]);
    
    thx = th_grid.theta(1,:);
    thy = th_grid.theta(2,:);
    true_post = vec_to_grid_matrix(sim_model.true_post_pdf2d, th_grid);
    
    figure(1);
    subplot(1,n,i);
    nr_contour_lines = 30;
    contour(thx, thy, true_post, nr_contour_lines); % true post
    title(sim_model.name);
    xlabel('\theta_1');
    ylabel('\theta_2');
end

set(gcf,'Position',[50 50 400*n 300]);

% save to file
if 1
    fn = '../results/synth_densities';
    my_export_fig(fn,'-transparent','-pdf');
end
end



