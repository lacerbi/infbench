function [norm_pdf,c] = normalise_pdf_in_grid(unn_pdf, th_grid)
% Numerically normalise unnormalised pdf computed in a grid so that it integrates to 1. 
% If 2d case, then the output is converted to matrix corresponding to the 2d grid.

d = th_grid.dim;
if d == 1
    c = trapz(th_grid.theta, unn_pdf);
    norm_pdf = unn_pdf/c;
    
elseif d == 2
    th1 = th_grid.theta(1,:);
    th2 = th_grid.theta(2,:);
    unn_pdf = vec_to_grid_matrix(unn_pdf(:), th_grid); % make it a matrix
    c = trapz(th2,trapz(th1,unn_pdf,2));
    norm_pdf = unn_pdf/c; % normalise to sum to one in the 2d grid
    
else
    % normalise each marginal given in 'unn_pdf'
    norm_pdf = NaN(d,size(th_grid.theta,2));
    grid_thi.dim = 1;
    c = NaN(1,d);
    for i = 1:th_grid.dim
        grid_thi.theta = th_grid.theta(i,:);
        [norm_pdf(i,:),c(i)] = normalise_pdf_in_grid(unn_pdf(i,:), grid_thi);
    end
end
end






