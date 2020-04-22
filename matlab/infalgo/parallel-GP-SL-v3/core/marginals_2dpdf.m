function margs = marginals_2dpdf(th_grid, dens)
% Compute the marginal densities from a 2d density evaluated at grid points and normalise
% them so that they integrate to 1.

margs = NaN(2,size(dens,1)); % assuming square grid
th_gridi.dim = 1;
for i = 1:2
    margs(i,:) = sum(dens,i);
    th_gridi.theta = th_grid.theta(i,:);
    margs(i,:) = normalise_pdf_in_grid(margs(i,:), th_gridi); % normalise
end
end



