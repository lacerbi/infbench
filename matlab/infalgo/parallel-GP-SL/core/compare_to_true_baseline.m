function res = compare_to_true_baseline(th_grid, estim_post, sim_model)
% Computes KL, TV and L^2 between the estimated and the exact posterior pdf. Computations
% are done using the grid (if dim <= 2) or using the MCMC samples (if dim > 2).  

%% Compute KL/TV/L2 between estimated and true posterior:
kl = []; tv = []; l2 = [];
d = th_grid.dim;
res.dim = d;
estim_pdf = estim_post.epost;
if d == 2 && isfield(sim_model,'true_post_pdf2d') && ~isempty(sim_model.true_post_pdf2d)
    % 2d case
    true_pdf = sim_model.true_post_pdf2d;
    [kl,tv,l2] = compute_dist(th_grid, true_pdf, estim_pdf);

elseif isfield(sim_model,'true_post_pdf') && ~isempty(sim_model.true_post_pdf)
    % 1d or >=3d case
    true_pdf = sim_model.true_post_pdf;
    [kl,tv,l2] = compute_dist(th_grid, true_pdf, estim_pdf);
    
end
res.kl = kl; res.tv = tv; res.l2 = l2;

% Note: sometimes the density is estimated to be very peaked, in this case the density may
% evaluate to be zero in each grid point and, consequently, the TV and L2 become NaN.
%if any(~isfinite(res.tv)) || any(~isfinite(res.l2))
%    % note: KL can be Inf even when computations itself succeed
%    keyboard;
%end

%% Compare point estimate to the true parameter:
r = [];
if isfield(sim_model,'true_theta')
    r = rel_error(th_grid, estim_pdf, sim_model);
end
res.rel_err = r;
end

%##########################################################################

function r = rel_error(th_grid, estim_post, sim_model)
% Compute the relative error between the point estimate of the posterior (mean) and the 
% true parameter value that was used to generate the data set. This is computed for each 
% component of the parameter vector separately. 

d = th_grid.dim;
if d == 1 
    estim_post = normalise_pdf_in_grid(estim_post, th_grid);
    estim_post = estim_post(:)';
elseif d == 2
    % 2d case we need to compute the marginals from the 2d joint posterior
    estim_post = vec_to_grid_matrix(estim_post, th_grid);
    estim_post = marginals_2dpdf(th_grid, estim_post);
else
    estim_post = normalise_pdf_in_grid(estim_post, th_grid);
end
r = NaN(1,d); e = NaN(1,d);
for i = 1:d
    e(i) = trapz(th_grid.theta(i,:), th_grid.theta(i,:).*estim_post(i,:));
    r(i) = abs((e(i) - sim_model.true_theta(i))/sim_model.true_theta(i));
end
end




