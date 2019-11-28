function [kl_div,tv,l2] = compute_dist(th_grid, approx_dist, true_dist)
% Computes the KL divergence, TV and L^2 distance between two (possibly unnormalized) and 
% discretized densities in 1d or 2d case. If dim > 2, then these values between all the 
% marginal densities are computed and returned as column vector.  
%
% In 1d case the densities must be vectors with equal length. In 2d case the input 
% densities must be either vectors or matrices as long as the elements are in 
% correspondence. If dim > 2, then each column of the input densities must contain the 
% marginal densities (e.g. obtained by KDE) and must be of the same size. 
%
% If the true_dist = p and approximate_dist = q, then approximately
% KL = sum_i p_i*log(p_i/q_i)
% TV = 1/2*sum_i |p_i-q_i|
% L2 = sum_i (p_i-q_i)^2/dx (Note the extra scaling with grid element volume dx)
% where p_i, q_i are bin probabilities of each discretisation grid points.

% tolerance for making sure that small bins do not cause numerical issues in KL case
tol = 1e-6;

% if 2d, computation goes similarly after concatenating the bins into vectors
d = th_grid.dim;
if d <= 2
    true_dist = true_dist(:);
    approx_dist = approx_dist(:);
    if length(true_dist) ~= length(approx_dist)
        warning('Sizes of the pdfs do not match. Returning NaN.');
        kl_div = NaN; tv = NaN; l2 = NaN;
        return;
    end
    
    % Make sure distributions are normalized and sum of bins == 1:
    true_dist = true_dist ./ sum(true_dist);
    approx_dist = approx_dist ./ sum(approx_dist);

    % KL: deal with small probability values in some bins i.e. if true_dist_i == 0, then 
    % the corresponding contribution to the sum is also 0 despite value of approximative distribution
    bins_to_include = find(true_dist > tol);
    if isempty(bins_to_include)
        kl_div = NaN;
    else
        true_dist4kl = true_dist(bins_to_include);
        approximate_dist4kl = approx_dist(bins_to_include);
        kl_div = sum(true_dist4kl .* log(true_dist4kl ./ approximate_dist4kl));
    end
    
    % TV and L2:
    tv = 0.5 * sum(abs(true_dist - approx_dist));
    l2 = sum((true_dist - approx_dist).^2) / th_grid.dx;
    
elseif d > 2 
    % We compute KL/TV/L2 between all marginal densities only
    if size(true_dist,1) ~= size(approx_dist,1) || size(true_dist,2) ~= size(approx_dist,2)
        warning('Sizes of the pdfs do not match. Returning NaN.');
        kl_div = NaN(d,1); tv = NaN(d,1); l2 = NaN(d,1);
        return;
    end

    % computing KL/TV/L2 for each dimension separately
    kl_div = zeros(d,1); tv = zeros(d,1); l2 = zeros(d,1);
    th_gridi.dim = 1;
    for i = 1:d
        th_gridi.theta = th_grid.theta(i,:);
        th_gridi.dx = (th_grid.range(i,2)-th_grid.range(i,1))/length(th_gridi.theta);
        [kl_div(i),tv(i),l2(i)] = compute_dist(th_gridi,true_dist(i,:),approx_dist(i,:));
    end
end
end




