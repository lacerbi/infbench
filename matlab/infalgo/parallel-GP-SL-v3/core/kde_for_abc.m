function dens_estim = kde_for_abc(grid_th,samples,use_ksdensity)
% Wrapper for KDE. Used when (ABC) samples need to be 'converted' to a density function 
% evaluated in a grid. This code works in 1d or 2d case. If dim >= 3, then marginal 
% densities are similarly estimated using 1d KDE.

if nargin < 3
    use_ksdensity = 1;
end

d = grid_th.dim;
if d == 1
    if length(samples) > 1
        if use_ksdensity
            % use kde from matlab toolbox
            dens_estim = kde_matlab(samples,grid_th.theta);
        else
            % use kde code from the kde toolbox
            % lcv = leave-one-out likelihood criterion 
            p = kde(samples', 'lcv'); 
            dens_estim = evaluate(p,grid_th.theta);
        end
    else
        dens_estim = ones(size(grid_th.theta));
    end
    dens_estim = dens_estim(:);
    
elseif d == 2
    % take kde for the joint 2d density
    
    if use_ksdensity
        % use kde from matlab toolbox
        dens_estim = kde_matlab(samples,grid_th.theta2d');
    else
        % use 2d kde code from the kde toolbox
        % lcv = leave-one-out likelihood criterion
        p = kde(samples', 'lcv');
        dens_estim = evaluate(p,grid_th.theta');
    end
    dens_estim = dens_estim(:);
    
else 
    % Compute kde for marginal densities only
    dens_estim = zeros(size(grid_th.theta,2),d);
    grid_i.dim = 1;
    for i = 1:d
        grid_i.theta = grid_th.theta(i,:);
        dens_estim(:,i) = kde_for_abc(grid_i,samples(:,i),use_ksdensity);
    end
end
end


function dens_estim = kde_matlab(samples,x)
% Wrapper for matlab kde.

% quick fix: needed because mcmc toolbox also contains function named 'mad'
%curf = pwd; 
%cd /opt/matlab2016b/toolbox/stats/stats/;
dens_estim = ksdensity(samples,x);
%cd(curf);
end





