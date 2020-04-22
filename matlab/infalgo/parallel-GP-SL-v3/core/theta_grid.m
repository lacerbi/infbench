function th_grid = theta_grid(th_range, nr_grid)
% Constructs a grid structure required for evaluations/plotting of the simulation model. 
% Each row of th_range gives the range of the corresponding dimension of the parameter 
% theta.

% th_range must be a dim x 2 matrix of ranges of the resulting grid
if size(th_range,2) ~= 2 || isempty(th_range)
    error('Dimensions of range of the grid are incorrect.');
end

th_grid.dim = size(th_range,1);
th_grid.range = th_range;

if th_grid.dim == 1
    th_grid.theta = linspace(th_range(1),th_range(2),nr_grid); 
    th_grid.dx = (th_range(2)-th_range(1))/nr_grid;
    
elseif th_grid.dim == 2
    
    % construct a 2d grid using meshgrid and save also "marginal" grids
    n1 = nr_grid;
    n2 = n1;
    th1 = linspace(th_range(1,1),th_range(1,2),n1);
    th2 = linspace(th_range(2,1),th_range(2,2),n2);
    th_grid.theta = [th1(:)';th2(:)'];
    [th1,th2] = meshgrid(th1,th2);
    th_grid.theta2d = [th1(:), th2(:)]';
    th_grid.dx = prod(th_grid.range(:,2)-th_grid.range(:,1))/nr_grid^2;
    
else % dim > 2
    
    % construct equidistant grid for marginals densities only (the whole 
    % n-dim grid would obviously require huge amount of points)
    th_grid.theta = zeros(th_grid.dim,nr_grid);
    th_grid.dx = zeros(th_grid.dim,1);
    for i = 1:th_grid.dim
        th_grid.theta(i,:) = linspace(th_range(i,1),th_range(i,2),nr_grid); 
        th_grid.dx(i) = (th_range(i,2)-th_range(i,1))/nr_grid;
    end
end
end




