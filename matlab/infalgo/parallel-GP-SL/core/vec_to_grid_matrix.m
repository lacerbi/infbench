function grid_matrix = vec_to_grid_matrix(v, th_grid)
% Converts a vector to matrix, basically the inverse operaration for matlab's (:) in the 
% case of where our grid structure is used.
% This is needed when the likelihood values computed and stored in a vector need to be
% plotted as a contour.

if(th_grid.dim ~= 2)
    %error('This should be used only when dim==2.');
    grid_matrix = v;
    return;
end

v = v(:);
n = size(th_grid.theta,2);
m = length(v)/n;
if mod(m,1) ~= 0 || n*m ~= length(v)
    error('Size mismatch.');
end
grid_matrix = reshape(v,n,m);

end




