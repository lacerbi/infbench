function p = normcdf_fast(z)
% This code is much faster than normcdf if called repeatedly

p = 0.5 * erfc(-z ./ sqrt(2));

end

