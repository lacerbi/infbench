function y = normpdf_fast(x,mu,sigma)
% This code is much faster than normpdf if called repeatedly
% Return NaN for out of range parameters.

sigma(sigma <= 0) = NaN;
y = exp(-0.5 * ((x - mu)./sigma).^2) ./ (sqrt(2*pi) .* sigma);

end

