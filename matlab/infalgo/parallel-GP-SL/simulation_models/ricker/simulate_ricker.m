function y = simulate_ricker(theta,N,T)
% Simulates one data set with length T from the Ricker model.
%
% INPUT:
% N - the starting population (equal to 1 in our application)
% T - the length of the data set
%
% OUTPUT:
% y - the simulated data set

r = exp(theta(1)); % the parameter we sample over is log(r)
phi = theta(2);
sigma_e = theta(3);

Ns = [N; zeros(T,1)];
srn = sigma_e*randn(T,1);
for t = 1:T
    Ns(t+1) = r*Ns(t)*exp(-Ns(t) + srn(t)); % population size
end
y = my_poissrnd(phi*Ns(2:end)); % the actual observed random variable

%y = y(51:end);
end


function r = my_poissrnd(lambda)
% Modified MATLAB poisson generator function: Removed input checkings etc. to speed up
% computations slightly. 
% TODO: This could make slightly faster by e.g. removing input checking from binornd.m
%
% Input lambda must be scalar or column vector with finite and non-negative values.
%
%   POISSRND uses a waiting time method for small values of LAMBDA,
%   and Ahrens' and Dieter's method for larger values of LAMBDA.

%   References:
%      [1]  Devroye, L. (1986) Non-Uniform Random Variate Generation,
%           Springer-Verlag.

%   Copyright 1993-2015 The MathWorks, Inc.

lambda(lambda < 0) = NaN;

%Initialize r to zero.
r = zeros(size(lambda));

r(isinf(lambda)) = Inf;

% For large lambda, use the method of Ahrens and Dieter as
% described in Knuth, Volume 2, 1998 edition.
k = find(15 <= lambda & lambda < Inf);
if ~isempty(k)
   alpha = 7/8;
   lk = lambda(k);
   m = floor(alpha * lk);

   % Generate m waiting times, all at once
   x = randg(m);
   t = (x <= lk);

   % If we did not overshoot, then the number of additional times
   % has a Poisson distribution with a smaller mean.
   r(k(t)) = m(t) + my_poissrnd(lk(t)-x(t));

   % If we did overshoot, then the times up to m-1 are uniformly
   % distributed on the interval to x, so the count of times less
   % than lambda has a binomial distribution.
   if ~all(t)
       r(k(~t)) = binornd(m(~t)-1, lk(~t)./x(~t));
   end
end

% For small lambda, generate and count waiting times.
j = find(lambda < 15);
p = zeros(numel(j),1,'like',lambda);
while ~isempty(j)
    p = p - log(rand(numel(j),1,'like',lambda));
    t = (p < lambda(j));
    j = j(t);
    p = p(t);
    r(j) = r(j) + 1;
end

% Return NaN if LAMBDA is negative.
r(isnan(lambda)) = NaN;
end




