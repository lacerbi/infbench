function y = simulate_ricker(theta,N,T)
% Simulates one data set with length T from the Ricker model.
%
% INPUT:
% theta - parameters
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
y = poissrnd_fast(phi*Ns(2:end)); % the actual observed random variable

end




