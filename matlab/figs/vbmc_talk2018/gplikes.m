%


mu = 5;
sigma = 0.1;
xx = linspace(-10,10,1e5);
yy = normpdf(xx,mu,sigma);
ly = normlogpdf(xx,mu,sigma);

subplot(1,3,1);
plot(xx,yy);

subplot(1,3,2);
plot(xx,sqrt(yy));

subplot(1,3,3);
plot(xx,ly);
