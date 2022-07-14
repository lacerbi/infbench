function xx = stratais(X,y,fun,Ns)
%STRATAIS Stratified importance sampling.

Nstrata = 20;
[N,D] = size(X);

Ns_perstratum = ceil(1e6/Nstrata);

% Sort by order
[~,ord] = sort(y,'descend');
X = X(ord,:);

ss = [];
ftilde = [];
for iter = 1:Nstrata
    Nmax = max(ceil(iter*N/Nstrata), D+1);
    idx = 1:Nmax;
    
    mu{iter} = mean(X(idx,:),1);
    Sigma{iter} = cov(X(idx,:));

    ss1 = mvnrnd(mu{iter},Sigma{iter},Ns_perstratum);
    
    ss = [ss; ss1];
    ftilde = [ftilde; fun(ss1)];
end

ypdf = zeros(size(ss));
for iter2 = 1:Nstrata
    ypdf = ypdf + mvnpdf(ss,mu{iter2},Sigma{iter2}) / Nstrata;
end
ww = ftilde ./ ypdf;
ww = ww / sum(ww);

y = randsample(size(ww,1),Ns,true,ww);

xx = ss(y,:);

end