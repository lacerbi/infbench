% krajbich2010_test

data = get_krajbich2010_data(13,1);
params = [1,1,0.3,0.05];

Ns = 1e6;
resp = krajbich2010_gendata(params,ones(Ns,1),data);
Nbins = data.totfixdurbin;

edges = (1:Nbins+1)-0.5;

Nleft = histcounts(resp(resp(:,1) == 1,2),edges)/Ns;
Nright = histcounts(resp(resp(:,1) == 2,2),edges)/Ns;

plot(1:Nbins,Nleft,'b:','LineWidth',2); hold on;
plot(1:Nbins,Nright,'g:','LineWidth',2);

[ll,Pr_mat] = krajbich2010_loglike(params,data);

plot(1:Nbins,Pr_mat(1,1:Nbins,1),'b-','LineWidth',2); hold on;
plot(1:Nbins,Pr_mat(1,1:Nbins,2),'g-','LineWidth',2);

[ll_ibs,var_ll] = ibslike(@krajbich2010_gendata,params,[data.choice,data.totfixdurbin],[],[],data);
ll