function resp = krajbich2010_gendata(params,idx,data)
%KRAJBICH2010_GENDATA Generate responses for Krajbich et al. model

% Assign parameters from parameter vector
sigma = params(1)*sqrt(data.dt);
d = params(2)*data.dt;
theta = params(3);
lambda = params(4);

Nt = numel(idx);
Nbins = size(data.fixmat,2);

% Random noise matrix
nn_mat = sigma*randn(Nt,Nbins);

leftval = data.leftval(idx);
rightval = data.rightval(idx);
fixmat = data.fixmat(idx,:);

% Drift matrix
drift_mat = d*(bsxfun(@times,leftval - theta*rightval,fixmat == 1) ...
    + bsxfun(@times,theta*leftval - rightval,fixmat == 2));

% Generate traces
trace = [cumsum(nn_mat + drift_mat,2),ones(Nt,1)];

% Find boundary passage
[~,crossed_idx] = max(abs(trace) >= 1,[],2);

% Get choice at passage
ind = sub2ind(size(trace),(1:Nt)',crossed_idx); 
choice = trace(ind);
choice(choice >= 0) = 1;
choice(choice < 0) = 2;

resp = [choice,crossed_idx];

% Add lapse trials
lapse_idx = rand(Nt,1) < lambda;
Nlapse = sum(lapse_idx);
resp_lapse = [randi(2,[Nlapse,1]),randi(Nbins,[Nlapse,1])];

resp(lapse_idx,:) = resp_lapse;


end