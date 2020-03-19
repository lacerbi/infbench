%KRAJBICH2010
% cols: subject,trial,fixnum,fixdur,leftval,rightval,rt,choice,roi,revfixnum
function data = get_krajbich2010_data(subjid,test_flag)

if nargin < 1 || isempty(subjid); subjid = 13; end
if nargin < 2 || isempty(test_flag); test_flag = false; end

mat = csvread('krajbichetal2010_data.csv',1);
mat = mat(mat(:,1) == subjid,:);

dt = 0.1;   % 100 ms time bins
data.dt = dt;

if test_flag
    trials_idx = mat(1,2);    
else
    trials_idx = unique(mat(:,2))';
end
Ntrials = numel(trials_idx);

data.leftval = zeros(Ntrials,1);
data.rightval = data.leftval;
data.choice = data.leftval;
data.totfixdur = data.leftval;
data.totfixdurbin = data.leftval;

ii = 1;
for iTrial = trials_idx
    mat_t = mat(mat(:,2) == iTrial,:);
    
    data.leftval(ii) = mat_t(1,5);
    data.rightval(ii) = mat_t(1,6);
    data.choice(ii) = mat_t(1,8);
    data.totfixdur(ii) = sum(mat_t(:,4))/1000;
    data.totfixdurbin(ii) = round(data.totfixdur(ii)/dt);
    
    ii = ii + 1;
end

data.choice(data.choice == 0) = 2;

Nbins = ceil(max(data.totfixdur)/dt);
data.fixmat = ones(Ntrials,Nbins);

ii = 1;
for iTrial = trials_idx
    mat_t = mat(mat(:,2) == iTrial,:);
    fixtime = cumsum(mat_t(:,4))/1000;
    
    % Assign time bins during trial to fixated item
    t = 0.5*dt;
    for iBin = 1:Nbins
        idx = find(t < fixtime,1);
        if isempty(idx); break; end        
        data.fixmat(ii,iBin) = mat_t(idx,9);       
        t = t + dt;
    end
    ii = ii + 1;
end



end