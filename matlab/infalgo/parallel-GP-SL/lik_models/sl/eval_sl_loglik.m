function log_val = eval_sl_loglik(obs_summary, simul_summaries, estimator)
% Evaluates log of SL
% d == dim of summary vector, first dim of simul_summaries
% N == dim of repeated samples at proposed theta, second dim of simul_summaries

[d,N] = size(simul_summaries);

if strcmp(estimator,'sl')
    % log of (standard) SL
    if d > N
        error('Dim. of summary vector > N in SL.');
    end
    
    use_robust_vcov = 0;
    
    if use_robust_vcov
        [Sigma,invSigma,muhat,logdetSd2] = robust_vcov(simul_summaries);
        dx = obs_summary(:) - muhat(:);
        log_val = -logdetSd2 - 0.5*(dx'*invSigma*dx + d*log(2*pi));
    else
        % naive implementation
        ml_mean = mean(simul_summaries');
        ml_cov = cov(simul_summaries');
        %ml_cov = robustcov(simul_summaries');
        dx = obs_summary(:) - ml_mean(:);
        
        jitter = 1e-12;
        ml_cov = ml_cov + jitter*eye(size(ml_cov));
        
        log_val = -0.5*(logdet(ml_cov) + dx'*(ml_cov\dx) + d*log(2*pi));
    end
    %if abs(log_val) > 10^6
    %    keyboard;
    %end
    
elseif strcmp(estimator,'ubsl')
    % log of uBSL as it is called in the Bayesian synthetic likelihood paper 2018
    % This uses the unbiased estimator of the normal density
    if d+3 > N
        error('Dim. of summary vector > N-3 in uBSL.');
    end

    ml_mean = mean(simul_summaries');
    ml_cov = cov(simul_summaries');
    dx = obs_summary(:) - ml_mean(:);
    
    jitter = 1e-12;
    ml_cov = ml_cov + jitter*eye(size(ml_cov));
    
    Mn = (N-1)*ml_cov;
    arg_phi = Mn - dx*dx'/(1-1/N);

    [~,p] = chol(arg_phi);

    if (p == 0) % positive definite
        log_val = -d/2*log(2*pi) + ckv(d,N-2) - ckv(d,N-1) - d/2*log(1-1/N) ...
            - (N-d-2)/2*logdet(Mn) + (N-d-3)/2*logdet(arg_phi);
    else % not positive definite => value 0, log value -inf
        warning('-inf value in UBSL');
        log_val = -Inf;
    end
else
    % As in Variational Bayes with synthetic likelihood paper 2018
    % This is the unbiased estimator for the log of normal density
    if d+2 > N
        error('Dim. of summary vector > N-2 in ulogSL.');
    end
    
    ml_mean = mean(simul_summaries');
    ml_cov = cov(simul_summaries');
    dx = obs_summary(:) - ml_mean(:);
    
    log_val = -d/2*log(2*pi) - 0.5*(logdet(ml_cov) + d*log((N-1)/2) ...
        - sum(psi((N-(1:d))/2)) + (N-d-2)/(N-1)*(dx'*(ml_cov\dx)) - d/N);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ld = logdet(A)
% Compute log of the determinant using Cholesky
U = chol(A);
ld = 2*sum(log(diag(U)));
end


function f = ckv(k,nu)
% Computes log of c(k,nu), see Ghurye & Olkin (1969)
f = -k*nu/2*log(2) - k*(k-1)/4*log(pi) - sum(gammaln(0.5*(nu-(1:k)+1)));
end


function [Sigma,invSigma,muhat,logdetSigmad2] = robust_vcov(S)
% Robust covariance matrix estimation as described in the METHODS section of the SL paper
% by Wood 2010.
% NOTE: quickly done, seems to work but not carefully tested

alpha=2; beeta=1.25;

% (i)
[Ns,Nr] = size(S);
dbar = sqrt(sum(S.^2,2)/Nr);
%Dbar = diag(dbar);
invDbar = diag(1./dbar);
[Qbar,Rbar] = qr((S'*invDbar)./sqrt(Nr-1));
%DRT = Dbar*Rbar';
%Sigmabar = DRT*DRT';
DiRi = invDbar/Rbar;
invSigmabar = DiRi*DiRi';
%Sigma = Sigmabar; invSigma = invSigmabar; muhat = []; logdetSigmad2 = [];

% (ii) & (iii)
muhat = mean(S')';
wjs = ones(1,Nr);
m0 = sqrt(Ns) + alpha/sqrt(2);
dx = bsxfun(@minus,S,muhat);
mjs = sum(dx'.*(invSigmabar*dx)',2);
ind = mjs > m0;
wjs(ind) = exp(-0.5*(mjs(ind) - m0).^2/beeta)*m0./mjs(ind);

% (iv)
muhat = sum(bsxfun(@times,S,wjs),2)/sum(wjs);
S = bsxfun(@minus,S,muhat(:));

% (v)
d = sqrt(sum(S.^2,2)/Nr);
D = diag(d);
W = diag(wjs);
invD = diag(1./d);
[Q,R] = qr((W*S'*invD)./sqrt(sum(wjs.^2)-1));

% (vi)
DRT = D*R';
Sigma = DRT*DRT';
E = (R')\invD;
invSigma = E'*E;
if size(R,2) == 1 % just one summary, needs to be handled separately
    logdetSigmad2 = log(abs(R(1))) + log(d);
else
    logdetSigmad2 = sum(log(abs(diag(R)))) + sum(log(d));
end
end




