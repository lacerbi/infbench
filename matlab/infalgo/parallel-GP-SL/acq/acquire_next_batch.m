function [ths,tm] = acquire_next_batch(th_grid, batch_size, iter, gp, y_tr, th_tr, ...
    sigma_tr, P, sim_model, gp_opt, acq_opt, other_opt, graphics_on)
% Select the next evaluation locations by optimising the acquisition function (aka design
% criterion).
% TODO: could remove duplicate code by refactoring, scaling of acq for better optimization?
%
% Acq methods implemented here:
% - unif                        (called 'RAND' in the paper)
% - MAXV   [seq, greedy]        (previously called EV)
% - MAXIQR [seq, greedy]        (previously called MIQR) 
% - EIV    [seq, joint, greedy] (previously called EIEV) 
% - IMIQR  [seq, joint, greedy] (previously called EIIQR/MIIQR) 

if iter < 1 || isempty(gp)
    error('Wrong iteration.');
end
if batch_size <= 0 || batch_size > 100
    error('Incorrect batch size.');
end
d = th_grid.dim;
interp_acq = 0;
vis_acq = (graphics_on >= 2);
seq_acq = (batch_size == 1);
greedy_batch = (batch_size > 1) && strcmp(acq_opt.batch_method,'greedy');
joint_batch = (batch_size > 1) && ~greedy_batch;
rep_th = (joint_batch==1)*batch_size+(joint_batch==0);
P_imiqr = []; f2 = [];

% record time used to acquire the next point(s) in the current iteration
start = tic;

if strcmpi(acq_opt.method,'unif')
    %% BASELINE METHOD: acquire uniformly randomly from the parameter space
    % ignores greedy batch option
    ths = acquire_init_batch(batch_size,th_grid);
    vis_acq = 0;

elseif strcmp(acq_opt.method,'MAXV') || strcmp(acq_opt.method,'EV')
    %% 'MAXVAR' AKA EXPONENTIATED VARIANCE/UNCERTAINTY SAMPLING (EVALUATE WHERE VARIANCE MAXIMISED) 
    if seq_acq
        % sequential case
        f = @(th_new) acq_MAXV(th_new, gp, th_tr, y_tr, P, sim_model);
        [ths,f_ths] = optimize_acq(f, th_grid, th_tr, acq_opt, rep_th);
        f = @(th_new) -acq_MAXV(th_new, gp, th_tr, y_tr, P, sim_model); % for plotting
        
    elseif greedy_batch
        % choose first point in a usual way and then always maximise expected acq function
        ths = zeros(batch_size,d);
        f_ths = zeros(batch_size,1);
        
        % 1) first point in batch:
        if strcmp(other_opt.display_type,'iter')
            disp(['greedy batch optimisation: 1/',num2str(batch_size)]);
        end
        f = @(th_new) acq_MAXV(th_new, gp, th_tr, y_tr, P, sim_model);
        [ths(1,:),f_ths(1)] = optimize_acq(f, th_grid, [], acq_opt, 1);
        
        % 2) later points:
        for i = 2:batch_size
            if strcmp(other_opt.display_type,'iter')
                disp(['greedy batch optimisation: ',num2str(i),'/',num2str(batch_size)]);
            end
            th_pend = ths(1:(i-1),:);
            P = acq_greedy_precompute(gp,th_tr,y_tr,sigma_tr,th_pend,P); % precompute
            fi = @(th_newi) acq_expected_V(th_pend, th_newi, gp, th_tr, y_tr, sigma_tr, P, sim_model);
            [ths(i,:),f_ths(i)] = optimize_acq(fi, th_grid, [], acq_opt, 1);
        end
        f = @(th_last) -acq_expected_V(th_pend, th_last, gp, th_tr, y_tr, sigma_tr, P, sim_model); % for plotting
        
    else
        % there is no 'non-greedy' MAXVAR method
        error('Not applicable.');
    end
    
elseif strcmp(acq_opt.method,'MAXIQR') || strcmp(acq_opt.method,'MIQR')
    %% IQR (EVALUATE WHERE THE LENGTH OF IQR IS MAXIMISED)
    if seq_acq
        % sequential case
        f = @(th_new) acq_MAXIQR(th_new, gp, th_tr, y_tr, P, sim_model);
        [ths,f_ths] = optimize_acq(f, th_grid, th_tr, acq_opt, rep_th);
        f = @(th_new) -acq_MAXIQR(th_new, gp, th_tr, y_tr, P, sim_model); % for plotting
        
    elseif greedy_batch
        % choose first point in a usual way and then always maximise median acq function
        ths = zeros(batch_size,d);
        f_ths = zeros(batch_size,1);
        
        % 1) first point in batch:
        if strcmp(other_opt.display_type,'iter')
            disp(['greedy batch optimisation: 1/',num2str(batch_size)]);
        end
        f = @(th_new) acq_MAXIQR(th_new, gp, th_tr, y_tr, P, sim_model);
        [ths(1,:),f_ths(1)] = optimize_acq(f, th_grid, [], acq_opt, 1);
        
        % 2) later points:
        for i = 2:batch_size
            if strcmp(other_opt.display_type,'iter')
                disp(['greedy batch optimisation: ',num2str(i),'/',num2str(batch_size)]);
            end
            th_pend = ths(1:(i-1),:);
            P = acq_greedy_precompute(gp,th_tr,y_tr,sigma_tr,th_pend,P); % precompute
            fi = @(th_newi) acq_median_iqr(th_pend, th_newi, gp, th_tr, y_tr, sigma_tr, P, sim_model);
            [ths(i,:),f_ths(i)] = optimize_acq(fi, th_grid, [], acq_opt, 1);
        end
        f = @(th_last) -acq_median_iqr(th_pend, th_last, gp, th_tr, y_tr, sigma_tr, P, sim_model); % for plotting
    else
        % there is no 'non-greedy' MAXIQR method
        error('Not applicable.');
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
elseif strcmp(acq_opt.method,'EIV') || strcmp(acq_opt.method,'EIEV')
    %% EXPECTED INTEGRATED VARIANCE (EVALUATE AT A POINT PRODUCING SMALLEST EXPECTED 
    %% INTEGRATED VARIANCE LOSS)
    if seq_acq || joint_batch
        % sequential or joint batch case
        [P_eiv,P] = exp_acq_precompute(th_grid, gp, th_tr, y_tr, P, sim_model, acq_opt);
        f = @(th_new) acq_eiv(th_new, gp, th_tr, y_tr, sigma_tr, P, P_eiv, acq_opt, 0);
        [ths,f_ths] = optimize_acq(f, th_grid, th_tr, acq_opt, rep_th);
        
    else
        % populate the batch in a greedy fashion
        % Choose first point in a usual way and then always maximise expected acq function
        ths = zeros(batch_size,d);
        f_ths = zeros(batch_size,1);
        
        % 1) first point in batch:
        if strcmp(other_opt.display_type,'iter')
            disp(['greedy batch optimisation: 1/',num2str(batch_size)]);
        end
        [P_eiv,P] = exp_acq_precompute(th_grid, gp, th_tr, y_tr, P, sim_model, acq_opt);
        f = @(th_new) acq_eiv(th_new, gp, th_tr, y_tr, sigma_tr, P, P_eiv, acq_opt, 0);
        [ths(1,:),f_ths(1)] = optimize_acq(f, th_grid, [], acq_opt, 1);
        
        % 2) later points:
        for i = 2:batch_size
            if strcmp(other_opt.display_type,'iter')
                disp(['greedy batch optimisation: ',num2str(i),'/',num2str(batch_size)]);
            end
            is_pend = 0;
            th_pend = ths(1:(i-1),:); % points kept 'fixed' while optimising the current one
            if 1
                % Additional GP precomputions to save some comput. time
                % In practice there doesn't seem to be much difference between this and 
                % naive computation, with small batch sizes at least
                is_pend = 1;
                P = exp_acq_greedy_precompute(gp, th_tr, y_tr, sigma_tr, P, P_eiv, th_pend, acq_opt);
            end
            th_pend = reshape(th_pend',1,[]);
            fi = @(th_new) acq_eiv([repmat(th_pend,size(th_new,1),1), th_new], gp, ...
                th_tr, y_tr, sigma_tr, P, P_eiv, acq_opt, is_pend);
            [ths(i,:),f_ths(i)] = optimize_acq(fi, th_grid, [], acq_opt, 1);
        end
        f = fi; % for plotting
    end
    interp_acq = 1;
    
elseif strcmp(acq_opt.method,'IMIQR') || strcmp(acq_opt.method,'EIIQR')
    %% INTEGRATED MEDIAN IQR (AS EIV BUT MORE ROBUST: IQR INSTEAD OF VARIANCE)
    %% (EVALUATE AT A POINT PRODUCING SMALLEST INTEGRATED MEDIAN LOSS)
    if seq_acq || joint_batch
        % sequential or joint batch case
        [P_imiqr,P] = exp_acq_precompute(th_grid, gp, th_tr, y_tr, P, sim_model, acq_opt);
        f = @(th_new) acq_imiqr(th_new, gp, th_tr, y_tr, sigma_tr, P, P_imiqr, acq_opt, 0);
        [ths,f_ths] = optimize_acq(f, th_grid, th_tr, acq_opt, rep_th);
        
        EXACT_MED_COMPUT = 0;
        if EXACT_MED_COMPUT
            f2 = @(th_new) acq_mediiqr_simul(th_new, gp, th_tr, y_tr, sigma_tr, P, P_imiqr, acq_opt);
        end
        
    else
        % Populate the batch in a greedy fashion
        % Choose first point in a usual way and then always maximise expected acq function
        ths = zeros(batch_size,d);
        f_ths = zeros(batch_size,1);
        
        % 1) first point in batch:
        if strcmp(other_opt.display_type,'iter')
            disp(['greedy batch optimisation: 1/',num2str(batch_size)]);
        end
        [P_imiqr,P] = exp_acq_precompute(th_grid, gp, th_tr, y_tr, P, sim_model, acq_opt);
        f = @(th_new) acq_imiqr(th_new, gp, th_tr, y_tr, sigma_tr, P, P_imiqr, acq_opt, 0);
        [ths(1,:),f_ths(1)] = optimize_acq(f, th_grid, [], acq_opt, 1);
        
        % 2) later points:
        for i = 2:batch_size
            if strcmp(other_opt.display_type,'iter')
                disp(['greedy batch optimisation: ',num2str(i),'/',num2str(batch_size)]);
            end
            is_pend = 0;
            th_pend = ths(1:(i-1),:); % points kept 'fixed' while optimising the current one
            if 1
                % Additional GP precomputions to save some comput. time
                % In practice there doesn't seem to be much difference between this and 
                % naive computation, with small batch sizes at least
                is_pend = 1;
                P = exp_acq_greedy_precompute(gp, th_tr, y_tr ,sigma_tr, P, P_imiqr, th_pend, acq_opt);
            end
            th_pend = reshape(th_pend',1,[]);
            fi = @(th_new) acq_imiqr([repmat(th_pend,size(th_new,1),1), th_new], gp, ...
                th_tr, y_tr, sigma_tr, P, P_imiqr, acq_opt, is_pend);
            [ths(i,:),f_ths(i)] = optimize_acq(fi, th_grid, [], acq_opt, 1);
        end
        f = fi; % for plotting
    end
    interp_acq = 1;
else
    error('Incorrect acq method.');
end

tm = toc(start);
if strcmp(other_opt.display_type,'iter')
    disp(' ');
    disp(['Time used for acquiring the next point(s) = ', num2str(tm), ' seconds.']);
    disp(' ');
end

%% make sure each row == one acquired location
if seq_acq || joint_batch
    ths = reshape(ths,d,[])';
end

%% Plot acq surface
if vis_acq
    viz_acq(f,f2,th_grid,f_ths,ths,joint_batch, batch_size, interp_acq, P_imiqr,th_tr, y_tr,...
        iter, acq_opt, other_opt, graphics_on);
end
end


function [] = viz_acq(f, f2, th_grid,f_ths,ths,joint_batch, batch_size, interp_acq, P_imiqr,...
    th_tr, y_tr, iter, acq_opt, other_opt, graphics_on)
% Wrapper for plotting the acq function surface

d = th_grid.dim;
acq_grid.opt_f = f_ths;
acq_grid.opt_th = ths;

if (graphics_on >= 2 && any(other_opt.viz_ind==iter)) || (graphics_on == 1 && iter == nr_iter)
    if ~joint_batch
        %% plotting in sequential or greedy batch case
        % acq surface is only to be plotted if dim <= 2
        if d == 1
            acq_grid.acq = f(th_grid.theta(:));
            if ~isempty(f2)
                acq_grid.acq_exact = f2(th_grid.theta(:));
            end
            if isfield(P_imiqr,'logIQR_cur')
                % compare current uncertainty value to the expected loss surface
                acq_grid.cur_acq = P_imiqr.logIQR_cur;
            end
        elseif d == 2
            % Plotting in full 2d grid may be costly for some acq functions so evaluate it in
            % a smaller 2d grid and then interpolate to the larger grid. Computing exactly
            % here is unnecessary because this computation is only for visualisation.
            
            if ~interp_acq
                acq_grid.acq = vec_to_grid_matrix(f(th_grid.theta2d'), th_grid);
                if ~isempty(f2)
                    acq_grid.acq_exact = vec_to_grid_matrix(f2(th_grid.theta2d'), th_grid);
                end
            else
                interp_method = 'cubic';
                grid_size = 50;
                th_s = theta_grid(th_grid.range, grid_size);
                acq_s = vec_to_grid_matrix(f(th_s.theta2d'), th_s);
                [th1_s,th2_s] = meshgrid(th_s.theta(1,:),th_s.theta(2,:));
                [th1,th2] = meshgrid(th_grid.theta(1,:),th_grid.theta(2,:));
                acq_grid.acq = interp2(th1_s, th2_s, acq_s, th1, th2, interp_method);
                if ~isempty(f2)
                    grid_size = 15;
                    th_s = theta_grid(th_grid.range, grid_size);
                    acq_s = vec_to_grid_matrix(f2(th_s.theta2d'), th_s);
                    [th1_s,th2_s] = meshgrid(th_s.theta(1,:),th_s.theta(2,:));
                    [th1,th2] = meshgrid(th_grid.theta(1,:),th_grid.theta(2,:));
                    acq_grid.acq_exact = interp2(th1_s, th2_s, acq_s, th1, th2, interp_method);
                end
            end
        else
            return;
        end
        figure(2);
        clf;
        plot_acq(th_grid, acq_grid, th_tr, y_tr, acq_opt);
        drawnow;
    else
        %% plotting in joint batch case
        if d == 1 && batch_size == 2
            % construct 2d grid from 1d parameter grid
            b_grid = theta_grid([th_grid.range(:)'; th_grid.range(:)'], 50);
            acq_grid.acq = vec_to_grid_matrix(f(b_grid.theta2d'), b_grid);
            
            figure(2);
            clf;
            plot_acq_batch(b_grid, acq_grid, th_tr, y_tr, acq_opt);
            drawnow;
        else
            return;
        end
    end
    
    %% save acq plot
    if other_opt.viz_save
        fn = [other_opt.output_folder,'/',other_opt.output_filename,'_acq_iter',num2str(iter)];
        my_export_fig(fn,'-transparent','-pdf');
    end
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ACQ: EVALUATE WHERE THE CURRENT MEASURE OF UNCERTAINTY IS HIGHEST

function nlogV = acq_MAXV(th_new, gp, th_tr, y_tr, P, sim_model)
% Evaluate the negative of the logarithm of the exponentiated variance (EV) acq function.

pr_val = sim_model.prior_eval(th_new);
[eft,varft] = gp_pred_fast(gp,th_tr,y_tr,th_new,P);

% nlogV = -2*log(pr_val) - (2*eft + varft) - log(exp(varft) - 1); % naive implementation
% % handle possible underflow of exp(varft)-term:
% nans = ~isfinite(nlogV);
% nlogV(nans) = -2*log(pr_val(nans)) - 2*(eft(nans) + varft(nans)); 

% robust implementation using a reformulated equation:
nlogV = -2*log(pr_val) - 2*(eft + varft) - log1p(-exp(-varft)); 

nlogV = real(nlogV); % ensure output is real valued just in case
if any(isnan(nlogV))
    warning('NaN values in acq function.');
	%keyboard;
end
end


function nlogIQR = acq_MAXIQR(th_new, gp, th_tr, y_tr, P, sim_model)
% Evaluate the negative of the logarithm of the IQR acq function.

pr_val = sim_model.prior_eval(th_new);
[eft,varft] = gp_pred_fast(gp,th_tr,y_tr,th_new,P);

iu = norminv(0.75);
il = norminv(0.25);
sft = sqrt(varft);
nlogIQR = -log(pr_val) - eft - iu*sft - log1p(-exp(sft*(il-iu))); 

nlogIQR = real(nlogIQR); % ensure output is real valued just in case
if any(isnan(nlogIQR))
    warning('NaN values in acq function.');
	%keyboard;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ACQ: EVALUATE WHERE THE CURRENT UNCERTAINTY IS HIGHEST - GREEDY BATCH VERSIONS

function P = acq_greedy_precompute(gp,th_tr,y_tr,sigma_tr,th_pend,P)
% Precompute some quantities for the greedy batch acquisition case where the median/expected
% value of the acquisition function is used to iteratively to fill the batch one point at
% a time. 

cov_new = gp_pred_cov_fast(gp,th_tr,y_tr,th_pend,th_pend,P,0);

% Sometimes the greedy variance method evaluates a lot near previous points and this appears to
% make the covariance matrix below almost non-positive definite. As a solution for this, we
% suppose new points are observed with negligible but higher than the initial default tol=
% 1e-9 noise value.
tols = [1e-9 1e-8 1e-7 1e-6 1e-5 1e-4 1e-3]; 
for tol = tols
    cov_new = cov_new + gp_noise_model_var(gp,th_tr,sigma_tr,1,tol) * eye(size(cov_new));
    [P.L_new,p] = chol(cov_new,'lower'); 
    if p == 0
        return; % chol succeeded
    end
end
error('Could not invert the L_new covariance matrix.');
end


function nlogexpV = acq_expected_V(th_pend, th_new, gp, th_tr, y_tr, sigma_tr, P, sim_model)
% Evaluate the negative of the logarithm of the expected variance acq function.
% 'th_pend' are the pending points i.e. points that have already been chosen to the batch
% in a greedy manner. Each row == one pending location. 
% NOTE: Here we cannot drop the first term because we are not integrating over the
% parameter space as in the corresponding EIV case. That is, the first term is not constant.

pr_val = sim_model.prior_eval(th_new);
[eft,varft] = gp_pred_fast(gp,th_tr,y_tr,th_new,P);

% Compute the deterministic change in GP variance
n_new = size(th_new,1);
tau_tpp = zeros(n_new,1);
for i = 1:n_new
    tau_tpp(i) = gp_lookahead_var_fast_pend(gp, varft(i), th_tr, y_tr, th_new(i,:), th_pend, P);
end
m_tpp = eft;
var_tpp = varft;

%nlogexpV = -2*log(pr_val) - 2*m_tpp + var_tpp - tau_tpp; % old, incorrect formula
nlogexpV = -2*log(pr_val) - 2*(m_tpp + var_tpp) - log1p(-exp(tau_tpp-var_tpp));

nlogexpV = real(nlogexpV); % ensure output is real valued just in case
if any(isnan(nlogexpV))
    warning('NaN values in acq function.');
	%keyboard;
end
end


function nlogmedIQR = acq_median_iqr(th_pend, th_new, gp, th_tr, y_tr, sigma_tr, P, sim_model)
% Evaluate the negative of the logarithm of the median (or expected) IQR acq function.
% 'th_pend' are the pending points i.e. points that have already been chosen to the batch
% in a greedy manner. Each row == one pending location. 

% NOTE: If MEDIAN == 0, then expected IQR is actually computed!
MEDIAN = 1; % This is exact because the median is taken (elementwise) and not of the integral formula

pr_val = sim_model.prior_eval(th_new);
[eft,varft] = gp_pred_fast(gp,th_tr,y_tr,th_new,P);

iu = norminv(0.75);
il = norminv(0.25);

% Compute the deterministic change in GP variance caused by the pending points
n_new = size(th_new,1);
tau_tpp = zeros(n_new,1); % change of var when pending points added to training data
s2_tpp = zeros(n_new,1); % var after pending points added to training data
for i = 1:n_new
    [tau_tpp(i),s2_tpp(i)] = gp_lookahead_var_fast_pend(gp, varft(i), th_tr, y_tr, ...
        th_new(i,:), th_pend, P);
end
s_tpp = sqrt(max(0,s2_tpp));
m_tpp = eft;

logtermi = log1p(-exp((il-iu)*s_tpp));
nlogmedIQR = -(log(pr_val) + m_tpp + (MEDIAN~=1)*0.5*tau_tpp + iu*s2_tpp + logtermi); % IQR

nlogmedIQR = real(nlogmedIQR); % ensure output is real valued just in case
if any(isnan(nlogmedIQR))
    warning('NaN values in acq function.');
    %keyboard;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ACQ: EVALUATE WHERE THE EXPECTED REDUCTION IN OVERALL UNCERTAINTY IS HIGHEST

function [P_acq,P] = exp_acq_precompute(th_grid, gp, th_tr, y_tr, P, sim_model, acq_opt)
% Precompute some terms for evaluating the EIV/MIIQR acq function. 

%% Set up the integration grid and weights
d = th_grid.dim;
if d <= 2
    % if dim <= 2, use a uniform grid
    if d == 1
        integr_grid = theta_grid(th_grid.range, acq_opt.exp.nr_grid.dim1);
        th_is = integr_grid.theta(:);
    else
        integr_grid = theta_grid(th_grid.range, acq_opt.exp.nr_grid.dim2);
        th_is = integr_grid.theta2d';
    end
    n_is = size(th_is,1); % number of integration grid points
    log_w_is = log(th_grid.dx*ones(n_is,1));
else
    % use importance sampling otherwise
    if strcmp(acq_opt.method,'EIV') || strcmp(acq_opt.method,'EIEV')
        f = @(th_new) acq_MAXV(th_new, gp, th_tr, y_tr, P, sim_model);
    elseif strcmp(acq_opt.method,'IMIQR') || strcmp(acq_opt.method,'EIIQR')
        f = @(th_new) acq_MAXIQR(th_new, gp, th_tr, y_tr, P, sim_model);
    else
        error('Incorrect acq method.');
    end
    [th_is, log_w_is1] = sample_from_acq(f, th_grid, th_tr, acq_opt);
    
    %% compute IS weights:
    m = max(log_w_is1);
    %log_w_is = log_w_is1 - m - log(sum(exp(log_w_is1 - m)));
    w_is1 = exp(log_w_is1 - m)/sum(exp(log_w_is1 - m)); % use logsumexp-trick
    %min_is = min(w_is1)
    
    % If the sampling fails, 0 values can occur here, this is an ad-hoc solution
    % for that. This is essentially the idea in the paper 'Truncated importance sampling'
    % by E. L. Ionides 2008
    is_tol = eps; 
    w_is1 = max(w_is1,is_tol);
    
    w_is = (1./w_is1)/(sum(1./w_is1));
    log_w_is = log(w_is);
    % ess = numel(w_is) / (1 + var(w_is))
    
    % FINITENESS CHECKING OF IS WEIGHTS
    if any(~isfinite(log_w_is)) || any(~isfinite(th_is(:)))
        %keyboard;
        error('IS in acq computation failed.');
    end
    
    if 0
        % simple debug plottings for inspecting IS weights
        figure(200); subplot(1,2,1); 
        histogram(w_is,40);
        if d == 1
            subplot(1,2,2); plot(th_is,w_is,'*k');
            xlim([th_grid.range(1),th_grid.range(2)]);
            ylim([0,0.05+max(w_is)]);
        elseif d == 2
            subplot(1,2,2); plot(th_is(:,1),th_is(:,2),'*k');
            xlim([th_grid.range(1,1), th_grid.range(1,2)]);
            ylim([th_grid.range(2,1), th_grid.range(2,2)]);
            set(gcf,'Position',[800 50 600 250]);
        else
            pts=th_is'
            %pause;
        end
        drawnow;
    end
end

%% Precompute
pr_val_is = sim_model.prior_eval(th_is);
% current GP mean and var at integration grid points:
[eft_is,varft_is] = gp_pred_fast(gp,th_tr,y_tr,th_is,P); 

if strcmp(acq_opt.method,'EIV') || strcmp(acq_opt.method,'EIEV')
    c_is = 2*eft_is + varft_is + 2*log(pr_val_is);
elseif strcmp(acq_opt.method,'IMIQR') || strcmp(acq_opt.method,'EIIQR')
    c_is = eft_is + log(pr_val_is); % IQR
    %c_is = 2*eft_is + 2*log(pr_val_is); % squared IQR
end
P_acq.log_pr_val = log(pr_val_is);
P_acq.eft_is = eft_is;
P_acq.varft_is = varft_is;
P_acq.th_is = th_is;
P_acq.lnw_is = log_w_is;
P_acq.c_is = c_is;

% Compute also current uncertainty value
if strcmp(acq_opt.method,'EIV') || strcmp(acq_opt.method,'EIEV')
    % TODO: Implement this
    %...
elseif strcmp(acq_opt.method,'IMIQR') || strcmp(acq_opt.method,'EIIQR')
    % This is similarly computed no matter which point is used of if expected/median
    % is used to select the next evaluation location(s)
    iu = norminv(0.75);
    s_cur = sqrt(P_acq.varft_is);
    logtermi = log1p(-exp(-2*iu*s_cur));
    v = P_acq.lnw_is + P_acq.c_is + iu*s_cur + logtermi; % IQR
    %v = P_acq.lnw_is + P_acq.c_is + 2*iu*s_cur + 2*logtermi; % squared IQR
    r_cur = max(v);
    P_acq.logIQR_cur = r_cur + log(sum(exp(v - r_cur)));
end

%% We also precompute some other GP terms.
% These are needed only for EIV/IMIQR acq function and are computed thus only here.
P.K_s1_tr = gp_cov(gp, th_is, th_tr, []);
P.K_s1_trK = (P.K_s1_tr/(P.L'))/P.L;
if isfield(gp,'meanf')
    P.HK_s1_trKT = P.H*P.K_s1_trK';
end
end


function [samples, acqs] = sample_from_acq(f, th_grid, th_tr, acq_opt)
% Generate samples from the acq surface (interpreted as a probability density function) 
% that is used as an IS proposal for computing the integral over the parameter space. 
% This code uses adaptive MCMC and is intended for dimensions > 2. If dim <= 2, then grid 
% approximation is used instead. 
% Other MCMC algorithms such as Hamiltonian Monte Carlo could be alternatively used but 
% HMC would require coding the gradients which would be somewhat tedious and might have 
% challenges because the acq surface can be (and often is) multimodal. 
%
% NOTE: f must a handle to function that computes the negative log of the acq surface 
% interpreted as a pdf. 

% Set up variables and options for MCMC 
display_type = acq_opt.display_type;
npar = th_grid.dim;
model.ssfun = @(th,data) 2*f(th);
data = [];
model.N = 1;

% The maximum acq value evaluated at training data points is taken as the initial point: 
% if the sampler can't move/gets trapped, this way we at least average over such a region 
% where the uncertainty is currently high.
f_all = f(th_tr);
[f_opt,opt_ind] = min(f_all); 
init = th_tr(opt_ind,:);

if ~strcmp(display_type,'off')
    disp('Initial point for MCMC:');
    init
end
params = cell(npar,1);
for i = 1:npar
    params{i} = {sprintf('\\theta_{%d}',i), init(i), th_grid.theta(i,1), th_grid.theta(i,end)};
end

% Additional MCMC settings
nfinal = acq_opt.exp.is_samples;
nchains = 5;
options.nsimu = 10000;
options.qcov = 1/10^2*diag((th_grid.theta(:,1)-th_grid.theta(:,end)).^2);
options.method = 'am';
options.updatesigma = 0;
options.verbosity = ~strcmp(display_type,'off'); % no printing from mcmc
options.waitbar = 0;
%options.burnintime = 1000; % what is this number actually here?

% Initialize results
samples_all = NaN(options.nsimu,npar,nchains);
acqs_all = NaN(options.nsimu,nchains);
results_all = cell(nchains,1);

% Run MCMC chains!
for i = 1:nchains
    if ~strcmp(display_type,'off')
        if i == 1
            disp('Running MCMC for sampling from acq...');
        end
        if nchains > 1
            disp(['Chain ', num2str(i), '/', num2str(nchains)]);
        end
    end
    [results,samples,~,acqs] = mcmcrun(model,data,params,options);
    results_all{i} = results;
    samples_all(:,:,i) = samples;
    acqs_all(:,i) = acqs;
    if i == nchains && ~strcmp(display_type,'off')
        disp('Done.');
    end
end

% Leave out burn-in (e.g. the first half of each chain)
cl = size(samples_all,1);
samples_all = samples_all(ceil(cl/2:cl),:,:);
acqs_all = acqs_all(ceil(cl/2:cl),:);

% Assess the convergence
% psrf is taken from GPstuff/diag
[R,neff,Vh,W,B,tau,thin] = psrf(samples_all);
mcmc_diag.R = R;
mcmc_diag.neff = neff;
mcmc_diag.is_converged = (max(abs(R-1)) < 0.1); % one number summary of convergence assessment
if ~strcmp(display_type,'off')
    disp(['nr chains = ', num2str(nchains)]);
    disp(['R = ',num2str(R)]);
    disp(['neff = ', num2str(neff)]);
end
if mcmc_diag.is_converged ~= 1
    warning('Convergence not reached when sampling from the acq pdf.');
end

% Add the chains together
cl = size(samples_all,1);
samples = NaN(nchains*cl,npar);
for i = 1:nchains
    samples(1 + (i-1)*cl:i*cl,:) = samples_all(:,:,i);
end
acqs = acqs_all(:);

% Thin to final length (to reduce the save size of the samples-matrix)
final_cl = min(nfinal,size(samples,1));
final_ind = floor(linspace(1,size(samples,1),final_cl));
samples = samples(final_ind,:);
acqs = acqs(final_ind);

% Print and plot results to visually examine whether convergence is reached
if ~strcmp(display_type,'off')
    figure(100);
    clf;
    mcmcplot(samples,1:npar,results_all{1}.names,'chainpanel');
    suptitle('acq pdf');
    set(gcf,'Position',[60 600 600 400]);
    drawnow;
end

% Finally, transform acq values back to the original domain (log pdf)
acqs = -0.5*acqs;
end


function P = exp_acq_greedy_precompute(gp, th_tr, y_tr, sigma_tr, P, P_acq, th_pend, acq_opt)
% Precompute some extra terms needed for evaluating the greedy batch version of the 
% EIV/IMIQR acq function. 
% These are needed on lines 52...55 of the function 'gp_lookahead_var_fast_pend'.
% Note: We could alternatively update the GP and then use it for the computations.

P.d_var_pend = gp_lookahead_var_fast(gp,NaN,th_tr,y_tr, sigma_tr,P_acq.th_is,th_pend,P,0);
cov_pend = gp_pred_cov_fast(gp,th_tr,y_tr,th_pend,th_pend,P,0);
cov_pend = cov_pend + gp_noise_model_var(gp, th_tr, sigma_tr) * eye(size(cov_pend));
P.Lpend = chol(cov_pend,'lower');
P.cov_spend = gp_pred_cov_fast(gp,th_tr,y_tr,P_acq.th_is,th_pend,P,1); 
end


function logEIV = acq_eiv(th_new, gp, th_tr, y_tr, sigma_tr, P, P_eiv, acq_opt, is_pend)
% Evaluate the expected integrated variance (EIV) acq function.
% IMPORTANT NOTE: We actually drop the constant term and compute only the second term that
% has negative sign in front of it. 
% NOTE: In batch case each row of 'th_new' is a concatenated vector of candidate batch
% locations with length dim of param space * batch_size

% Compute the deterministic change in GP variance
d = size(th_tr,2);
[n_new,fd] = size(th_new);
n_is = length(P_eiv.lnw_is);
tau = zeros(n_is,n_new);
for i = 1:n_new
    tau(:,i) = gp_lookahead_var_fast(gp, P_eiv.varft_is, th_tr, y_tr, sigma_tr, P_eiv.th_is, ...
        reshape(th_new(i,:),d,fd/d)', P, is_pend);
end

% For each candidate point in th_new, compute the log acq function value
% Use logsumexp trick for stable log computations
logEIV = zeros(n_new,1);
for i = 1:n_new
    v_i = P_eiv.lnw_is + P_eiv.c_is + tau(:,i);
    %logEIV(i) = log(sum(exp(v_i))); % naive implementation
    r_newi = max(v_i);
    logEIV(i) = r_newi + log(sum(exp(v_i - r_newi)));
end
logEIV = -logEIV;

logEIV = real(logEIV); % ensure output is real valued just in case
if any(isnan(logEIV))
    warning('NaN values in acq function.');
    %keyboard;
end
end


function logIMIQR = acq_imiqr(th_new, gp, th_tr, y_tr, sigma_tr, P, P_imiqr, acq_opt, is_pend)
% Evaluate the integrated median IQR (IMIQR) acq function.
% NOTE: We compute the whole quantity (also the const part) and its value is then to be minimised. 
% NOTE: In batch case each row of 'th_new' is a concatenated vector of candidate batch
% locations with length dim of param space * batch_size

% testing: compute by simulation the MIIQR
%if 1
%    logMIIQR = acq_mediiqr_simul(th_new, gp, th_tr, y_tr, sigma_tr, P, P_miiqr, acq_opt);
%    logIMIQR = logMIIQR; 
%    return;
%end

APPROX_MEDIAN = 1; % if 1, we compute IMIQR, otherwise EIIQR

% Compute the deterministic change in GP variance
d = size(th_tr,2);
[n_new,fd] = size(th_new);
n_is = length(P_imiqr.lnw_is);
tau = zeros(n_is,n_new);
s2_new = zeros(n_is,n_new);
for i = 1:n_new
    [tau(:,i),s2_new(:,i)] = gp_lookahead_var_fast(gp, P_imiqr.varft_is, th_tr, y_tr, sigma_tr,...
        P_imiqr.th_is, reshape(th_new(i,:),d,fd/d)', P, is_pend);
end
s_new = sqrt(max(0,s2_new));

% For each candidate point in th_new, compute the log acq function value
% Use logsumexp trick for stable log computations
iu = norminv(0.75); 
logIMIQR = zeros(n_new,1);
for i = 1:n_new
    %logtermi = real(log1p(-exp(-2*iu*s_new(:,i))));
    logtermi = (log1p(-exp(-2*iu*s_new(:,i))));
    v_i = P_imiqr.lnw_is + P_imiqr.c_is + (APPROX_MEDIAN~=1)*0.5*tau(:,i) + iu*s_new(:,i) + logtermi; % IQR
    %v_i = P_miiqr.lnw_is + P_miiqr.c_is + (APPROX_MEDIAN~=1)*2*tau(:,i) + 2*iu*s_new(:,i) + 2*logtermi; % squared IQR
    %v_i = P_miiqr.lnw_is + P_miiqr.c_is + 0.5*sqrt(tau(:,i)) + iu*s_new(:,i) + logtermi; % TESTING, extra sqrt(tau)-term
    r_newi = max(v_i);
    logIMIQR(i) = r_newi + log(sum(exp(v_i - r_newi)));
end

logIMIQR = real(logIMIQR); % ensure output is real valued just in case
if any(isnan(logIMIQR))
    warning('NaN values in acq function.');
    %keyboard;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function log_MIIQR = acq_mediiqr_simul(th_new, gp, th_tr, y_tr, sigma_tr, P, P_miiqr, acq_opt)
% Computes the exact median integrated loss MIIQR instead of IMIQR using costly simulations 
% from the GP. Works in 1d or 2d case.
%
% NOTE: This is a quickly made implementation intended only for testing purposes and could
% be improved various ways. Computations would nevertheless remain very costly and impractical.

d = size(th_tr,2);
n_new = size(th_new,1);
if d > 2
    error('Not implemented, needs tight discretisation.');
end
n_is = length(P_miiqr.lnw_is);

% how many samples from the GP to approx. the density of the integral:
if d==1
    N = 5000;
else
    N = 2000;
end
iu = norminv(0.75); % NOTE: here we consider the case with iu == -il

log_MIIQR = zeros(n_new,1);
for i = 1:n_new
    taucov = gp_lookahead_cov_fast(gp, th_tr, y_tr, sigma_tr, P_miiqr.th_is, th_new(i,:), P);
    s2_new = P_miiqr.varft_is - diag(taucov);
    s_new = sqrt(max(0,s2_new));
    
    % sample potential m_(t+1) value evaluated in the grid from the GP model
    mt = P_miiqr.eft_is; % current GP mean function
    if d==1
        % often taucov is not pd but only psd, this is handled gracefully here: mvnrnd
        % uses cholcov which uses eigendecomposition whenever standard chol fails
        % but this can be slow for large taucov matrix!
        mtpps_j = mvnrnd(mt,taucov,N); 
    else
        % this is much faster than above alternative, but requires some jitter to make the
        % matrix psd so that standard cholesky can be used.
        try
            jitter = 1e-7;
            tauL = chol(taucov + jitter*eye(size(taucov)),'lower');
        catch
            keyboard;
        end
        mtpps_j = bsxfun(@plus, mt, tauL*randn(n_is,N))';
    end
    
    ilv = NaN(N,1); % integrated losses saved here
    logtermi = log1p(-exp(-2*iu*s_new)); 
    c_i = P_miiqr.lnw_is + P_miiqr.log_pr_val + iu*s_new + logtermi;
    for j = 1:N
        % TODO: over/underflow possible, could use logsumexp trick:
        vals = mtpps_j(j,:) + c_i';
        ilv(j) = sum(exp(vals));
    end
    log_MIIQR(i) = log(median(ilv));
    %log_MIIQR(i) = log(quantile(ilv,0.1));
    %log_MIIQR(i) = log(mean(ilv));
end

log_MIIQR = real(log_MIIQR); % ensure output is real valued just in case
if any(isnan(log_MIIQR))
    keyboard;
end
end




