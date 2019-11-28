function [] = demo_sl_gaussianity()
% Plots *Figure 10* in the paper. 
%
% Run each benchmark simulation model several times to examine the 'Gaussianity' of
% the resulting synthetic loglik-values.
% Computing takes ~5...30minutes depending on N values and amount of repeats.
% NOTE: There is also a similar function explicitly designed for ricker case: 'test_ricker' 

clc; close all; format('compact');
rng(1234567);

subplot = @(m,n,p) subtightplot (m, n, p, [0.1 0.05], [0.1 0.1], [0.1 0.1]); % better subplots

% Use the following test models, one realisation of data and the true parameter value. 
% Of course, it could be that Gaussianity assumptions are heavily violated only in some
% specific parts of the parameter space or with some specific data etc. but these all 
% cases cannot be tested here reasonably. 

test_models = {'gaussian2d','ricker','gk_model','ma2','lorenz','cell_model'};
%Ns = [100, 100, 50, 20, 50, 500]; % OLD, INITIAL
Ns = [50, 80, 20, 20, 50, 500];
titles = {'Gaussian','Ricker','g-and-k','Moving average','Lorenz','Cell biology'};

reps = 300; % how many repeated set of N summaries for each test model
sl_method = 'sl'; %'sl','ubsl','ublogsl'

load_old_res = 1; % LOAD ALREADY COMPUTED SUMMARIES
root = '../results/sl_gaussianity/'; % results saved here

% get current test problem settings etc.
nmodels = length(test_models);
sim_models = cell(nmodels,1);
for t = 1:nmodels
    [~,sim_models{t}] = get_test_model(test_models{t},[],Ns(t));
end

%% 1) Compute summaries for each model considered (this is comput. costly)
tic;
sim_summaries = cell(nmodels,1);
for t = 1:nmodels
    data_name = [root,'/summary_data_',test_models{t},'_',num2str(Ns(t)),'.mat'];
    if load_old_res
        % already computed, load it
        load(data_name); % summaries_th
        sim_summaries{t} = summaries_th;
    else 
        summaries_th = NaN(Ns(t),reps,sim_models{t}.summary_dim); % N x reps x summary dim
        
        % repeated sampling for SL (this could be parallelized)
        theta = sim_models{t}.true_theta;
        for k = 1:reps
            for j = 1:Ns(t)
                data_j = sim_models{t}.gen(theta, sim_models{t}.n_data);
                summaries_th(j,k,:) = sim_models{t}.comp_summaries(data_j, sim_models{t}.data);
            end
        end
        sim_summaries{t} = summaries_th;
        save(data_name, 'summaries_th');
        disp(['Simulations done for model ', num2str(t), '/', num2str(nmodels)]);
    end
end
toc


%% 2) Compute SL values using the already computed summaries
for t = 1:nmodels
    %% Extract values
    t
    sl_vals = NaN(reps,1);
    for j = 1:reps % repeats
        simul_summs = squeeze(sim_summaries{t}(:,j,:))';
        sl_vals(j) = eval_sl_loglik(sim_models{t}.summary_true, simul_summs, sl_method);
    end
    if any(~isfinite(sl_vals))
        disp('Nonfinite values in SL computation occurred.');
    end
    
    sls_all = sl_vals;
    sls_vals = sls_all(isfinite(sls_all)); % ignore possible NaNs
    if length(sls_all) > length(sls_vals) % check how many NaNs was computed
        disp(num2str(length(sls_all) - length(sls_vals)));
    end
    
    me = mean(sls_vals);
    va = var(sls_vals);
    stdev = sqrt(va);
    mi = min(sls_vals);
    ma = max(sls_vals);
    
    nbins = 25;
    
    figure(1);
    subplot(2,ceil(nmodels/2),t);
    hold on;
    h=histogram(sls_vals,nbins,'Normalization','pdf');
    % PLOT GAUSSIAN
    lims = xlim;
    x = linspace(lims(1),lims(2),1000);
    gauss = exp(-(x-me).^2./(2*stdev^2))./(stdev*sqrt(2*pi));
    plot(x,gauss,'LineWidth',1.5);
    hold off;
    box on;
    stdev_r = sprintf('%.2f',stdev);
    %title([titles{t}, ', stdev: ',stdev_r], 'interpreter', 'none');
    title([titles{t}, ' (stdev: ',stdev_r, ', N: ', num2str(Ns(t)),')'], 'interpreter', 'none');
    set(gcf,'Position',[50 500 ceil(nmodels/2)*400 500]);
    %set(h,'FaceAlpha',1);
    xlim([mi-0.5,ma+0.5]);
    drawnow;
    
    %% save the histogram(s) to file
    if t == nmodels
        img_name = [root, '/hist_all_', sl_method];
        %my_export_fig(img_name,'-transparent','-pdf'); 
        my_export_fig(img_name,'-dpdf','-pdf');
        % doesn't seem to work properly, as a fix, use the code below and manual cropping
        
        img_name2 = [root, '/hist_all2_', sl_method];
        set(gcf, 'Units', 'Inches');
        pos = get(gcf, 'Position');
        set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);
        print(gcf, img_name2, '-dpdf','-r0');
    end
end




