function [] = demo_gp_delta_s2()
% Plots *Figure 2* in the paper.

subplot = @(m,n,p) subtightplot (m, n, p, [0.07 0.07], [0.2 0.1], [0.1 0.1]); % better subplots
close all;

root = '../results/demo_delta_s2/';
n_grid = 1000;
x_grid = linspace(-4,12,n_grid)';

for i = 1:2
    if i == 1
        n = 3;
        x_tr0 = -200+[0,1,2]'; % code doesn't work without some initial evals so avoid the difficulty this way
        x_tr1 = 0;
        x_tr2 = 0.8;
        ms2 = 5^2;
        l = 2.5;
        sn2 = .005^2;
    elseif i == 2
        n = 4;
        x_tr0 = [-1,1,2, -2]';
        x_tr1 = 0.25; 
        x_tr2 = 4; 
        ms2 = 5^2;
        l = 2;
        sn2 = 2^2;
    end
    y_tr0 = ones(n,1); % these are redundant actually
    
    %% construct GP
    gpcf = gpcf_sexp('lengthScale', l, 'magnSigma2', ms2, ...
        'lengthScale_prior', prior_fixed(), 'magnSigma2_prior', prior_fixed());
    lik = lik_gaussian('sigma2', sn2, 'sigma2_prior', prior_fixed());
    gp_opt = [];
    gp = gp_set('lik', lik, 'cf', gpcf, 'jitterSigma2', 1e-9);
    
    %% compute variances
    % current gp variance surface at time 0
    P = precompute_gp_pred(gp, x_tr0, y_tr0, gp_opt);
    [~,s20] = gp_pred_fast(gp, x_tr0, y_tr0, x_grid, P);
    
    % gp variance surface when obs. 1 added
    P.K_s1_tr = gp_cov(gp, x_grid, x_tr0, []);
    P.K_s1_trK = (P.K_s1_tr/(P.L'))/P.L;
    [ds2_1,s2_1] = gp_lookahead_var_fast(gp,s20,x_tr0,y_tr0,[],x_grid,x_tr1,P,0);
    
    % gp variance surface when obs. 2 added
    [ds2_2,s2_2] = gp_lookahead_var_fast(gp,s20,x_tr0,y_tr0,[],x_grid,x_tr2,P,0);
    
    % gp variance surface when both observations added
    [ds2_12,s2_12] = gp_lookahead_var_fast(gp,s20,x_tr0,y_tr0,[],x_grid,[x_tr1;x_tr2],P,0);
    
    % gp variance surface when both observations added and handled SEPARATELY i.e. as if
    % there would be no correlations between the points
    s2_12_sep = max(0, s20 - (ds2_1 + ds2_2));
    
    
    jitter = 1; % add some jitter so that the curves can be better seen
    dj = 0.03;
    lw = 1; % linewidth, matlab default 0.5
    
    %% plot
    figure(1);
    subplot(1,2,i);
    hold on;
    % lines:
    y_col = [93 69 12]/100; % better yellow color
    h1(5) = plot(x_grid, s2_12_sep, 'color', y_col,'Linewidth',lw);
    h1(1) = plot(x_grid, s20 + jitter*4*dj, '-k','Linewidth',lw);
    h1(4) = plot(x_grid, s2_12 + jitter*1*dj, '-b','Linewidth',lw);
    h1(2) = plot(x_grid, s2_1 + jitter*2*dj, '-r','Linewidth',lw);
    h1(3) = plot(x_grid, s2_2 + jitter*3*dj, '--r','Linewidth',lw);
    lm = [0,max(s20)*1.05];
    ylim(lm);
    
    % gray vertical lines and evaluation locations:
    gray_col = 0.75*[1 1 1]; % gray color
    plot(x_tr1*[1 1], lm, 'color',gray_col,'Linewidth',lw);
    plot(x_tr2*[1 1], lm, 'color',gray_col,'Linewidth',lw);
    plot(x_tr1, lm(1), '*k' ,'MarkerSize',7,'Linewidth',1.2); % obs. 1
    plot(x_tr2, lm(1), '*k' ,'MarkerSize',7,'Linewidth',1.2); % obs. 2
    plot(x_tr0, lm(1)*ones(size(x_tr0)), 'xk','MarkerSize',7,'Linewidth',1.2); % initial observations
    xlim([min(x_grid),max(x_grid)]);
    hold off;
    ylabel('GP variance');
    xlabel('\theta');
    box on;
    if i == 1
        legend(h1,{'s^2_0(\theta)',...
            's^2_1(\theta;\theta_1^{\ast})',...
            's^2_1(\theta;\theta_2^{\ast})',...
            's^2_{1:2}(\theta;\theta_{1:2}^{\ast})',...
            's_0^2(\theta){}-{}\tau_0^2(\theta;\theta_1^{\ast}){}-{}\tau_0^2(\theta;\theta_2^{\ast})'...
            },'Location','Southeast');
        %'s_0^2(\theta){}-{}\Delta{}s_0^2(\theta;\theta_1^{\ast}){}-{}\Delta{}s_0^2(\theta;\theta_2^{\ast})'...
    end
    set(gcf,'Position',[50 50 1000 340]);
    if i == 1
        title('(a)');
    elseif i == 2
        title('(b)');
    end
end

% save the generated figure
if 1
    fn = [root,'fig_s2_gp'];
    my_export_fig(fn,'-transparent','-pdf');
end
end




