%PLOT_MARGINALS Plot marginals of infbench problem.

close all;

if ~exist('problem','var') || isempty(problem)
    problem = 'acerbi2012';
end

filename = [problem '_marginals.mat'];
load(filename);

N = numel(MarginalBounds);

for n = 1:N    
    if isempty(MarginalBounds{n}); continue; end
    figure(n);

    D = size(MarginalBounds{n},2);
    LB = MarginalBounds{n}(1,:);
    UB = MarginalBounds{n}(2,:);
    nx = size(MarginalPdf{n},2);

    for d = 1:D
        switch D
            case 1; subplot(1,1,d);
            case 2; subplot(1,2,d);
            case 3; subplot(1,3,d);
            case 4; subplot(2,2,d);
            case {5,6}; subplot(2,3,d);
            case {7,8,9}; subplot(3,3,d);        
        end
        xx = linspace(LB(d),UB(d),nx);
        pdf = MarginalPdf{n}(d,:);
        plot(xx, pdf, 'k-','LineWidth',2);
        xlim([LB(d),UB(d)]);
        set(gca,'TickDir','out');
        box off;
        xlabel(['x_' num2str(d)])
        ylabel('pdf');
    end
    set(gcf,'Color','w');

end