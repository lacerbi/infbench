function MTV = ComputeMarginalTotalVariation(xx,probstruct)
%COMPUTEMARGINALTOTALVARIATION Compute MTV between samples and PROBSTRUCT.

D = size(xx,2);
nkde = 2^13;
MTV = zeros(1,D);

% Compute marginal total variation
for i = 1:D
    [~,yy1,xmesh] = kde1d(xx(:,i),nkde);
    yy1 = yy1/(qtrapz(yy1)*(xmesh(2)-xmesh(1))); % Ensure normalization
    yy_true = probstruct.Post.MarginalPdf(i,:);
    xx_true = linspace(probstruct.Post.MarginalBounds(1,i),probstruct.Post.MarginalBounds(2,i),size(yy_true,2));    
    f = @(x) abs(interp1(xmesh,yy1,x,'spline',0) - interp1(xx_true,yy_true,x,'spline',0));    
    bb = sort([xx_true([1,end]),xmesh([1,end])]);
    for j = 1:3
        xx_range = linspace(bb(j),bb(j+1),1e6);
        MTV(i) = MTV(i) + 0.5*qtrapz(f(xx_range))*(xx_range(2)-xx_range(1));
    end
end

end