function Hs = mean_prep_xs(gp,xs)
% Calculates the term Hs needed in inference with mean function

num_mf=length(gp.meanf);          % number of base functions used
Hapu2 = cell(1,length(gp.meanf));

for i=1:num_mf
    gpmf=gp.meanf{i};
    Hapu2{i}=gpmf.fh.geth(gpmf,xs);
end

Hs = cat(1,Hapu2{1:end});
end


