function [H,b,B] = mean_prep1(gp,x)
% Calculates the terms H, b and B needed in inference with mean function
% This is a modified version of a similar function in GPstuff. 

% Prepare variables
Hapu = cell(1,length(gp.meanf));
dimcount=0;
num_mf=length(gp.meanf);          % number of base functions used
meanf_dim=zeros(num_mf,2);        % col 1: base func nro i, col 2: amount of relevant input dimensions

for i=1:num_mf
    gpmf=gp.meanf{i};
    % base functions' values
    Hapu{i}=gpmf.fh.geth(gpmf,x);
    [dim,~] = size(Hapu{i});
    dimcount=dimcount+dim;          % amount of input dimensions total
    meanf_dim(i,:)=[i dim];
    
    
    % Gather prior mean for base functions into one vector
    if i==1                         % starting round?
        if dim==1
            b=gpmf.b;
        else
            if length(gpmf.b)==1
                b=repmat(gpmf.b,dim,1);
            else
                b=gpmf.b';
            end
        end
    else
        if dim==1
            b=cat(1,b,gpmf.b');
        else
            if length(gpmf.b)==1
                bvec=repmat(gpmf.b,dim,1);
                b=cat(1,b,bvec);
            else
                b=cat(1,b,gpmf.b');
            end
        end
    end
    
end

[nm mm]=size(meanf_dim);
% Gather base functions' values in one matrix
H = cat(1,Hapu{1:end});
% Gather prior covariances in one matrix B
B=zeros(dimcount,dimcount);
i1=1;
for i=1:nm
    if meanf_dim(i,2)==1
        if length(gp.meanf{i}.B)==1
            B(i1,i1)=gp.meanf{i}.B;
            i1=i1+1;
        else
            B(i1,:)=gp.meanf{i}.B;
            i1=i1+1;
        end
    else
        if length(gp.meanf{i}.B)>1 && length(gp.meanf{i}.B{i})==1
            for j=1:meanf_dim(i,2)
                B(i1+j-1,i1+j-1)=gp.meanf{i}.B{j};
                i1=i1+1;
            end
        elseif length(gp.meanf{i}.B)>1 && length(gp.meanf{i}.B{i})>1
            for j=1:meanf_dim(i,2)
                B(i1+j-1,:)=gp.meanf{i}.B{i};
            end
            i1=i1+meanf_dim(i,2);
        else
            for j=1:meanf_dim(i,2)
                B(i1+j-1,i1+j-1)=gp.meanf{i}.B;
            end
            i1=i1+meanf_dim(i,2);
        end
    end
end
end


