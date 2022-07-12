function sgRNA=divider_30mer(gene)

%% Divide into the 30mers
b=length(gene)-29; %Number of 30mers
genes_matrix=char(zeros((2*b),30));
for a=1:b
    genes_matrix(a,:)=gene(a:a+29);
end

%% Add reverse gene sequences
for a=1:b
    x=seqrcomplement(genes_matrix(a,:));
    genes_matrix(a+b,:)=x;
end
% It would be possible to specify which strand the scores came off of by
% using a for loop.

%% Check to see if sgRNA
y=size(genes_matrix);
b=0;
sgRNA=char(zeros(1,30));
for a=1:y(1)
    if genes_matrix(a,26)=='G'
        if genes_matrix(a,27)=='G'
            b=b+1;
            sgRNA(b,:)=genes_matrix(a,:);
        end
    end
end 

if double(sgRNA(1,1))==0
    for a=1:30
        sgRNA(1,a)='X';
    end
end
end