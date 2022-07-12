function [sgRNA,x]=BLAST_Redux(sgRNA_scored,x,sgRNA,score)


%% Divide into 15mers
for a=1:4
    sgRNA_15mer(a,:)=sgRNA_scored(13:27);
    if a==1
    sgRNA_15mer(1,13)='A';
    elseif a==2
    sgRNA_15mer(2,13)='T';
    elseif a==3
    sgRNA_15mer(3,13)='C';
    else
    sgRNA_15mer(4,13)='G';
    end
end

%% Compare 15mer to gene
z=load('KM_17555.mat');
matrix_15mers=z.chromosome_matrix;
num=zeros(6,1);
for a=1:6
    dim=size(matrix_15mers{a});
    for b=1:4
        for c=1:dim(1)
            if strcmp(sgRNA_15mer(b,:),matrix_15mers{a}(c,:))==1
                num(a)=num(a)+1;
            end
        end
    end
end
num=sum(num);

%% Include
if num==1
    x=x+1;
    sgRNA{x,1}=sgRNA_scored;
    sgRNA{x,2}=score;
end
end