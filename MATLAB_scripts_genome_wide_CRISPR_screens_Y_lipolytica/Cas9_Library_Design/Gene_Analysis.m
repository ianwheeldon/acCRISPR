function [Output,blastscore]=Gene_Analysis(gene,name,chromosome)

%% Blast against selected chromosome
[gene,blastscore]=initial_blast(gene,chromosome);

%% Seperate into 30mers
genes=divider_30mer(gene);

%% Rate each Gene
a=size(genes);
score=NaN(a(1),1);
for b=1:a(1)
    score(b,1)=on_target_score_calculator(genes(b,:));
end

%% Combine scores rank them
sgRNA_scored=score_ranking(genes,score);

%% BLAST
n=0;
x=0;
[i j]=size(sgRNA_scored);
sgRNA_number=4; %This determines how many outputs there are per gene
sgRNA=cell(sgRNA_number,2);
while x<sgRNA_number
    n=n+1;
    [sgRNA,x]=BLAST_Redux(sgRNA_scored{n,1},x,sgRNA,sgRNA_scored{n,2});
    if n==i
        x=sgRNA_number;
    end
end

%% Format Output
Output=FormatOutput(name,sgRNA);
end