function ordered_score=score_ranking(genes,score)
%This function receives 2 inputs and combines them into a cell matrix. 
%It then orders them from highest to lowest by the second column.
format compact

%% Create Cells
y=size(genes);
sgRNA=cell(y(1),2);
for a=1:y(1)
    sgRNA{a,1}=genes(a,:);
end

for c=1:length(score)
    sgRNA{c,2}=score(c);
end

%% Sort 
sgRNA_sorted=sortrows(sgRNA,2);
if y(1)>1
    ordered_score=flip(sgRNA_sorted);
else
    ordered_score=sgRNA_sorted;
end
end