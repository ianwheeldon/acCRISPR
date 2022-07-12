 %% Cory's CRISPR Program

%% Admin

clear all
close all
clc
format compact

%% Determine test type
valz=0;
while valz==0
    answer=input('Are you running a nucleotide sequence or an amino sequence? ','s');
    if strcmp(answer,'nucleotide')==0
        if strcmp(answer,'amino')==0
            disp('Enter valid answer - "nucleotide" or "amino" ')
        else
            valz=1;
        end
    else 
        valz=1;
    end
end
        
%% Call Genome
if strcmp(answer,'nucleotide')==1
    [gene,chromosome]=gene_import();
elseif strcmp(answer,'amino')==1
    [gene,chromosome]=amino_import();
end

%% Analyze Genes
dimgene=size(gene);
number_complete=3
range=1;
if strcmp(answer,'nucleotide')==1
    for a=1:dimgene(1)
        for b=number_complete+1:number_complete+range
            [Output{a}{b},blastscore{a}{b}]=Gene_Analysis(gene{a,1}{b,1},gene{a,2}{b,1},chromosome);
            number_complete=number_complete+1
        end
    end
elseif strcmp(answer,'amino')==1
    for a=1:dimgene(1)
        for b=number_complete+1:number_complete+range 
            [Output{a}{b},blastscore{a}{b}]=Gene_Analysis_Amino(gene{a,1}{b,1},gene{a,2}{b,1},chromosome);
            number_complete=number_complete+1
        end
    end
end

%% Output to Excel
sgRNA_to_Excel(Output,blastscore)
winopen('Analyzed Data.xlsx')