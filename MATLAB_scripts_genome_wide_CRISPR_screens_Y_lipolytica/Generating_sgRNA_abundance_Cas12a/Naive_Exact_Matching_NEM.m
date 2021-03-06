% This script determines the naive exact matching counts of every sgRNAs in
% the final Cas12a sgRNA library, for a given NGS fastq read file. It does
% this by taking each sgRNA of the Cas12a library, and looking for its
% presence in the output file of Import_fastq.m.

% Input dependencies:
% 1. Run Import_fastq.m first (UniqueSeqs output variable is required).
% 2. The variable file All_sgRNA_pos_Cas12a.mat which also contains the list of sgRNA in
% the Cas12a sgRNA library is required. This file is included in the
% folder.

% Output Variables:
% NEM_counts_1: Gives the NEM counts of Cas12a library sgRNA for the NGS
% fastq file 'Examplefile1_Cas12a.fastqsanger'.


%Author: Adithya Ramesh
%PhD Candidate, Wheeldon Lab
%UC Riverside, 900 University Ave
%Riverside, CA-92507, USA
%Email: arame003@ucr.edu
%% Generate counts of sgRNA by Naive Exact Matching (NEM)
tic
clear
load('All_sgRNA_pos_Cas12a.mat'); %A cell type variable that contains all sgRNAs in the final Cas12a sgRNA library as a single column.
load('Counts_of_all_unique_reads1.mat'); %Output of Import_fastq.m
parfor i=1:length(All_sgRNA_pos_Cas12a)
    a=cell2mat(All_sgRNA_pos_Cas12a(i,2));
    z1=strfind(UniqueSeqs1,a);
    w1=find(~cellfun(@isempty,z1));
    if isempty(w1)
        ya=0;
    else
        ya=Counts1(w1);
    end
    
    y1=(sum(ya));
    NEM_counts_1(i,1)=y1;
end
save NEM_counts_1.mat NEM_counts_1
toc