% This script is used to determine the position of every sgRNA in the genome
% along with identifiers such as the strand it targets and the chromosome it
% aligns to. This information is also critical for generating counts using
% inexact matching.

% A BAM file is generated by mapping the list of all sgRNA in the library
% (including nontargeting) to a genome fasta file (This should ideally have
% an extra chromosome that contain the nontargeting sgRNA). This file is
% read and sgRNA position, strand and chromosome information is recorded as
% mentioned in the following sections.

% Inputs: 
% 1. 'Library_sgRNA_Cas9.bam': An example BAM file obtained by running
% bowtie on the Cas9 library against the genome file
% 'CLIB89+control_Cas9.fasta'.
% 2. 'accessory_files_sgRNApos.mat': Due to certain errors such as similar
% names and sequences of sgRNA during library design, some accessory files
% were needed for properly determining sgRNA positions. They are all
% provided as this single matlab variable.

% Outputs: 
% All_sgRNA_pos_Cas9.mat: variable file that contains all sgRNA names,
% sequences and their start positions.


%Author: Adithya Ramesh
%PhD Candidate, Wheeldon Lab
%UC Riverside, 900 University Ave
%Riverside, CA-92507, USA
%Email: arame003@ucr.edu
%       wheeldon@ucr.edu
%% Reading BAM file and storing alignments against all chromosomes
%There will be an extra chromosome here as all the nontargeting sgRNA were
%added as a seperate chromosome to allow for the aligment of the entire
%library.
tic
bamFilename = 'Library_sgRNA_Cas9.bam'; %This file will be the result of running bowtie on the Cas9 library against the genome file 'CLIB89+control_Cas9.fasta' 
info = baminfo(bamFilename,'ScanDictionary',true);
bmA = BioMap(bamFilename, 'SelectReference', info.ScannedDictionary{1});
bmB = BioMap(bamFilename, 'SelectReference', info.ScannedDictionary{2});
bmC = BioMap(bamFilename, 'SelectReference', info.ScannedDictionary{3});
bmD = BioMap(bamFilename, 'SelectReference', info.ScannedDictionary{4});
bmE = BioMap(bamFilename, 'SelectReference', info.ScannedDictionary{5});
bmF = BioMap(bamFilename, 'SelectReference', info.ScannedDictionary{6});
if length(info.ScannedDictionary)==9
    bmH = BioMap(bamFilename, 'SelectReference', info.ScannedDictionary{8});
else
    bmH = BioMap(bamFilename, 'SelectReference', info.ScannedDictionary{7});
end
load('accessory_files_sgRNApos.mat');
%% Concatenate Start Positions, Strand information and Chromosome information

% Start positions of every sgRNA in the genome is recorded. Other
% information that is appended to the start position of an sgRNA is the
% chromosome it aligned to, and the strand that it was found on.
% For example: a start position of 2500018 is read this way: 25000|1|8.
% 8 refers to the chromosome, a flag of 0 or 1 mean Top and Bottom strands
% respectively, and the remaining is the actual start position of the
% sgRNA. So this sgRNA is found at the 25000 position on the bottom strand
% in the 8th chromosome.

Astart=bmA.Start;
Aflag=bmA.Flag;

Bstart=bmB.Start;
Bflag=bmB.Flag;

Cstart=bmC.Start;
Cflag=bmC.Flag;

Dstart=bmD.Start;
Dflag=bmD.Flag;

Estart=bmE.Start;
Eflag=bmE.Flag;

Fstart=bmF.Start;
Fflag=bmF.Flag;

Hstart=bmH.Start;
Hflag=bmH.Flag;

Aflag(Aflag==16)=1;
Bflag(Bflag==16)=1;
Cflag(Cflag==16)=1;
Dflag(Dflag==16)=1;
Eflag(Eflag==16)=1;
Fflag(Fflag==16)=1;
Hflag(Hflag==16)=1;

for i=1:length(Astart)
as=sprintf('%1d%1d%1d ',Astart(i),Aflag(i),1);
a1=str2num(as);
Afn(i,1)=a1;
end
for i=1:length(Bstart)
as=sprintf('%1d%1d%1d ',Bstart(i),Bflag(i),2); 
a1=str2num(as); 
Bfn(i,1)=a1;
end
for i=1:length(Cstart)
as=sprintf('%1d%1d%1d ',Cstart(i),Cflag(i),3); 
a1=str2num(as); 
Cfn(i,1)=a1;
end
for i=1:length(Dstart)
as=sprintf('%1d%1d%1d ',Dstart(i),Dflag(i),4); 
a1=str2num(as); 
Dfn(i,1)=a1;
end
for i=1:length(Estart)
as=sprintf('%1d%1d%1d ',Estart(i),Eflag(i),5); 
a1=str2num(as); 
Efn(i,1)=a1;
end
for i=1:length(Fstart)
as=sprintf('%1d%1d%1d ',Fstart(i),Fflag(i),6); 
a1=str2num(as); 
Ffn(i,1)=a1;
end
for i=1:length(Hstart)
as=sprintf('%1d%1d%1d ',Hstart(i),Hflag(i),8); 
a1=str2num(as); 
Hfn(i,1)=a1;
end

A(:,1)=bmA.Header;
A(:,2)=num2cell(Afn);
A(:,3)=bmA.Sequence;

B(:,1)=bmB.Header;
B(:,2)=num2cell(Bfn);
B(:,3)=bmB.Sequence;

C(:,1)=bmC.Header;
C(:,2)=num2cell(Cfn);
C(:,3)=bmC.Sequence;

D(:,1)=bmD.Header;
D(:,2)=num2cell(Dfn);
D(:,3)=bmD.Sequence;

E(:,1)=bmE.Header;
E(:,2)=num2cell(Efn);
E(:,3)=bmE.Sequence;

F(:,1)=bmF.Header;
F(:,2)=num2cell(Ffn);
F(:,3)=bmF.Sequence;

load('H.mat');
H(:,2)=num2cell(Hfn);
H(:,3)=bmH.Sequence;

sgAll=vertcat(A,B,C,D,E,F,H);
sgAll3=sgAll(:,3);
%% Add Start Position information to accessory information variables

for i=1:length(ChrA9)
    a=cell2mat(ChrA9(i,2));
    z=strfind(A(:,1),a);
    w=find(~cellfun(@isempty,z));
    if w>0
        ChrA9(i,4)=A(w,2);
    end
end


for i=1:length(ChrB9)
    a=cell2mat(ChrB9(i,2));
    z=strfind(B(:,1),a);
    w=find(~cellfun(@isempty,z));
    if w>0
        ChrB9(i,4)=B(w,2);
    end
end

for i=1:length(ChrC9)
    a=cell2mat(ChrC9(i,2));
    z=strfind(C(:,1),a);
    w=find(~cellfun(@isempty,z));
    if w>0
        ChrC9(i,4)=C(w,2);
    end
end

for i=1:length(ChrD9)
    a=cell2mat(ChrD9(i,2));
    z=strfind(D(:,1),a);
    w=find(~cellfun(@isempty,z));
    if w>0
        ChrD9(i,4)=D(w,2);
    end
end

for i=1:length(ChrE9)
    a=cell2mat(ChrE9(i,2));
    z=strfind(E(:,1),a);
    w=find(~cellfun(@isempty,z));
    if w>0
        ChrE9(i,4)=E(w,2);
    end
end

for i=1:length(ChrF9)
    a=cell2mat(ChrF9(i,2));
    z=strfind(F(:,1),a);
    w=find(~cellfun(@isempty,z));
    if w>0
        ChrF9(i,4)=F(w,2);
    end
end

for i=1:length(ChrExt9)
    a=cell2mat(ChrExt9(i,2));
    z=strfind(H(:,1),a);
    w=find(~cellfun(@isempty,z));
    if w>0
        ChrExt9(i,4)=H(w,2);
    end
end
%% Find any sgRNA for which position information might be missing

temp1=find(cellfun(@isempty,ChrA9(:,4)));
temp2=find(cellfun(@isempty,ChrB9(:,4)));
temp3=find(cellfun(@isempty,ChrC9(:,4)));
temp4=find(cellfun(@isempty,ChrD9(:,4)));
temp5=find(cellfun(@isempty,ChrE9(:,4)));
temp6=find(cellfun(@isempty,ChrF9(:,4)));
temp7=find(cellfun(@isempty,ChrExt9(:,4)));
Missing=vertcat(ChrExt9(temp7,1:3),ChrA9(temp1,1:3),ChrB9(temp2,1:3),ChrC9(temp3,1:3),...
    ChrD9(temp4,1:3),ChrE9(temp5,1:3),ChrF9(temp6,1:3));

for i=1:length(Missing)
    a=cell2mat(Missing(i,2));
    z=strfind(sgAll(:,1),a);
    w=find(~cellfun(@isempty,z));
    if w>0
        Missing(i,4)=sgAll(w,2);
    end
end

MissingAll=Missing;

ch=cell2mat(ChrExt9(:,1));
for i=1:length(ChrExt9)
    a=ch(i);
    mi=cell2mat(Missing(:,1));
    b=find(mi==a);
    if b>0
       Missing(b,:)=[];
    end
end
%% Remove missing/repeated sgRNA and tabulate the position for the final sgRNA list

All=vertcat(ChrExt9,Missing,ChrA9,ChrB9,ChrC9,ChrD9,ChrE9,ChrF9);
Ap=cell2mat(All(:,4));

[S,N]=count_unique(Ap);
test2=find(N>1);
test3=S(test2);

All(find(cell2mat(All(:,4))==test3(1)),:)=[];
All(find(cell2mat(All(:,4))==test3(2)),:)=[];
All(find(cell2mat(All(:,4))==test3(3)),:)=[];
All(find(cell2mat(All(:,4))==test3(4)),:)=[];
All(find(cell2mat(All(:,4))==test3(5)),:)=[];
All(find(cell2mat(All(:,4))==test3(6)),:)=[];
All(find(cell2mat(All(:,4))==test3(7)),:)=[];

All2=All(~cellfun('isempty',All(:,4)),:);
All_sgRNA_pos_Cas9=All2;
save All_sgRNA_pos_Cas9.mat All_sgRNA_pos_Cas9
toc