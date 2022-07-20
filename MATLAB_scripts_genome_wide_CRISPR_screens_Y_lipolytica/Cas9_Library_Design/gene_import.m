function [gene,chromosome]=gene_import()
x=input('Input number of genes to be imported ');

for a=1:x
    [myfilename, mypathname] = uigetfile('.fasta','Pick a file');
    myPathFile = [mypathname,myfilename];
    [header,sequence] =  fastaread(myPathFile);
    if ischar(sequence)==1
        seq=sequence;
        sequence=cell(1);
        sequence{1}=seq;
    end
    
    %% Read the sequences
    gene_cell=cell(length(sequence),1);
    for b=1:length(sequence)
        if length(sequence{1,b})<300
            gene_cell{b}=sequence{1,b};
        else
            gene_cell{b}=sequence{1,b}(1:300);
        end
    end
    gene{a,1}=gene_cell;
    
    %% Place headers into the cell
    name=cell(length(header),1);
    for b=1:length(header)
        if b<10
            name{b}=header{b}(44:56);
        elseif b>=10&&b<100
            name{b}=header{b}(45:57);
        else
            name{b}=header{b}(46:58);
        end
    end
    gene{a,2}=name;
end

%% Select Chromosome to use as database for run
y=input('For these genes, compare them to Chromosome: ','s');

if y=='a'
    chromosome='ChromosomeA.fasta';
elseif y=='A'
    chromosome='ChromosomeA.fasta';
elseif y=='b'
    chromosome='ChromosomeB.fasta';
elseif y=='B'
    chromosome='ChromosomeB.fasta';
elseif y=='c'
    chromosome='ChromosomeC.fasta';
elseif y=='C'
    chromosome='ChromosomeC.fasta';
elseif y=='d'
    chromosome='ChromosomeD.fasta';
elseif y=='D'
    chromosome='ChromosomeD.fasta';
elseif y=='e'
    chromosome='ChromosomeE.fasta';
elseif y=='E'
    chromosome='ChromosomeE.fasta';
elseif y=='f'
    chromosome='ChromosomeF.fasta';
elseif y=='F'
    chromosome='ChromosomeF.fasta';
end

end

