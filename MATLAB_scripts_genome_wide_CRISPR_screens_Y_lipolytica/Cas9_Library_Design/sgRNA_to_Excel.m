function []=sgRNA_to_Excel(Output,blastscore)

%% Format to write in Excel
n=0;
[x,s]=size(Output);
for b=1:x
    [i j]=size(Output{1,b});
    for a=1:j
        [i2 j2]=size(Output{1,b}{1,a});
        for c=1:i2
            Outputs2{n+a,1}=Output{1,b}{1,a}{c,1};
            Outputs2{n+a,2}=Output{1,b}{1,a}{c,2};
            Outputs2{n+a,3}=Output{1,b}{1,a}{c,3};
            [n m]=size(Outputs2);
        end
    end
end

Outputs3=cell(1,3);
[i j]=size(Outputs2);
n=1;
for a=1:i
    if ischar(Outputs2{a,2})==1
        Outputs3{n,1}=Outputs2{a,1};
        Outputs3{n,2}=Outputs2{a,2};
        Outputs3{n,3}=Outputs2{a,3};
        n=n+1;
    end
end

n=1;
[x,s]=size(blastscore);
for b=1:s
    [i j]=size(blastscore{b});
    for a=1:j
        blastscore2{n}=blastscore{b}{a};
        n=n+1;
    end
end

[i j]=size(blastscore2);
n=1;
for a=1:j
    if isempty(blastscore2{1,a})==0
        blastscore3{n}=blastscore2{1,a};
        n=n+1;
    end
end

Outputs=cell(1000,10);
[i j]=size(Outputs3);
n=1;
for a=1:i
    if strcmp(Outputs{n,1},Outputs3{a,2})==1
        Outputs{n,x+1}=Outputs3{a,1};
        Outputs{n,x+2}=Outputs3{a,3};
        x=x+2;
    else
        if a~=1
            n=n+1;
        end
        Outputs{n,1}=Outputs3{a,2};
        Outputs{n,2}=blastscore3{n};
        Outputs{n,3}=Outputs3{a,1};
        Outputs{n,4}=Outputs3{a,3};
        x=4;
    end
end
        
%% Output into Excel
xlswrite('Analyzed Data.xlsx',Outputs)

end