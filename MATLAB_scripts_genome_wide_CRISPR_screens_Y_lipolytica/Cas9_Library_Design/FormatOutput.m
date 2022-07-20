function Output=FormatOutput(name,sgRNA)
dim=size(sgRNA);
for a=1:dim(1);
    Output{a,2}=name;
    x=sgRNA{a,1};
    if x==1
        Output{a,1}=[];
    elseif isempty(sgRNA{a,1})==1
        Output{a,1}=[];
    else
        Output{a,1}=x(5:24);
    end
    Output{a,3}=sgRNA{a,2};
end
end