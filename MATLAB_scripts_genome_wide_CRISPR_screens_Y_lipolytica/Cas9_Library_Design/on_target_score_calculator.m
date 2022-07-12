function final_score=on_target_score_calculator(s)
s_20mer=s(5:24);

%Create nuc_hash key
nuc_hash_key={'A' 'T' 'C' 'G'};
nuc_hash_val={0 1 2 3};
nuc_hash=containers.Map(nuc_hash_key,nuc_hash_val);

%Define Score
score=.597636154;
gc=sum(s_20mer=='G')+sum(s_20mer=='C');
gclow=-.202625894;
gchigh=-.166587752;
if gc<10
    gc_val=abs(gc-10);
    score=score+(gc_val*gclow);
elseif gc>10
    gc_val=gc-10;
    score=score+(gc_val*gchigh);
end

%Create sing_nuc_hash
sing_nuc_key={'G2' 'A3' 'C3' 'C4' 'C5' 'G5' 'A6' 'C6' 'C7' 'G7' 'A12'...
    'A15' 'C15' 'A16' 'C16' 'T16' 'A17' 'G17' 'C18' 'G18' 'A19' 'C19'...
    'G20' 'T20' 'G21' 'T21' 'C22' 'T22' 'T23' 'C24' 'G24' 'T24' 'A25'...
     'C25' 'T25' 'G28' 'T28' 'C29' 'G30'};
sing_nuc_val={-0.275377128 -0.323887456 0.172128871 -0.100666209...
    -0.20180294 0.245956633 0.036440041 0.098376835 -0.741181291...
    -0.393264397 -0.466099015 0.085376945 -0.013813972 0.272620512...
    -0.119022648 -0.285944222 0.097454592 -0.17554617 -0.345795451...
    -0.678096426 0.22508903 -0.507794051 -0.417373597 -0.054306959...
    0.379899366 -0.090712644 0.057823319 -0.530567296 -0.877007428...
    -0.876235846 0.278916259 -0.403102218 -0.077300704 0.287935617...
    -0.221637217 -0.689016682 0.117877577 -0.160445304 0.386342585};
sing_nuc_hash=containers.Map(sing_nuc_key,sing_nuc_val);

%Create dinuc_hash
dinuc_key={'GT2' 'GC5' 'AA6' 'TA6' 'GG7' 'GG12' 'TA12' 'TC12' 'TT12'...
    'GG13' 'GA14' 'GC14' 'TG17' 'GG19' 'TC19' 'CC20' 'TG20' 'AC21'...
    'CG21' 'GA21' 'GG21' 'TC22' 'CG23' 'CT23' 'AA24' 'AG24' 'AG25'...
    'CG25' 'TG25' 'GT27' 'GG29'};
dinuc_val={-0.625778696 0.300043317 -0.834836245 0.760627772 -0.490816749...
    -1.516907439 0.7092612 0.496298609 -0.586873894 -0.334563735 ...
    0.76384993 -0.53702517 -0.798146133 -0.66680873 0.353183252...
    0.748072092 -0.367266772 0.568209132 0.329072074 -0.836456755...
    -0.782207584 -1.029692957 0.856197823 -0.463207679 -0.579492389...
    0.649075537 -0.077300704 0.287935617 -0.221637217 0.117877577...
    -0.697740024};
dinuc_hash=containers.Map(dinuc_key,dinuc_val);
    
for i=1:length(s)
    key=strcat(s(i),num2str(i));
    if ismember(key,sing_nuc_key)==1
        nuc_score=sing_nuc_hash(key);
        score=score+nuc_score;
    else
        nuc_score=0;
    end
    if i<30
        dinuc=strcat(s(i),s(i+1),num2str(i));
        if ismember(dinuc,dinuc_key)==1
            dinuc_score=dinuc_hash(dinuc);
            score=dinuc_score+score;
        else
            dinuc_score=0;
        end
    end
end

    partial_score=exp(-score);
    final_score=1/(1+partial_score);
end