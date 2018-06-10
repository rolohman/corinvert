decide_ints
nd    = length(dates);

rlooks = 7;
alooks = 3;
newnx  = floor(nx/rlooks)
newny  = floor(ny/alooks);


%open time inversion files

for i=1:nd
    fidr(i)=fopen(['results7/rel_' dates(i).name '.cor'],'r');
end
fids    = fopen('results_TS/shift','r');
for i=1:nd
    fidro(i)=fopen(['results7/rel_' dates(i).name '_shifted.cor'],'w');
end

%load data and models
for i=1:newny    
    shift=fread(fids,newnx,'real*4');
    for i=1:nd  
        rels=fread(fidr(i),newnx,'real*4');
        fwrite(fidro(i),rels-shift,'real*4');
    end
end
fclose('all');


    