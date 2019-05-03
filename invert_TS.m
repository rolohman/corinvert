pol='_VV';
decide_ints_stack
slcdir=['merged/SLC' pol '/'];
skips=1;  %1=sequential, larger=longer time pairs
for i=1:nd
    dates(i).name    = files(i).name(1:8);
    dates(i).slc     = [slcdir dates(i).name '/' dates(i).name '.slc.full'];
    dates(i).dn      = datenum(dates(i).name,'YYYYmmdd');
end
dn=[dates.dn];

ddir=['dates' pol '/'];
if(~exist(ddir,'dir'))
    mkdir(ddir)
end

for i=1:nd-1
    j=i+1;
    intdir         = (['intdir' pol '/' dates(i).name '/']);
    ints(i).name   = [dates(i).name '_' dates(j).name '_' num2str(rlooks) 'rlk_' num2str(alooks) 'alk'];
    ints(i).unw    = [intdir ints(i).name '_highpass.unw'];  

end

for i=1:nd
    dates(i).unw  = [ddir dates(i).name '_simple.unw'];
end


%make simple time series
G           = -eye(nd);
G           = G-circshift(G,[0,1]);
G           = G(1:end-1,:);
Gr          = G;
Gr(end+1,:) = 1;
Gg          = inv(Gr'*Gr)*G';

for i=1:nd
    fido(i)=fopen(dates(i).unw,'w');
end
for i=1:nd-1
    fidi(i)=fopen(ints(i).unw,'r');
end
fidr(i)=fopen([ddir 'simpleres.r4']);

for j=1:newny
    tmp=zeros(nd-1,newnx);
    for i=1:nd-1
        tmp(i,:)=fread(fidi(i),newnx,'real*4');
    end
    model=Gg*tmp;
    synth=G*model;
    res=tmp-synth;
    res=sqrt(mean(res.^2,1,'omitnan'));
    for i=1:nd
        fwrite(fido(i),model(i,:),'real*4');
    end
    fwrite(fidr,res,'real*4');
end
fclose('all');
