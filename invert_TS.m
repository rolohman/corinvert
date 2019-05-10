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
    ints(i).infill = [intdir ints(i).name '_highpass_infill.unw'];
    ints(i).mask   = [intdir ints(i).name '.msk'];
end

for i=1:nd
    dates(i).unw     = [ddir dates(i).name '_simple.unw'];
    dates(i).infill  = [ddir dates(i).name '_infill.unw'];
end


%make simple time series
G           = -eye(nd);
G           = G-circshift(G,[0,1]);
G           = G(1:end-1,:);
Gr          = G;
Gr(end+1,:) = 1;
Gg          = inv(Gr'*Gr)*G';

%slope/int
dn          = [dates.dn]';
dn          = dn-min(dn);
G2          = [ones(nd,1) dn];
Gg2         = inv(G2'*G2)*G2';

if(~exist(dates(1).unw,'file'))
    for i=1:nd
        fido(i)=fopen(dates(i).unw,'w');
    end
    for i=1:nd-1
        fidi(i)=fopen(ints(i).unw,'r');
    end
    fida=fopen([ddir 'simpleslope.r4'],'w');
    fidb=fopen([ddir 'simpleint.r4'],'w');
    fidc=fopen([ddir 'simpleres.r4'],'w');
    
    for j=1:newny
        tmp=zeros(nd-1,newnx);
        for i=1:nd-1
            tmp(i,:)=fread(fidi(i),newnx,'real*4');
        end
        model=Gg*tmp;
        synth=G*model;
        for i=1:nd
            fwrite(fido(i),model(i,:),'real*4');
        end
        mod2=Gg2*model;
        synth=G2*mod2;
        res2=model-synth;
        res2=sqrt(mean(res2.^2,1,'omitnan'));
        fwrite(fida,mod2(2,:),'real*4');
        fwrite(fidb,mod2(1,:),'real*4');
        fwrite(fidc,res2,'real*4');
    end
    fclose('all');
end
%option 2 - filter the masked parts, then simple time series

for i=1:nd-1
    j=i+1;
    if(~exist(ints(i).infill,'file'))
        disp(['making ' ints(i).infill])
        myfilt(ints(i).unw,ints(i).mask,'tmpfilt',150,150,newnx,newny,2,3,4,'/dev/null');
        %replace masked regions with filtered (in matlab since there is
        %some issue with integer*1 for imagemath
        fidi=fopen(ints(i).unw,'r');
        fidm=fopen(ints(i).mask,'r');
        fidf=fopen('tmpfilt','r');
        fido=fopen(ints(i).infill,'w');
        for k=1:newny
            a=fread(fidi,newnx,'real*4');
            b=fread(fidm,newnx,'integer*1');
            c=fread(fidf,newnx,'real*4');
            mask=and(b==0,isfinite(a));
            a(mask)=c(mask);
            fwrite(fido,a,'real*4');
        end
    else
        disp([ints(i).infill ' already made'])
    end
end

if(~exist(dates(1).infill,'file'))
    for i=1:nd
        fido(i)=fopen(dates(i).infill,'w');
    end
    for i=1:nd-1
        fidi(i)=fopen(ints(i).infill,'r');
    end
    fida=fopen([ddir 'infillslope.r4'],'w');
    fidb=fopen([ddir 'infillint.r4'],'w');
    fidc=fopen([ddir 'infillres.r4'],'w');
    
    for j=1:newny
        tmp=zeros(nd-1,newnx);
        for i=1:nd-1
            tmp(i,:)=fread(fidi(i),newnx,'real*4');
        end
        model=Gg*tmp;
        for i=1:nd
            fwrite(fido(i),model(i,:),'real*4');
        end
        mod2=Gg2*model;
        synth=G2*mod2;
        res2=model-synth;
        res2=sqrt(mean(res2.^2,1,'omitnan'));
        fwrite(fida,mod2(2,:),'real*4');
        fwrite(fidb,mod2(1,:),'real*4');
        fwrite(fidc,res2,'real*4');
    end
    fclose('all');
    
end
