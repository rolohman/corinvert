im=sqrt(-1);
pol='_VV';
ddir=['dates' pol '/'];
decide_ints_stack

for i=1:nd
    dates(i).name    = files(i).name(1:8);
    dates(i).dn      = datenum(dates(i).name,'YYYYmmdd');
end
dn=[dates.dn];
dn=dn-min(dn);

if(~exist(ddir,'dir'))
    mkdir(ddir)
end

for i=1:nd-1
    j=i+1;
    intdir         = (['intdir' pol '/' dates(i).name '/']);
    ints(i).name   = [dates(i).name '_' dates(j).name '_' num2str(rlooks) 'rlk_' num2str(alooks) 'alk'];
    ints(i).unw    = [intdir ints(i).name '_highpass.unw'];
    ints(i).fix    = [intdir ints(i).name '_highpass_fix.unw'];
    ints(i).mask   = [intdir ints(i).name '.msk'];
end

for i=1:nd
    dates(i).unw     = [ddir dates(i).name '_simple.unw'];
    dates(i).infill  = [ddir dates(i).name '_infill.unw'];
end

dti         = diff(dn);
[bp,intbp]  = read_baselines;

%make simple time seriestest
Gi          = -eye(nd);
Gi          = Gi-circshift(Gi,[0,1]);
Gi          = Gi(1:end-1,:);

G           = [Gi]; %datanoise intercept rate demerr

Gr          = G;
Gr(end+1,1) = 1; %datanoise at pt 1=0 (will shift next

Gr=Gr*[ones(nd,1) dn' bp'];

Gg          = inv(Gr'*Gr)*Gr';



if(~exist(dates(1).unw,'file'))
    for i=1:nd
        fido(i)  = fopen(dates(i).unw,'w');
    end
    for i=1:nd-1
        fidi(i)  = fopen(ints(i).unw,'r');
        fido2(i) = fopen(ints(i).fix,'w');
    end
    fida=fopen([ddir 'simplevel.r4'],'w');
    fidb=fopen([ddir 'simpleres.r4'],'w');
    fidc=fopen([ddir 'simpledemerr.r4'],'w');
    fidd=fopen([ddir 'sigR2.r4'],'w'); %after fit
    
    for j=1:newny
        tmp=zeros(nd-1,newnx);
        for i=1:nd-1
            tmp(i,:)=fread(fidi(i),newnx,'real*4');
        end
        model     = Gg*[tmp;zeros(1,newnx)];
        intercept = model(1,:);
        vel       = model(2,:);
        demerr    = model(3,:);
  
        
        synth     = Gr*model;
        res       = tmp-synth(1:end-1,:);
        resd      = [zeros(1,newnx);cumsum(res)]-intercept;
        for i=1:nd
            fwrite(fido(i),resd(i,:),'real*4');
        end
        for i=1:nd-1
            fwrite(fido2(i),res(i,:),'real*4');
        end
        resv  = sqrt(mean(res.^2,1,'omitnan'));
        orig  = sqrt(mean(tmp.^2,1,'omitnan'));
        R2    = 1-resv./orig;
        fwrite(fida,vel,'real*4');
        fwrite(fidb,resv,'real*4');
        fwrite(fidc,demerr,'real*4');
        fwrite(fidd,R2,'real*4');
    end
    fclose('all');
end
