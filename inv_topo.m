im=sqrt(-1);
pol  = '_VV';
ddir = ['dates' pol '/'];
decide_ints_stack

for i=1:nd
    dates(i).name    = files(i).name(1:8);
    dates(i).dn      = datenum(dates(i).name,'YYYYmmdd');
end
dn=[dates.dn];
dn=dn-min(dn);

clear ints
for i=1:nd-1
    j=i+1;
    ints(i).dir   = (['intdir' pol '/' dates(i).name '/']);
    ints(i).name  = [dates(i).name '_' dates(j).name '_' num2str(rlooks) 'rlk_' num2str(alooks) 'alk'];
    ints(i).flat   = [ints(i).dir ints(i).name '_highpass.unw'];
    ints(i).msk   = [ints(i).dir ints(i).name '.msk'];
    ints(i).fix   = [ints(i).dir ints(i).name '_highpass_fix.unw'];
    fidi(i)=fopen(ints(i).flat,'r');
    fid(i)=fopen(ints(i).fix,'w');
end

fido1=fopen('demerr.r4','w');
fido2=fopen('sigstd.r4','w'); %after fit

[bp,intbp]=read_baselines;
intdn=diff(dn);
G=[intdn' intbp'];
Gg=inv(G'*G)*G';

data=nan(nd-1,newnx);
for j=1:newny
    for i=1:nd-1
        data(i,:)=fread(fidi(i),newnx,'real*4');
    end
    mods=Gg*data;
    slope=mods(2,:);
    synth=intbp'*slope;

    a=sqrt(mean(data.^2,1));
    b=sqrt(mean((data-G*mods).^2,1));
    %c=1-b./a; R^2

    fwrite(fido1,slope,'real*4');
    fwrite(fido2,b,'real*4');
    for i=1:nd-1
        fwrite(fid(i),data(i,:)-synth(i,:),'real*4');
    end
end

  