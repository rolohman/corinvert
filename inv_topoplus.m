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
   cordir=(['cordir2' pol '/' dates(i).name '/']);
    ints(i).dir   = (['intdir' pol '/' dates(i).name '/']);
    ints(i).name  = [dates(i).name '_' dates(j).name '_' num2str(rlooks) 'rlk_' num2str(alooks) 'alk'];
    ints(i).flat   = [ints(i).dir ints(i).name '_highpass.unw'];
    ints(i).msk   = [ints(i).dir ints(i).name '.msk'];
    ints(i).fix   = [ints(i).dir ints(i).name '_highpass_fixplus.unw'];
    ints(i).cor   = [cordir ints(i).name '.cor'];
    fidi(i)=fopen(ints(i).flat,'r');
    fidc(i)=fopen(ints(i).cor,'r');

    fid(i)=fopen(ints(i).fix,'w');
end

fido1=fopen('a.r4','w');
fido2=fopen('b.r4','w');
fido3=fopen('c.r4','w'); %after fit
fido4=fopen('d.r4','w');
fido5=fopen('sigstdplus.r4','w'); %after fit
fido6=fopen('R2.r4','w'); %after fit


[bp,intbp]=read_baselines;
intdn=diff(dn);
%G=[intdn' intbp'];
%Gg=inv(G'*G)*G';

intd=dn(1:end-1)+intdn/2;

for j=1:newny
    j
    data=nan(nd-1,newnx);
    cors=data;
    synth=data;
    
    for i=1:nd-1
        data(i,:)=fread(fidi(i),newnx,'real*4');
        cors(i,:)=fread(fidc(i),newnx,'real*4');
    end
    var     = -2*log(cors);
    weights = 1./var;
    good=and(isfinite(weights),isfinite(data));
    count=sum(good,1);
    
    a=nan(1,newnx);
    b=a;
    c=a;
    d=a;
    
    R2=a;
    rms=a;
    
    for i=find(count>100)
        ids=find(good(:,i));
        [fitresult, gof] = crazyfit(intd(ids), intbp(ids), data(ids,i)', weights(ids,i)');
        a(i)=fitresult.a;
        b(i)=fitresult.b;
        c(i)=fitresult.c;
        d(i)=fitresult.d;
        R2(i)=gof.rsquare;
        rms(i)=gof.rmse;
        synth(:,i)=fitresult(intd,intbp);
    end

 

    fwrite(fido1,a,'real*4');
    fwrite(fido2,b,'real*4');
    fwrite(fido3,c,'real*4');
    fwrite(fido4,d,'real*4');
    fwrite(fido5,rms,'real*4');
    fwrite(fido6,R2,'real*4');
    for i=1:nd-1
        fwrite(fid(i),data(i,:)-synth(i,:),'real*4');
    end
end

  