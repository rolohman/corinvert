pol  = '_VV';
ddir = ['dates' pol '/'];
decide_ints_stack
for i=1:nd
    dates(i).name    = files(i).name(1:8);
    dates(i).dn      = datenum(dates(i).name,'YYYYmmdd');
    dates(i).unw     = [ddir dates(i).name '_simple.unw'];
    dates(i).fixunw  = [ddir dates(i).name '_simple_fix.unw'];
end
dn=[dates.dn];


for i=1:nd
    fid(i)=fopen(dates(i).unw,'r');
    fido(i)=fopen(dates(i).fixunw,'w');
end
def=zeros(nd,newnx);
for j=1:newny
    for i=1:nd
        def(i,:)=fread(fid(i),newnx,'real*4');
    end
    for i=1:nd
        win=max(i-10,1):i;
        smoo=median(def(win,:),1,'omitnan');
        cycles=sign(fix((def(i,:)-smoo)/pi)); %pos or neg if abs(dif)>pi
        def(i,:)=def(i,:)+2*pi*cycles;
    end
    for i=1:nd
        fwrite(fido(i),def(i,:),'real*4');
    end
    
end
    