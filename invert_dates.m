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

maskfile='wgtdir/average.1alks_3rlks.amp';
maskthresh=0.14;


tsdir='tsdir/';
if(~exist(tsdir,'dir'))
    mkdir(tsdir);
end
ddir=['tsdir/fixunw/'];
if(~exist(ddir,'dir'))
    mkdir(ddir)
end
fdir=['tsdir/filt/'];
if(~exist(fdir,'dir'))
    mkdir(fdir)
end

fid=fopen(maskfile,'r');
mask=fread(fid,[newnx newny],'real*4');
mask=mask>maskthresh;
fclose(fid);
geomfile='merged/geom_master/hgt.rdr.1alks_3rlks.full';
fid=fopen(geomfile,'r');
for i=1:newny
    a=fread(fid,newnx,'real*8');
    mask(a<0,i)=0;
end
fclose(fid);


for i=1:nd-1
    j=i+1;
    intdir=(['intdir' pol '/' dates(i).name '/']);
    ints(i).name=[dates(i).name '_' dates(j).name '_' num2str(rlooks) 'rlk_' num2str(alooks) 'alk'];
    ints(i).unw=[intdir ints(i).name '.unw'];  
    ints(i).filt=[fdir ints(i).name '.unw'];
    ints(i).fixunw=[ddir ints(i).name '.unw'];
end

for i=1:nd
    dates(i).unw=[ddir dates(i).name '.unw'];
    dates(i).filt=[fdir dates(i).name '.unw'];
end


ry=500;
windx=ones(1,ry*2+1);
windy=ones(1,ry*2+1);
windx=windx/sum(windx);
windy=windy/sum(windy);
for i=1:nd-1
    i
    if(~exist(ints(i).filt,'file'))
        
        z       = zeros(1,newnx);     
        fidi    = fopen(ints(i).unw,'r');
        fido    = fopen(ints(i).filt,'w');
        
        mskunw=fread(fidi,[newnx,ry],'real*4');
        mskunw(~mask(:,1:ry))=0; %mask on fly
        mskunw = [flipud(mskunw');nan(ry+1,newnx)];
        mskunw(isnan(mskunw))=0;
        
        %now go to end (passing end by ry lines, filling in with zeros. "active"
        %line is at ry+1th row
        
        for j=1:newny
            
            mskunw=circshift(mskunw,1);
            [a,count1]=fread(fidi,newnx,'real*4');
            a(isnan(a))=0;
            if(count1==newnx)
                a(~mask(:,ry+j))=0;
                mskunw(1,:)=a;
            else
                mskunw(1,:)=z;
            end
            good = mskunw~=0;
            a    = good;
            c    = mskunw;
            
            a    = windy*a;
            c    = windy*c;
            
            asum = conv(a,windx,'same');
            csum = conv(c,windx,'same');
            out  = csum./asum;
            
            fwrite(fido,out,'real*4'); %1000pixel filtered product, at all pixels, even masked ones.
        end
        fclose(fidi);
        fclose(fido);
        fidi1   = fopen(ints(i).unw,'r');
        fidi2   = fopen(ints(i).filt,'r');
        fido    = fopen(ints(i).fixunw,'w');
        for j=1:newny
            a=fread(fidi1,newnx,'real*4');
            b=fread(fidi2,newnx,'real*4');
            fwrite(fido,a-b,'real*4');
        end
    end 
end


return


%make simple time series
G = -eye(nd);
G = G-circshift(G,[0,1]);
G = G(1:end-1,:);
Gr=G;
Gr(end+1,:)=1;
Gg = inv(Gr'*Gr)*G';

for i=1:nd
    fido(i)=fopen(dates(i).unw,'w');
end
for i=1:nd-1
    fidi(i)=fopen(ints(i).unw,'r');
end

for j=1:newny
    tmp=zeros(nd-1,newnx);
    for i=1:nd-1
        tmp(i,:)=fread(fidi(i),newnx,'real*4');
    end
    model=Gg*tmp;
    
    for i=1:nd
        fwrite(fido(i),model(i,:),'real*4');
    end
end
fclose('all');
