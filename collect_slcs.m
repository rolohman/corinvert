decide_ints
rx=7;
ry=3;

fid=fopen('commands','w');
d1=datestr(dn(testid),'yyyymmdd');
datedir=[d1 '/'];


d2=datestr(dn(2),'yyyymmdd');
intname=['int_' d1 '_' d2];
intdir=[datedir intname '/'];
primary    = [intdir 'merged/master.slc.full'];
if(~exist(primary))
    command = ['looks.py -i ' primary ' -o ' primary ' -r 1 -a 1\n'];
    fprintf(fid,command);
end
dates(testid).slc=primary;


for i=[1:testid-1 testid+1:nd]
    d2         = datestr(dn(i),'yyyymmdd');
    intname    = ['int_' d1 '_' d2];
    intdir     = [datedir intname '/'];
    secondary  = [intdir 'merged/slave.slc.full'];
    dates(i).slc=secondary;
    dates(i).flatint=[intdir 'merged/topophase.flat.full'];
    dates(i).fullint=[intdir 'merged/fullint.full'];
    dates(i).ramp=[intdir 'merged/ramp.full'];
    dates(i).flatslc=[intdir 'merged/' d2 '.slc.full'];
    
    if(~exist(secondary))
        command = ['looks.py -i ' secondary ' -o ' secondary ' -r 1 -a 1\n'];
        fprintf(fid,command);
    else
        disp([secondary ' already made']);
    end
    if(~exist(dates(i).flatint))
        command=['looks.py -i ' dates(i).flatint ' -o ' dates(i).flatint ' -r 1 -a 1\n'];
        fprintf(fid,command);
    else
        disp([dates(i).flatint ' already made'])
    end
    if(~exist(dates(i).ramp))
        command=['imageMath.py -e=''a*conj(b)'' -o tmpint -t cfloat --a=' dates(testid).slc ' --b=' dates(i).slc '\n'];
        fprintf(fid,command);
        command=['imageMath.py -e=''a*b'' -o ' dates(i).ramp ' -t cfloat --a=tmpint --b=' dates(i).flatint '\n'];
        fprintf(fid,command);
    else
        disp([dates(i).ramp ' already made'])
    end
    if(~exist(dates(i).flatslc))
        command=['imageMath.py -e=''a*b'' -o ' dates(i).flatslc ' -t cfloat --a=' dates(i).slc ' --b=' dates(i).ramp '\n'];
        fprintf(fid,command);
    else
        disp([dates(i).flatslc ' already made'])
    end
    
end


fclose('all')
return

for i=2:nd
    d2         = datestr(dn(i),'yyyymmdd');
    intname    = ['int_' d1 '_' d2];
    intdir     = [datedir intname '/'];
    
end
fclose(fid);


return


for i=1:length(ints2make)
    d1=datestr(ints2make(i,1),'yyyymmdd');
    d2=datestr(ints2make(i,2),'yyyymmdd');
    intname = ['int_' d1 '_' d2];
    datedir = [d1 '/'];
    intdir  = [datedir intname '/'];
    
    infile=([intdir 'merged/topophase.flat.full']);
    outfile=([intdir 'merged/mycor.cor.full']);
    
    if(exist(infile,'file'))
        if(~exist(outfile,'file'))
            a=xml2struct([infile '.xml']);
            n=a.imageFile.property;
            for j=1:length(n)
                tmp=n{j}.Attributes.name;
                if(regexp(tmp,'length'))
                    ny=str2num(n{j}.value.Text);
                end
                if(regexp(tmp,'width'))
                    nx=str2num(n{j}.value.Text);
                end
            end
            
        end
    end
end

