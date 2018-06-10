decide_ints
rx=7;
ry=3;
fid1=fopen('lookfile','w');
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
            
            
            mycor4(infile,[],outfile,nx,ny,rx,ry,0,0);
        end
    else
        disp(['did not find ' infile])
        if(exist([infile '.vrt']))
            disp('trying to look down')
            fprintf(fid1,['looks.py -i ' infile ' -o tmp -r 1 -a 1\n']);
            fprintf(fid1,['mv tmp ' infile '\n']);
        else
            disp('no vrt file found either')
        end
    end
end
fclose(fid1);