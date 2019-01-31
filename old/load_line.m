function cors = load_line(line)

decide_ints
for i=1:ni
    cordir=(['cordir/' dates(id1(i)).name '/']);
    intcordir=[cordir dates(id2(i)).name '/'];
    corfile_small=[intcordir dates(id1(i)).name '_' dates(id2(i)).name '_' num2str(rlooks) 'rlk_' num2str(alooks) 'alk.cor'];
    if(exist(corfile_small,'file'))
        fidi(i)=fopen(corfile_small,'r');
    else
        disp(['where is ' corfile_small '?'])
    end
end

for i=1:ni
    fseek(fidi(i),newnx*(line-1)*4,-1);
end

cors=zeros(ni,newnx);

for i=1:ni
    [tmp,count]=fread(fidi(i),newnx,'real*4');
    if(count>0)
        cors(i,1:count)=tmp;
    end
end
    