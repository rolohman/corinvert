function [bp,intbp] = read_baselines
%reads all baselines
pol='_VV';
decide_ints_stack
basedir='baselines/';
for i=1:nd
    dates(i).name    = files(i).name(1:8);
    dates(i).dn      = datenum(dates(i).name,'YYYYmmdd');
    if(i>1)
        name=[dates(1).name '_' dates(i).name];
        dates(i).basefile=[basedir name '/' name '.txt'];
        tmp=dlmread(dates(i).basefile,':',1,1);
        bp(i)=tmp(1);
    end
end


intbp=diff(bp);