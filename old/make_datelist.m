function make_datelist(pth)
%pth=28 for Saudi
%latlon.txt must exist, with one line of text, i.e., [18.2,21.5,52.5,55]


%This also requires that you have downloaded jsonlab (has loadjson)

fid=fopen('latlon.txt','r');
a=fgetl(fid);
fclose(fid);
a=a(2:end-1); %get rid of brackets
b=strsplit(a,',');
minlat=str2num(b{1});
maxlat=str2num(b{2});
minlon=str2num(b{3});
maxlon=str2num(b{4});

jnk=[minlon minlat maxlon minlat maxlon maxlat minlon maxlat minlon minlat];
polygon='polygon=';
for i=1:9
    polygon=[polygon num2str(jnk(i)) ','];
end
polygon=[polygon num2str(jnk(end))];

rootadr='https://api.daac.asf.alaska.edu/services/search/param?';
apicall=[rootadr polygon '&platform=Sentinel-1a,Sentinel-1b&processingLevel=SLC&start=Jan 1, 2017&relativeOrbit=' num2str(pth) '&output=json'];
c=webread(apicall);

%now pull the whole thing back into a structure
d=loadjson(c);
d=d{1};
e(1)=d{1};
for i=1:length(d)
    e(i)=d{i};
    dn(i)=datenum(e(i).sceneDate(1:10));
end

dn=unique(dn);
ds=datestr(dn,'yyyymmdd');
fid=fopen('alldates.txt','w');
for i=1:length(ds)
    fprintf(fid,'%s\n',ds(i,:));
end
fclose(fid);

%edit this file to remove or add dates if necessary


