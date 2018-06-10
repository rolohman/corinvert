
if(exist('dnr','var'))
for i=1:length(dnr)
    plot([dnr(i) max(dn)],[dnr(i) dnr(i)],'m')
    plot([dnr(i) dnr(i)],[min(dn) dnr(i)],'m')
end
end
if(badi)
    for i=1:length(badi)
        plot([dni(badi(i)) dni(badi(i))],[min(dn) dn(badi(i)+1)],'-','color',[.5 .5 .5])
    end
end
datetick('x','mmmYY')
datetick('y','mmmYY')
axis tight
if(exist('ticks','var'))
set(gca,'xtick',ticks,'ytick',ticks);
set(gca,'XtickLabel',datestr(ticks,'mmmyy'))
set(gca,'YtickLabel',datestr(ticks,'mmmyy'))
end
caxis([0 1])


