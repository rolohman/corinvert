function cr=est_cr(cres,synth,nd,cidl)

crs              = zeros(nd);
crs(cidl)        = cres;
crs              = crs+crs';
crs(crs==0)      = NaN;

   
if(isempty(synth))
    cerr=0.3083*(1-crs)-.2*(1-crs).^2;
    wgts=1./cerr;
else
    wgts             = zeros(nd);
    wgts(cidl)       = synth;  %0 to 1
    wgts             = wgts+wgts';
    wgts(wgts==0)    = NaN;
    wgts             = exp((wgts-1)*5); %scalar makes near 0.1 at synth=0.5;
end
cr=mymed(crs,wgts);
% crsnew=-(-abs(repmat(cr,nd,1)-repmat(cr',1,nd))-repmat(cr',1,nd)'-crs);
% cr=mymed(crsnew,wgts);


%
% return
% figure
% plot(cr)
% crs2=crs-repmat(cr',1,nd);
% for i=1:nd
%     cr(i)=weightedMedian(crs2(:,i),wgts(:,i));
% end
% hold on
% plot(cr)
% 
% % wgtsum1          = sum(crs.*wgts,1,'omitnan');
% % wgtsum2          = sum(wgts,1,'omitnan');
% % cr               = wgtsum1./wgtsum2;
% % 
% % %iterate, downweighting strong signals
% % crs2             = crs-repmat(cr',1,nd);
% % wgtsum1          = sum(crs2.*wgts,1,'omitnan');
% % cr               = wgtsum1./wgtsum2;
% % 
% % cr               = cr-mymax(cr,10);
% % cr(isnan(cr))    = 0;
% 
% if(~islog)
%     cr=cr+1;
% end