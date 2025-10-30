% calculate average FFT from FoV and other measures
% we expect constant dt, so low photon threshold
% parameters
p.minlen=5000; %minimum number of localizations in track
p.maxfreq=350;
p.skipfirst=5;
p.dfreq=0.1;
p.ymax=0.25;

figure(88); clf
ax1=subplot(2,1,1);
plotfft("xnm",p,g,ax1);
ax2=subplot(2,1,2);
plotfft("ynm",p,g,ax2);

function [xfp,freq,peaks]=plotfft(field,p,g,ax)
if nargin<4
    ax=gca;
end
usefields={'xnm','ynm','znm','time','groupindex','numberInGroup','filenumber','efo','cfr','eco','ecc','efc','tid','fbg','phot'};
locs=g.locData.getloc(usefields,'layer',find(g.getPar('sr_layerson')),'Position','all','removeFilter',{'filenumber','time'},'grouping','ungrouped');
mtid=max(double(locs.tid));
tidf=double(locs.tid)+double(locs.filenumber)*mtid;
indlong=locs.numberInGroup>=p.minlen;
trackids=unique(tidf(indlong));
fall=0:p.dfreq:p.maxfreq;
ftxa=0;

x=locs.(field);t=locs.time;
for k=1:length(trackids)
    ig=tidf==trackids(k);
    igf=find(ig);
    [xfp,freq]=getfft(x(igf(p.skipfirst:end)),t(igf(p.skipfirst:end)));
    xfpbh=bindata(freq,xfp,fall);
    plot(ax,fall,xfpbh,'c','LineWidth',0.1); hold(ax,'on') 
    xfpbh(isnan(xfpbh))=0;
    ftxa=xfpbh*sum(ig)+ftxa; 
end
% ftxa=ftxa/length(trackids);
ftxa=ftxa/length(ig);
plot(ax,fall,ftxa,'k','LineWidth',2);
xlabel(ax,'frequency (Hz)')
ylabel(ax,'Amplitude (nm), 2*fft (x)');
axis(ax,'tight')
ax.YLim(2)=p.ymax;

[pks,pkind]=findpeaks(ftxa,"MinPeakHeight",2*mean(ftxa),"MinPeakDistance",5);
% plot(fall(pkind),pks,'ro')
for k=1:length(pks)
    text(fall(pkind(k))+1,pks(k),string(fall(pkind(k))))
end
peaks=horzcat(fall(pkind)',pks');
end

function [xfp,freq]=getfft(x,t)
dt=mode(diff(t));
Fs=1e3/dt;
L=length(x);
xf=abs(fft(x)/L);
xfp=2*xf(1:floor(L/2)+1);
xfp(1)=0;
freq=(Fs*(0:(L/2))/L)';
% hold(axfx,'off')
% semilogy(axfx,freq,xfp,'r');


% if length(freq)>2500

    % hold(axfx,'on')
    % semilogy(axfx,f2,xfpb,'b');
    % ylim(axfx,[1e-5 1])
% end
% xlabel(axfx,'frequency (Hz)')
% ylabel(axfx,'Amplitude (nm), 2*fft (x)');
% axis(axfx,'tight')
% axfx.YLim(1)=quantile(xfp,0.01);
end
