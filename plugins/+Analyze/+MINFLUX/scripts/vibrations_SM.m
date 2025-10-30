% calculate average FFT from FoV and other measures
% we expect constant dt, so low photon threshold
% parameters
minlen=5000; %minimum number of localizations in track
maxfreq=350;
skipfirst=5;
dfreq=0.1;

usefields={'xnm','ynm','znm','time','groupindex','numberInGroup','filenumber','efo','cfr','eco','ecc','efc','tid','fbg','phot'};
locs=g.locData.getloc(usefields,'layer',find(g.getPar('sr_layerson')),'Position','all','removeFilter',{'filenumber','time'});
mtid=max(double(locs.tid));
tidf=double(locs.tid)+double(locs.filenumber)*mtid;
indlong=locs.numberInGroup>=minlen;
trackids=unique(tidf(indlong));
figure(88); clf
ax1=subplot(2,1,1);
ax2=subplot(2,1,2);
% dfreq=1/min(diff(locs.time(indlong)))
fall=0:dfreq:maxfreq;
ftxa=0;ftya=0;


for k=1:length(trackids)
    ig=tidf==trackids(k);
    igf=find(ig);
    [yfp,freq]=getfft(locs.ynm(igf(skipfirst:end)),locs.time(igf(skipfirst:end)));
    
    yfpbh=bindata(freq,yfp,fall);
    plot(ax2,fall,yfpbh,'c'); hold(ax1,'on') 
    yfpbh(isnan(yfpbh))=0;
    ftya=yfpbh*sum(ig)+ftya;
    % end

    [xfp,freq]=getfft(locs.xnm(igf(skipfirst:end)),locs.time(igf(skipfirst:end)));
    xfpbh=bindata(freq,xfp,fall);
    plot(ax1,fall,xfpbh,'c'); hold(ax2,'on') 
    xfpbh(isnan(xfpbh))=0;
    ftxa=xfpbh*sum(ig)+ftxa; 
    % end
end
% ftxa=ftxa/length(trackids);
ftxa=ftxa/length(ig);
ftya=ftya/length(ig);
plot(ax1,fall,ftxa,'k','LineWidth',2);
xlabel(ax1,'frequency (Hz)')
ylabel(ax1,'Amplitude (nm), 2*fft (x)');
axis(ax1,'tight')
ax1.YLim(2)=0.25;

plot(ax2,fall,ftya,'k','LineWidth',2);
xlabel(ax2,'frequency (Hz)')
ylabel(ax2,'Amplitude (nm), 2*fft (x)');
axis(ax2,'tight')
ax2.YLim(2)=0.25;

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
