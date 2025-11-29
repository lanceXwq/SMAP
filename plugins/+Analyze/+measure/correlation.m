classdef correlation<interfaces.DialogProcessor
   % also look at FFT
   % empty / add button to add up FFT, AC
   % include 2D analysis? put 1D to ROI?
   % PERPL: calculate 2D g(r)? fit with lines or grid (2 distances)?
    methods
        function obj=correlation(varargin)        
            obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={'sr_layerson','linewidth_roi','znm_min','znm_max','sr_pixrec','layernames'};
            obj.showresults=true;
        end
        
        function out=run(obj,p)    
            if p.linecorrelation
                axcorr=obj.initaxis('corr1D');
                hold(axcorr,'off')
                axprof=obj.initaxis('prof1D');
                hold(axprof,'off')             
                axfft=obj.initaxis('FFT');
                hold(axfft,'off')     
            end
            layers=find(p.sr_layerson);
            k=1;
            [locs,~, hroi]=obj.locData.getloc({'xnm','ynm','znm','locprecnm','locprecznm','xnmline','ynmline'},'layer',layers(k),'position','roi');
            linew=p.linewidth_roi/2;
            if isa(hroi,'imline')
                x=locs.xnmline;
                y=locs.ynmline;
                mx=0;
                my=0;
                pos=getPosition(hroi);
                linel=sqrt((pos(2,1)-pos(1,1))^2+(pos(2,2)-pos(1,2))^2)/2*1000;
            else
                x=locs.xnm;
                y=locs.ynm;
                pos=getPosition(hroi);
                mx=mean(pos(:,1));
                my=mean(pos(:,2));
            end
            if ~p.setbinwidth
                p.binwidth=p.sr_pixrec;
            end
            if p.linecorrelation
                [prof,n]=lineprof(x-mx,y-my,p,axprof);
                linecorr(prof,n,p,axcorr)
                fftprof(prof,n,p,axfft)
                
            end


        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end

function ach=linecorr(prof,x,p,axcorr)
% n=min(x):p.binwidth:max(x);
% nplot=n(1:end-1)+p.binwidth/2;
% prof = histcounts(x,  n);
ac=xcorr(prof,prof);
ach=ac(ceil(end/2):end)/max(ac);
nac=(0:length(ach)-1)'*p.binwidth;
[yRPeaks,xPeaksIdx]=findpeaks(detrend(ach,'linear'));
xPeaksIdx(xPeaksIdx<3)=[]; yRPeaks(xPeaksIdx<3)=[]; 

[~,imxp]=min(abs(nac(xPeaksIdx)-p.period));
period=nac(xPeaksIdx(imxp));
% [yRPeaks,xRPeaks] = refinepeaks(ach,xPeaksIdx,nac);
h1=plot(axcorr,nac,ach,'DisplayName',"correlation, P="+num2str(period));


acs=ac/max(ac);
% figure(88); hold off
for k=1:3
    acs=acs(ceil(end/2):end);
    acs=detrend(acs,'linear');
    acs=xcorr(acs-mean(acs),acs-mean(acs));
    acs=acs/max(acs);
    % plot(acs)
    % hold on

end


acmulticorr=0.5*(1+acs(ceil(end/2):ceil(end/2)+length(nac)-1));
acmulticorr=acmulticorr-min(acmulticorr);
acmulticorr=acmulticorr/mean(acmulticorr)*mean(ach);

hold(axcorr,'on');
[yRPeaks,xPeaksIdx]=findpeaks(acmulticorr);
[yRPeaks,xRPeaks] = refinepeaks(acmulticorr,xPeaksIdx,nac);

[~,imxp]=min(abs(xRPeaks-p.period));
% period=mean(diff(ppos));
period=xRPeaks(imxp);
ff="%2.1f";

h2=plot(axcorr,nac,acmulticorr,'--','DisplayName',"multi-c: P="+num2str(period,ff));
legend(axcorr)
end

function [prof,n]=lineprof(x,y,p,axprof)
n=min(x):p.binwidth:max(x);
nplot=n(1:end-1)+p.binwidth/2;
prof = histcounts(x,  n);
proffilt=smoothdata(prof,'gaussian',12);
plot(axprof,nplot,prof,nplot,proffilt)
[~,xpos]=findpeaks(proffilt,nplot);
dx=diff(xpos);
tpos=(xpos(1:end-1)+xpos(2:end))/2;
text(axprof,tpos,max(proffilt)+0*tpos,string(dx))
end

function fftprof(prof,nx,p,axfft)
[f, mag1] = fft_one_sided(nx,[prof zeros(1,0)]);


[yRPeaks,xPeaksIdx]=findpeaks(mag1);
[yRPeaks,xRPeaks] = refinepeaks(mag1,xPeaksIdx,f);
periods=1./xRPeaks

[~,imxp]=min(abs(periods-p.period));
period=periods(imxp);
plot(axfft, f, mag1,'DisplayName',"P: "+num2str(period,"%2.1f"))%,
hold(axfft,'on')
plot(axfft,xRPeaks(imxp),yRPeaks(imxp),'kx')
legend(axfft)
end

function [f, mag1] = fft_one_sided(t, x)
% Returns one-sided frequency axis and magnitude |X(f)| for real x(t)
t = t(:); x = x(:);
dt = diff(t);
if max(abs(dt - mean(dt))) > 1e-6*mean(dt)
    error('Time vector is not uniformly sampled.');
end
Fs = 1/mean(dt);
N = numel(x);
X = fft(x);
mag1 = abs(X(1:floor(N/2)+1))/N;
if N > 2
    mag1(2:end-1) = 2*mag1(2:end-1);
end
f = (0:floor(N/2))*(Fs/N);
end

function pard=guidef(obj)
pard.inputParameters={'sr_layerson','linewidth_roi','sr_pixrec'};

 p(1).value=0; p(1).on={}; p(1).off={'binwidth'};
p(2).value=1; p(2).on={'binwidth'}; p(2).off={};
pard.setbinwidth.object=struct('String','set binwidth (nm) (otherwise: pixelsize):','Style','checkbox','Value',1,'Callback',{{@obj.switchvisible,p}});
pard.setbinwidth.position=[1,1];
pard.setbinwidth.Width=3;

pard.binwidth.object=struct('String','10','Style','edit');
pard.binwidth.position=[1,3.5];
pard.binwidth.Width=0.5;
pard.binwidth.TooltipString='Binwidth for profiles. If not checked, use pixel size of reconstruction';
pard.setbinwidth.TooltipString=pard.binwidth.TooltipString;

pard.linecorrelation.object=struct('String','1D correlation','Style','checkbox','Value',1);
pard.linecorrelation.position=[2,1];
pard.linecorrelation.Width=1;

pard.c2d.object=struct('String','2D correlation','Style','checkbox','Value',0);
pard.c2d.position=[2,2];
pard.c2d.Width=1;


pard.periodt.object=struct('String','approximate period (nm)','Style','text');
pard.periodt.position=[3,1];
pard.periodt.Width=1.5;
pard.period.object=struct('String','200','Style','edit');
pard.period.position=[3,2.5];
pard.period.Width=0.5;

% pard.text2.object=struct('String','fitmodel:','Style','text');
% pard.text2.position=[2,1];
% 
% pard.fitmodel.object=struct('String','Gauss|Flat|Disk|Ring|Distance','Style','popupmenu');
% pard.fitmodel.position=[2,2];
% 
% pard.restrictsigma.object=struct('String','sigma=<locp>','Style','checkbox');
% pard.restrictsigma.position=[2,3];
% 
% 
% p(1).value=0; p(1).on={}; p(1).off={'linelength'};
% p(2).value=1; p(2).on={'linelength'}; p(2).off={};
% pard.linelengthcheck.object=struct('String','set length (nm)','Style','checkbox','Callback',{{@obj.switchvisible,p}});
% pard.linelengthcheck.position=[3,1];
% 
% pard.linelength.object=struct('String','250','Style','edit');
% pard.linelength.position=[3,2];
% pard.linelength.TooltipString='This overrides the length of the ROI and uses a well-defined ROI. Useful for direct comparison.';
% pard.linelengthcheck.TooltipString=pard.linelength.TooltipString;
pard.plugininfo.name='Correlation';
pard.plugininfo.description=sprintf('correlation analysis');
pard.plugininfo.type='ProcessorPlugin';
end