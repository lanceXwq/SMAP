classdef correlation<interfaces.DialogProcessor
    % LINEPROFILE Calculates profiles along a linear ROI and fits it with a
    % model of choice. Flat: step function convolved with Gaussian
    % (=Erf). Disk: Projection of a homogeneously filled disk, convolved
    % with Gaussian. Ring: Projection of a ring, convolved with
    % Gaussian. Distance: Two Gaussians in a distance d.
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
                prof=lineprof(x-mx,y-my,p,axprof);
                linecorr(x-mx,y-my,p,axcorr)
                
            end


        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end

function linecorr(x,y,p,axcorr)
n=min(x):p.binwidth:max(x);
nplot=n(1:end-1)+p.binwidth/2;
prof = histcounts(x,  n);
ac=xcorr(prof,prof);
ach=ac(ceil(end/2):end)/max(ac);
nac=(0:length(ach)-1)'*p.binwidth;
[yRPeaks,xPeaksIdx]=findpeaks(detrend(ach,'linear'));
xPeaksIdx(xPeaksIdx<3)=[]; yRPeaks(xPeaksIdx<3)=[]; 
[~,imxp]=max(yRPeaks);
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
% period=mean(diff(ppos));
period=xRPeaks(1);
ff="%2.1f";

h2=plot(axcorr,nac,acmulticorr,'DisplayName',"multi-c: P="+num2str(period,ff));
legend(axcorr)



end

function prof=lineprof(x,y,p,axprof)
n=min(x):p.binwidth:max(x);
nplot=n(1:end-1)+p.binwidth/2;
prof = histcounts(x,  n);
plot(axprof,nplot,prof)
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