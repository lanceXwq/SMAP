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
            out=[];
            if p.linecorrelation
                roih=obj.getPar('sr_roihandle');
                if ~contains(class(roih),'imline')
                    warning('ROI needs to be imline. cannot execute correlation analysis');
                else
                    axcorr=obj.initaxis('corr1D');
                    hold(axcorr,'off')
                    axprof=obj.initaxis('prof1D');
                    hold(axprof,'off')             
                    axfft=obj.initaxis('FFT');
                    hold(axfft,'off')     
                end
            end
            if p.pc
                axpcf=obj.initaxis('PCF');
                hold(axpcf,'off')     
            end
            if ~p.setbinwidth
                p.binwidth=p.sr_pixrec;
            end
            layers=find(p.sr_layerson);
            k=1;
            [locs,~, hroi]=obj.locData.getloc({'xnm','ynm','znm','locprecnm','locprecznm','xnmline','ynmline'},'layer',layers(k),'position','roi');


            % paircorrelationfunction(locs.xnm,locs.ynm,roihandle=obj.getPar('sr_roihandle'))
            if p.pc
                [grn,xx]=paircorrelationfunction(locs.xnm,locs.ynm,roihandle=roih,dr=p.binwidth,rmax=p.maxr);
                plot(axpcf,xx(2:end),grn(2:end)); 
                xlabel(axpcf,'r (nm)')
                ylabel(axpcf,'pair correlation function (norm)')
            end
            if p.linecorrelation && contains(class(roih),'imline')
                ca=correlationtools(locs,p.binwidth,periodguess=p.period,maxcorr=p.maxr);
                ca.plot('profile',axis=axprof,color='r')
                ca.plot('autocorrelation',axis=axcorr,color='r');
                ca.plot('fft',axis=axfft,color='r');
            end

           

        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end



function pard=guidef(obj)
pard.inputParameters={'sr_layerson','linewidth_roi','sr_pixrec'};

 p(1).value=0; p(1).on={}; p(1).off={'binwidth'};
p(2).value=1; p(2).on={'binwidth'}; p(2).off={};
pard.setbinwidth.object=struct('String','binwidth (nm) (empty: pixelsize):','Style','checkbox','Value',1,'Callback',{{@obj.switchvisible,p}});
pard.setbinwidth.position=[1,1];
pard.setbinwidth.Width=3;

pard.binwidth.object=struct('String','10','Style','edit');
pard.binwidth.position=[1,3];
pard.binwidth.Width=0.5;
pard.binwidth.TooltipString='Binwidth for profiles. If not checked, use pixel size of reconstruction';
pard.setbinwidth.TooltipString=pard.binwidth.TooltipString;

pard.linecorrelation.object=struct('String','1D correlation','Style','checkbox','Value',1);
pard.linecorrelation.position=[2,1];
pard.linecorrelation.Width=1;

pard.pc.object=struct('String','2D pair correlation','Style','checkbox','Value',0);
pard.pc.position=[2,2];
pard.pc.Width=1;


pard.periodt.object=struct('String','approximate period (nm)','Style','text');
pard.periodt.position=[3,1];
pard.periodt.Width=1.5;
pard.period.object=struct('String','200','Style','edit');
pard.period.position=[3,2.5];
pard.period.Width=0.5;

pard.maxrt.object=struct('String','max correlation (nm)','Style','text');
pard.maxrt.position=[3,3];
pard.maxrt.Width=1.5;
pard.maxr.object=struct('String','1000','Style','edit');
pard.maxr.position=[3,4.5];
pard.maxr.Width=0.5;
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