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
                axcorr=obj.initaxis('corr1D');
                hold(axcorr,'off')
                axprof=obj.initaxis('prof1D');
                hold(axprof,'off')             
                axfft=obj.initaxis('FFT');
                hold(axfft,'off')     
            end
            if ~p.setbinwidth
                p.binwidth=p.sr_pixrec;
            end
            layers=find(p.sr_layerson);
            k=1;
            [locs,~, hroi]=obj.locData.getloc({'xnm','ynm','znm','locprecnm','locprecznm','xnmline','ynmline'},'layer',layers(k),'position','roi');


            paircorrelationfunction(locs.xnm,locs.ynm,roihandle=obj.getPar('sr_roihandle'))


            ca=correlationtools(locs,p.binwidth,periodguess=p.period);
            ca.plot('profile',axis=axprof,color='r')
            ca.plot('autocorrelation',axis=axcorr,color='r');
            ca.plot('fft',axis=axfft,color='r');
            paircorrelationfunction(locs.xnm,locs.ynm,roihandle=getPar('sr_roihandle'))
           

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