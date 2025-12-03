classdef correlationtools<handle
    % ct=correlationtools(x,bin), then all is calculated, stored at
    % properties
    % ct.fft.f, ct.fft.mag ct.fft.peaks, ct.fft.mainpeak
    % ct.plot('fft',arguments)
    % cross-correlation: function taht takes two correlation objects (or input: additional one),
    % outputs cc structure
    %AC change to two-sided, then same for CC
    properties
        x
        bin
        fft
        profile
        autocorrelation
        periodguess
        sigmafilter
        maxcorr
        
    end
    methods
        function obj = correlationtools(posin,bin,options)
            arguments
                posin
                bin
                options.periodguess=[];
                options.sigmafilter=2;
                options.maxcorr=[];
            end
            obj.periodguess=options.periodguess;
            obj.sigmafilter=options.sigmafilter;
            obj.maxcorr=options.maxcorr;
            if isa(posin, 'struct')
                if isfield(posin,'xnmline')
                    obj.x=posin.xnmline;
                else
                    obj.x=posin.xnm;
                end
            else
                obj.x=posin;
            end
            obj.bin=bin;
            obj.calculateprofile;
            obj.calculatefft
            obj.calculateautocorrelation
            
        end
        function calculateprofile(obj)
            obj.profile.counts = histcounts(obj.x, obj.nx);  

            proffilt=smoothdata(obj.profile.counts,'gaussian',6*obj.sigmafilter);
            [ypos,xPeaksIdx]=findpeaks(proffilt);
            [yRPeaks,xRPeaks] = refinepeaks(proffilt,xPeaksIdx,obj.nxplot);

            obj.profile.nx=obj.nxplot;
            obj.profile.n=obj.nx;
            obj.profile.filter.counts=proffilt;
            obj.profile.filter.nx=obj.nxplot;
            obj.profile.filter.peaks.x=xRPeaks;
            obj.profile.filter.peaks.y=yRPeaks;
        end
        function calculateautocorrelation(obj)
            ac=xcorr(obj.profile.counts,obj.profile.counts)/mean(obj.profile.counts)^2;
            ach=ac(ceil(end/2):end)/max(ac);
            nac=(0:length(ach)-1)'*obj.bin;
            [yRPeaks,xPeaksIdx]=findpeaks(detrend(ach,'linear'));
            xPeaksIdx(xPeaksIdx<3)=[]; 
            [obj.autocorrelation.period,obj.autocorrelation.periodmag]=obj.selectpeak(nac(xPeaksIdx),ach(xPeaksIdx));

            obj.autocorrelation.mag=ach;
            obj.autocorrelation.lag=nac;
            obj.autocorrelation.nx=nac;
            obj.autocorrelation.peaks.x=nac(xPeaksIdx);
            obj.autocorrelation.peaks.y=ach(xPeaksIdx);

            acs=obj.autocorrelation.mag;
            for k=1:3
                acs=detrend(acs,'linear');
                acs=xcorr(acs-mean(acs),acs-mean(acs));
                acs=acs/max(acs);
                acs=acs(ceil(end/2):end);
            end

            [yRPeaks,xPeaksIdx]=findpeaks(acs);
            [yRPeaks,xRPeaks] = refinepeaks(acs,xPeaksIdx,nac);
            [obj.autocorrelation.filter.period,obj.autocorrelation.filter.periodmag]=obj.selectpeak(xRPeaks,yRPeaks);

            % obj.autocorrelation.filter.period=periodfilt;
            obj.autocorrelation.filter.peaks.x=xRPeaks;
            obj.autocorrelation.filter.peaks.y=yRPeaks;
            obj.autocorrelation.filter.mag=acs;
            obj.autocorrelation.filter.nx=nac;
            

            % acmulticorr=0.5*(1+acs(ceil(end/2):ceil(end/2)+length(nac)-1));
            % acmulticorr=acmulticorr-min(acmulticorr);
            % acmulticorr=acmulticorr/mean(acmulticorr)*mean(ach);
        end
        function calculatefft(obj)
            [f, mag] = fft_one_sided(obj.nx,obj.profile.counts);
            [yRPeaks,xPeaksIdx]=findpeaks(mag);
            [yRPeaks,xRPeaks] = refinepeaks(mag,xPeaksIdx,f);
            periods=1./xRPeaks;
            [obj.fft.period,obj.fft.periodmag]=obj.selectpeak(periods,yRPeaks);
            
            obj.fft.peaks.x=xRPeaks;
            obj.fft.peaks.y=yRPeaks;
            obj.fft.mag=mag;
            obj.fft.f=f;
            obj.fft.xplot=f;
        end
        function crosscorrelationi(obj)
        end
        function nx=nx(obj)
            nx=double(min(obj.x):obj.bin:max(obj.x))';
        end
        function nxp=nxplot(obj)
            nx=obj.nx;
            nxp=nx(1:end-1)+obj.bin/2;
        end
        function plot(obj,prop,options)
            arguments
                obj
                prop
                % options.filter=0;
                options.plotpeaks=true;
                options.axis=gca;
                options.color='b';
            end
            ff="%2.1f";
            switch prop
                case 'profile'
                    plot(options.axis,obj.profile.nx,obj.profile.counts,'Color',options.color,'DisplayName',"profile")
                    hold(options.axis,'on')
                    plot(options.axis,obj.profile.filter.nx,obj.profile.filter.counts,'--','Color',options.color,'DisplayName','filt')
                    if options.plotpeaks
                        xpos=obj.profile.filter.peaks.x;
                        dx=diff(xpos);
                        tpos=(xpos(1:end-1)+xpos(2:end))/2;
                        text(options.axis,tpos,max(obj.profile.filter.counts)+0*tpos,compose(ff,dx))
                        plot(options.axis,obj.profile.filter.peaks.x,obj.profile.filter.peaks.y,'o','Color',options.color)
                    end
                     xlabel(options.axis,"position (nm)")
                     ylabel(options.axis,"counts")
                case 'fft'
                    ft=obj.fft;
                    plot(options.axis, ft.f, ft.mag,'Color',options.color,'DisplayName',"P: "+num2str(ft.period,ff));
                    hold(options.axis,'on')
                    if options.plotpeaks
                        plot(options.axis,1/ft.period,ft.periodmag,'o','Color',options.color,'HandleVisibility','off')
                        plot(options.axis,ft.peaks.x,ft.peaks.y,'x','Color',options.color,'HandleVisibility','off')

                    end
                    legend(options.axis)
                    xlim(options.axis,[0 10/ft.period])
                    xlabel(options.axis,"frequency (1/nm)")
                    ylabel(options.axis,"magnitude")

                case 'autocorrelation'
                    ac=obj.autocorrelation;
                    h1=plot(options.axis,ac.nx,ac.mag,'Color',options.color,'DisplayName',"correlation, P="+num2str(ac.period,ff));
                    hold(options.axis,'on')


                    acf=ac.filter.mag;
                    off=min(acf);
                    acf=acf-off;
                    fac=mean(acf)*mean(ac.mag);
                    acf=acf*fac;
                    h2=plot(options.axis,ac.nx,acf,'--','Color',options.color,'DisplayName',"multi-c: P="+num2str(ac.filter.period,ff));
                    legend(options.axis)
                    if ~isempty(obj.maxcorr)
                        xlim(options.axis,[0 obj.maxcorr])
                    end

                    if options.plotpeaks
                        plot(options.axis,ac.peaks.x,(ac.peaks.y),'x','Color',options.color,'HandleVisibility','off')
                        plot(options.axis,ac.period,ac.periodmag,'o','Color',options.color,'HandleVisibility','off')
                        plot(options.axis,ac.filter.peaks.x,(ac.filter.peaks.y-off)*fac,'x','Color',options.color,'HandleVisibility','off')
                        plot(options.axis,ac.filter.period,(ac.filter.periodmag-off)*fac,'o','Color',options.color,'HandleVisibility','off')
                    end

                    xlabel(options.axis,"dx (nm)")
                    ylabel(options.axis,"g(r)")
                otherwise
                    disp('not implemented')
            end
           

        end
        function [peak,peaky]=selectpeak(obj,xpeaks,ypeaks)
            if ~isempty(obj.periodguess)
                [~,imxp]=min(abs(xpeaks-obj.periodguess));
            else
                [~,imxp]=max(ypeaks);
            end
            peak=xpeaks(imxp);
            peaky=ypeaks(imxp);
        end
           
    end
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