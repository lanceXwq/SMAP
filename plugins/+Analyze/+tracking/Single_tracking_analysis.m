classdef Single_tracking_analysis<interfaces.DialogProcessor
    %Links molecules in consecutive frames for SPT analysis
    methods
        function obj=Single_tracking_analysis(varargin)        
            obj@interfaces.DialogProcessor(varargin{:}) ;
            obj.showresults=true;
        end
        function out=run(obj,p)
            
            out=runi(obj,p);
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function loadtiff(obj,a,b,field,sstr)
            oldf=obj.getSingleGuiParameter(field);
            if isempty(oldf)
                oldf=sstr;
            end

            [newf,newpath]=uigetfile(oldf);
            if newf
                obj.setGuiParameters(struct(field,[newpath newf]));
            end
        
        end
    end
end
function out=runi(obj,p)
out=[];
%filters
minlenframes=p.minlenframes;
maxd=p.maxd; % nm max distance to associate localizations

% cotracklength=p.cotracklength; %minimum data poits associated between tracks
% cotrackfraction = p.cotrackfraction; %minimum fraction of localizations associated

% test directed movement
aspectratio=p.aspectratio; %should be below this for directed motion.
lennmstartend=p.lennmstartend; % minimum endpoint- startpoint
lennmmin=0; %minimum of largest extenstion standard deviation (not so useful)
minvelocity=p.velocitymin; %nm/s
maxvelocity=p.velocitymax;

% visualizing co-tracks
% onlyprogressivecotracks=true; %only when both partners are progressive
% markintiffile=false;
% tiffile=fout;
% tiffile='/Users/ries/datalocal/2color_kinesin/25_50ms_561nm01_640nm02_600w52_676w37_1_MMStack_Default_combined.tif';
layers=find(obj.getPar('sr_layerson'));
% obj=g;
obj.locData.filter; %does this fix the bug?
[locs,indin]=obj.locData.getloc({'xnm','ynm','znm','xpix','ypix','frame','track_id','track_length','layer','filenumber'},'layer',layers,'position','roi','grouping','ungrouped');

% unique(locs.channel)

exposuretime=obj.locData.files.file(locs.filenumber(1)).info.exposure;

% XXX if exposure is wront, overwrite here
% exposuretime= 10

% test directed movement, calculate statistics
usetracks=unique(locs.track_id(locs.track_id>0));
trackstat.lennmstartend=zeros(max(usetracks),1);
trackstat.stdshort=zeros(max(usetracks),1);
trackstat.stdlong=zeros(max(usetracks),1);
trackstat.lenframe=zeros(max(usetracks),1);
% trackstat.partnertrackid=zeros(max(usetracks),1);
trackstat.velocity=zeros(max(usetracks),1);
locs.track_length_new=zeros(size(locs.xnm));
% trackstat.partnerids=zeros(max(usetracks),1);

for k=1:length(usetracks)
    iduset=usetracks(k);
    tind=find(locs.track_id==iduset);
    trackstat.lenframe(iduset)=length(tind);
    locs.track_length_new(tind)=trackstat.lenframe(iduset);
    if trackstat.lenframe(iduset)<2
        continue
    end

    xh=locs.xnm(tind);yh=locs.ynm(tind);fh=locs.frame(tind);
    trackstat.lennmstartend(iduset)=sqrt((xh(end)-xh(1)).^2+(yh(end)-yh(1)).^2);
        
    c = cov(xh-mean(xh), yh-mean(yh));
    [a, ev] = eig(c);
    [ev,ind] = sort(diag(ev));
    [xa, ya] = deal(a(1,ind(end)), a(2,ind(end)));
    trackstat.angle(iduset)=cart2pol(xa, ya);
    trackstat.stdshort(iduset)=real(sqrt(ev(1)));
    trackstat.stdlong(iduset)=real(sqrt(ev(2)));
    trackstat.velocity(iduset)=trackstat.lennmstartend(iduset)/(exposuretime/1000*trackstat.lenframe(iduset)); %nm/s
           
end

intrack=(locs.track_id>0);

longtracks=locs.track_length_new>minlenframes;



% indr=find((intrack&longtracks));
% locr.x=locs.xnm(indr);
% locr.y=locs.ynm(indr);
% locr.frame=locs.frame(indr);
% locr.track_id=locs.track_id(indr);
if isempty(locs.track_id)
    disp('please use the SimpleTracking plugin first')
end


usetracks=unique(locs.track_id(intrack&longtracks));

roi=obj.locData.files.file(locs.filenumber(1)).info.roi;
pixelsize=obj.locData.files.file(locs.filenumber(1)).info.cam_pixelsize_um*1000;


%look for progressive movement
validstats=true(size(trackstat.velocity));
validstats = validstats & trackstat.stdshort./trackstat.stdlong<aspectratio;
validstats = validstats & trackstat.lennmstartend > lennmstartend;
validstats = validstats & trackstat.stdlong > lennmmin;
validstats = validstats & trackstat.velocity > minvelocity;
validstats = validstats & trackstat.velocity < maxvelocity;

trackstat.progressive=validstats;

% plot tracks

ax=obj.initaxis('xy');
cols=[1 0 1
      0 1 1];
      % 1 0 0 
      % 0 0 1
      % .5 0.2 0.2
      % 0.2 0.2 .5
      % .7 0 0.7
      % 0 0.7 .7];

colind=(trackstat.progressive)+1;   

for k=1:length(usetracks)  
    idh=usetracks(k);
    indtr=locs.track_id==idh;
    
    msize=3;
    lw=1;
    symb='-';

    hp=plot(ax,locs.xnm(indtr)/pixelsize(1)-roi(1),locs.ynm(indtr)/pixelsize(2)-roi(2),symb,'Color',cols(colind(idh),:),'LineWidth',lw,'Tag','test','MarkerSize',msize);
    hold(ax,"on")
    % pidlabel=0*locs.track_id(indtr)+trackstat.partnerids(idh);
    dtRows = [dataTipTextRow("frame",double(locs.frame(indtr))),...
    dataTipTextRow("ID",double(locs.track_id(indtr)))];
    % dataTipTextRow("partnerID",double(pidlabel))];
    alldatatip=vertcat(hp.DataTipTemplate.DataTipRows,dtRows');
    %hp.DataTipTemplate.DataTipRows(end+1:end+3) = dtRows;   
    hp.DataTipTemplate.DataTipRows=alldatatip;
end


    axis(ax,'ij');
    axis(ax,'equal');

% Plot cotracks vs time
%%
%only both good
goodpairs=find(trackstat.progressive);
if contains(p.showtraces.selection,'progressive') 
    figure;
   
end
    % numrows=ceil(length(goodpairs)/5);
f=0;
v=zeros(length(goodpairs),1);runlength=v; runtime=v;
goodv=true(length(goodpairs),1);
plotind=1;
for k=1:length(goodpairs)
    % subplot(5,6,2*k-1-f)
    % hold off
    id1=locs.track_id==goodpairs(k);
    statv=analyse_progressive_track(locs.xnm(id1),locs.ynm(id1),locs.frame(id1));


        tmin=(min(locs.frame(id1)));
        tmax=(max(locs.frame(id1)));
        % 
        % [x1,y1]=rotcoord(locs.xnm(id1)-mean(locs.xnm(id1)),locs.ynm(id1)-mean(locs.ynm(id1)),trackstat.angle(goodpairs(k)));
        % % [x2,y2]=rotcoord(locs.xnm(id2)-mean(locs.xnm(id2)),locs.ynm(id2)-mean(locs.ynm(id2)),trackstat.angle(pid));
        % plot(locs.frame(id1),x1,'.-')
        % hold on
        
        % xlabel('time(frame)')
        % ylabel('xrot (nm)')
        % subplot(5,6,2*k-f)
        % plot(locs.xnm(id1)-mean(locs.xnm(id1)),locs.ynm(id1)-mean(locs.ynm(id1)),'.-')
        % axis equal
        % xlabel('x (nm)')
        % ylabel('y (nm)')
    exposuretimes=exposuretime/1000;
    % statv
    if statv.gof>3 || statv.numpoints<minlenframes || statv.v/exposuretimes<minvelocity
        goodv(k)=false;
    end
    
    v(k)=statv.v/exposuretimes;
    runlength(k)=statv.len;
    runtime(k)=statv.numpoints;
    ff='%2.0f';
    if contains(p.showtraces.selection,'progressive') && goodv(k)
        if 2*plotind-f>30
            f=f+30;
            figure
        end
        ax1=subplot(5,6,2*plotind-1-f);
        ax2=subplot(5,6,2*plotind-f);    
        plotind=plotind+1;

        plot(ax1,statv.plot.time,statv.plot.xr,statv.plot.time,statv.plot.xfitr);
        xlabel(ax1,'time')
        ylabel(ax1,'xrot')
        plot(ax2,statv.plot.x,statv.plot.y,statv.plot.xfit,statv.plot.yfit);
        axis(ax2,'equal')
        xlabel(ax2,'xrot')
        ylabel(ax2,'yrot')
% end
        title(ax2,['f: ' num2str(tmin) ':', num2str(tmax),', x: ' num2str(mean(locs.xpix(id1)),'%3.0f') ', y: ' num2str(mean(locs.ypix(id1)),'%3.0f')])  
        title(ax1,"Id:"+goodpairs(k)+ ", v: " + num2str(v(k),ff))
    end        
end

axv=obj.initaxis('v');
histogram(axv,v(goodv))
title(axv, "v nm/s, mean "+mean(v(goodv)) + ", median " + median(v(goodv)));
axr=obj.initaxis('runlength');
histogram(axr,runlength(goodv))
title(axr,"length nm, mean "+mean(runlength(goodv)) + ", median " + median(runlength(goodv)));
axrt=obj.initaxis('runtime');
histogram(axrt,runtime(goodv))
title(axrt,"length time frames, mean "+mean(runtime(goodv)) + ", median " + median(runtime(goodv)));

% extracts filename from file path:
filePath = string(obj.getPar('lastSMLFile'));
% Find all occurrences of the substring
% slash_indices = strfind(filePath, '/');
% % If the substring is found
% if ~isempty(slash_indices)
%     % Get the last occurrence index
%     lastSlashIndex = slash_indices(end);
%     % Extract the substring after the last occurrence
%     fileName = extractAfter(filePath, lastSlashIndex);
% else
%     fileName = filePath;
% end
goodvi=find(goodv);
% fn=locs.filenumber(locs.track_id==goodpairs(goodvi(1)));
% fileslm=obj.getPar('filelist_long').String{fn};
tablename=strrep(filePath,'_sml.mat','_tracks.csv'); %if this does not work, replace fileslm by filePath
[path,fileName]=fileparts(filePath); 

% old output:
    % output=(sprintf([num2str(length(trackstat.progressive)), '\t' num2str(sum(trackstat.progressive))]));
% new output:
    % outputo=(table(fileName,length(trackstat.lenframe), sum(trackstat.progressive), 'VariableNames', {'Filename','Total', 'Progressive'}))
    output=(table(fileName,sum(trackstat.lenframe>=minlenframes), sum(goodv), mean(v(goodv)), mean(runlength(goodv)), mean(runtime(goodv)),'VariableNames', {'Filename','Total', 'Progressive','Velocity','runlength','runtime'}));

disp(output)
out.summary=output;

outputtracks=(table(repmat(string(fileName),sum(goodv),1), goodpairs(goodv),(v(goodv)), (runlength(goodv)), (runtime(goodv)),'VariableNames', {'Filename','ID','Velocity','runlength','runtime'}));
writetable(outputtracks,tablename)
out.tracks=outputtracks;
% disp(outputtracks)
% clipboard('copy',output)

end
             
function pard=guidef(obj)

pard.minlenframest.object=struct('String','Min length track (frames)','Style','text');
pard.minlenframest.position=[1,1];
pard.minlenframest.Width=1.5;

pard.minlenframes.object=struct('String','5','Style','edit');
pard.minlenframes.position=[1,2.5];
pard.minlenframes.Width=.5;

pard.maxdt.object=struct('String','Max distance (nm)','Style','text');
pard.maxdt.position=[1,3];
pard.maxdt.Width=1.5;

pard.maxd.object=struct('String','200','Style','edit');
pard.maxd.position=[1,4.5];
pard.maxd.Width=.5;

pard.dmt.object=struct('String','directed movement:','Style','text');
pard.dmt.position=[2,1];

pard.aspectratiot.object=struct('String','aspect ratio <','Style','text');
pard.aspectratiot.position=[2,2];
pard.aspectratiot.Width=1.;
pard.aspectratio.object=struct('String','1','Style','edit');
pard.aspectratio.position=[2,2.7];
pard.aspectratio.Width=.5;


pard.lennmstartendt.object=struct('String','start-end (nm) >','Style','text');
pard.lennmstartendt.position=[2,3.5];
pard.lennmstartendt.Width=1.;
pard.lennmstartend.object=struct('String','300','Style','edit');
pard.lennmstartend.position=[2,4.5];
pard.lennmstartend.Width=.5;

% pard.cotracklenght.object=struct('String','co track length (frames) >','Style','text');
% pard.cotracklenght.position=[3,1];
% pard.cotracklenght.Width=1.5;
% pard.cotracklength.object=struct('String','4','Style','edit');
% pard.cotracklength.position=[3,2.5];
% pard.cotracklength.Width=.5;
% 
% pard.cotrackfractiont.object=struct('String','co track fraction >','Style','text');
% pard.cotrackfractiont.position=[3,3];
% pard.cotrackfractiont.Width=1.5;
% pard.cotrackfraction.object=struct('String','0','Style','edit');
% pard.cotrackfraction.position=[3,4.5];
% pard.cotrackfraction.Width=.5;

pard.velocityt.object=struct('String','velocity (nm/s) min','Style','text');
pard.velocityt.position=[4,1];
pard.velocityt.Width=1.5;
pard.velocitymin.object=struct('String','10','Style','edit');
pard.velocitymin.position=[4,2.5];
pard.velocitymin.Width=.5;
pard.velocitymt.object=struct('String','max','Style','text');
pard.velocitymt.position=[4,4];
pard.velocitymt.Width=0.5;
pard.velocitymax.object=struct('String','1500','Style','edit');
pard.velocitymax.position=[4,4.5];
pard.velocitymax.Width=.5;

pard.showt.object=struct('String','Show:','Style','text');
pard.showt.position=[5,1];
pard.showt.Width=0.5;
pard.showtraces.object=struct('String',{{'progressive','none'}},'Style','popupmenu');
pard.showtraces.position=[5,1.5];
pard.showtraces.Width=1.5;

% pard.makemovie.object=struct('String','make movie','Style','checkbox');
% pard.makemovie.position=[6,1];
% pard.makemovie.Width=1.5;
% 
% pard.addtracksmovie.object=struct('String','add tracks to movie','Style','checkbox');
% pard.addtracksmovie.position=[6,3];
% pard.addtracksmovie.Width=1.3;
% 
% pard.tiffilet.object=struct('String','tif:','Style','text');
% pard.tiffilet.position=[7,1];
% pard.tiffilet.Width=0.5;
% 
% pard.tiffile.object=struct('String','','Style','edit');
% pard.tiffile.position=[7,1.5];
% pard.tiffile.Width=3.;
% 
% pard.tiffload.object=struct('String','load','Style','pushbutton','Callback',{{@obj.loadtiff,'tiffile','*.tif'}});
% pard.tiffload.position=[7,4.5];
% pard.tiffload.Width=0.5;
% 
% pard.Tfilet.object=struct('String','trafo:','Style','text');
% pard.Tfilet.position=[8,1];
% pard.Tfilet.Width=0.5;
% 
% pard.Tfile.object=struct('String','','Style','edit');
% pard.Tfile.position=[8,1.5];
% pard.Tfile.Width=3.;
% 
% pard.Tfload.object=struct('String','load','Style','pushbutton','Callback',{{@obj.loadtiff,'Tfile','*.mat'}});
% pard.Tfload.position=[8,4.5];
% pard.Tfload.Width=0.5;
% 

pard.plugininfo.description=sprintf('co-tracking analysis');
pard.plugininfo.type='ProcessorPlugin';
end

