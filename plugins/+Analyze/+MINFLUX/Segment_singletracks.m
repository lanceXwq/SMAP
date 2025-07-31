%% This pluging finds trajectories and list to the ROI manager. Modified from Segment_cotracks.m by Lasse Drengk.

classdef Segment_singletracks<interfaces.DialogProcessor
    %Links molecules in consecutive frames for SPT analysis
    methods
        function obj = Segment_singletracks(varargin)        
            obj@interfaces.DialogProcessor(varargin{:}) ;
            obj.showresults = true;
        end
        function out = run(obj,p)
            
            out = runi(obj,p);
        end
        function pard = guidef(obj)
            pard = guidef(obj);
        end
       
    end
end

function out = runi(obj,p)
out = [];

% Parameters
minlenlocs = p.minlenlocs;
lennmstartend = p.lennmstartend; 
layers = find(obj.getPar('sr_layerson'));

[locs,~] = obj.locData.getloc({'xnm','ynm','time','tid','filenumber','thi'},'layer',layers,'position','all','grouping','ungrouped','removefilter','filenumber');

SE = obj.locData.SE;

% Initialize cell grid if empty
if isempty(SE.cells)
    cg = ROIManager.Segment.makeCellGrid;
    cg.attachLocData(obj.locData);
    cg.attachPar(obj.P);
    p = cg.getAllParameters;
    cg.run(p);
end

cellp = vertcat(SE.cells.pos);

for c = length(SE.cells):-1:1
    cellfn(c,1) = SE.cells(c).info.filenumber;
end

fnums = unique(locs.filenumber);

for ff = 1:length(fnums)
    fhh = fnums(ff);

    % get tracks in the channel of interest, e.g. thi==0
    tids0 = locs.tid(locs.thi == 0 & locs.filenumber == fhh);
    numloc = histcounts(tids0,1:max(tids0)+1);
    id0 = find(numloc > minlenlocs);

    for k = 1:length(id0)
        idh = locs.tid==id0(k) & locs.filenumber == fhh;
        xh = locs.xnm(idh); 
        yh = locs.ynm(idh);
        time = locs.time(idh);

        len = sqrt((xh(end)-xh(1)).^2+(yh(end)-yh(1)).^2);
        if len < lennmstartend
            continue
        end

        currentsite = interfaces.SEsites;
        currentsite.pos = [mean(xh), mean(yh), 0];  

        thisf = ~(cellfn == fhh);
        [~,cind] = min(sum((currentsite.pos(1:2)-cellp).^2,2)+thisf*1e9);

        currentsite.info.cell = SE.cells(cind).ID;
        currentsite.info.filenumber = fhh;
        currentsite.annotation.comments = num2str(time([1 end])/1000);

        SE.addSite(currentsite);
    end
end

SE.processors.preview.updateSitelist;
end
             
function pard = guidef(obj)

pard.minlenlocst.object = struct('String','Min length track (locs)','Style','text');
pard.minlenlocst.position = [2,1];
pard.minlenlocst.Width = 1.5;

pard.minlenlocs.object = struct('String','5','Style','edit');
pard.minlenlocs.position = [2,2.5];
pard.minlenlocs.Width = .5;

pard.lennmstartendt.object = struct('String','start-end (nm) >','Style','text');
pard.lennmstartendt.position = [2,3.5];
pard.lennmstartendt.Width = 1.;

pard.lennmstartend.object = struct('String','300','Style','edit');
pard.lennmstartend.position = [2,4.5];
pard.lennmstartend.Width = .5;

pard.plugininfo.description = sprintf('co-tracking analysis');
pard.plugininfo.type = 'ProcessorPlugin';
end 

