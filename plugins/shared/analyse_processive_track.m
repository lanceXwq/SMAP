function out=analyse_processive_track(x,y,time)

x=double(x);y=double(y);time=double(time-time(1));%angle=-double(angle)+pi;
angle=atan2(y(end)-y(1),x(end)-x(1));
fitfun=@(par,time) horzcat(par(1)+par(2)*cos(par(3))*time, par(4)+par(2)*sin(par(3))*time);
v0=sqrt((x(end)-x(1))^2+(y(end)-y(1))^2)/(time(end)-time(1));
starpar=[x(1),v0,angle,y(1)];
startxy=fitfun(starpar,time);

huber_loss = @(r, delta) (abs(r) <= delta) .* (0.5 * r.^2) + ...
                         (abs(r) > delta) .* (delta * (abs(r) - 0.5 * delta));
delta=1;
% c=1;
% cauchy_loss = @(r) (c^2 / 2) * log(1 + (r / c).^2);


data=horzcat(x,y);
dataf=data; timef=time;
% objective = @(params) cauchy_loss(fitfun(params, time) - data);

% Use lsqnonlin to solve the robust nonlinear least squares problem
options = optimoptions('lsqnonlin', 'Display', 'off');
for k=1:3
    objective = @(params) huber_loss(fitfun(params, timef) - dataf, delta);
    [fitp,r,residuals] = lsqnonlin(objective, starpar, [], [], options);
    r2d=sqrt(sum(residuals.^2,2));
    [m,stdrob,~,outlier]=robustMean(r2d(:));
    if length(timef)-length(outlier) > 7 && ~isempty(outlier)
        dataf(outlier,:)=[]; timef(outlier,:)=[];
    else 
        break
    end

end

out.x0=fitp(1);out.v=fitp(2);out.angle=fitp(3); out.y0=fitp(4);out.rmse=sqrt(mean(residuals(:).^2));
dr=diff(residuals)/sqrt(2);
[m,stdrob,~,outlier]=robustMean(dr(:));
out.gof=out.rmse/stdrob;
out.numpoints=length(timef);
out.fraction=length(timef)/length(time);
out.len=sqrt(sum((dataf(end,:)-dataf(1,:)).^2));


fitxy=fitfun(fitp,time);
    [xr,yr]=rotcoord(x,y,fitp(3));
    [xf,yf]=rotcoord(fitxy(:,1),fitxy(:,2),fitp(3));
out.plot.x=x; out.plot.y=y; out.plot.xfit=fitxy(:,1); out.plot.yfit=fitxy(:,2);
out.plot.xr=xr; out.plot.yr=yr; out.plot.xfitr=xf; out.plot.yfitr=yf;
out.plot.time=time;

end
