function out=analyse_progressive_track(x,y,time,showfit)

x=double(x);y=double(y);time=double(time-time(1));%angle=-double(angle)+pi;
angle=atan2(y(end)-y(1),x(end)-x(1))
fitfun=@(par,time) horzcat(par(1)+par(2)*cos(par(3))*time, par(4)+par(2)*sin(par(3))*time);
v0=sqrt((x(end)-x(1))^2+(y(end)-y(1))^2)/(time(end)-time(1));
x0=[x(1),v0,angle,y(1)];
startxy=fitfun(x0,time);
[fitp,r,residuals]=lsqcurvefit(fitfun, x0, time, horzcat(x,y));
out.x0=fitp(1);out.v=fitp(2);out.angle=fitp(3); out.y0=fitp(4);out.rmse=sqrt(mean(residuals.^2));

if showfit
    figure(99)
    subplot(2,2,1);
    fitxy=fitfun(fitp,time);
    plot(time,x,time,fitxy(:,1),time,startxy(:,1));
    subplot(2,2,2);
    plot(time,y,time,fitxy(:,2),time,startxy(:,2));
    subplot(2,2,3);
    plot(x,y,fitxy(:,1),fitxy(:,2),startxy(:,1),startxy(:,2));
    axis equal
    [xr,yr]=rotcoord(x,y,fitp(3));
    [xf,yf]=rotcoord(fitxy(:,1),fitxy(:,2),fitp(3));
    subplot(2,2,4);
    plot(time,xr,time,xf);
    waitforbuttonpress
end

end
