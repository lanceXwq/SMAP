if isfield(results, 'Single_tracking_analysis')
    var = results.Single_tracking_analysis;
else
    if isfield(results, 'Dual_tracking_analysis')
        var = results.Dual_tracking_analysis;
    else
        disp('Neither Single nor Dual tracking analysis field exists.');
    end
end

track_summary=[];
track_details=[];
for k=1:length(var)
    track_summary=vertcat(track_summary,var{k}.summary);
    track_details=vertcat(track_details,var{k}.tracks);
end