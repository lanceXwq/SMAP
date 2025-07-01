if isfield(results, 'Single_tracking_analysis')
    stvar = results.Single_tracking_analysis;
    track_summary=[];
    track_details=[];
    for k=1:length(stvar)
        track_summary=vertcat(track_summary,stvar{k}.summarytable);
        track_details=vertcat(track_details,stvar{k}.trackstable);
    end
end
if isfield(results, 'Dual_tracking_analysis')
    dtvar = results.Dual_tracking_analysis;
    track_summary_dual=[];
    for k=1:length(dtvar)
        track_summary_dual=vertcat(track_summary_dual,dtvar{k});
    end
end