if isfield(results, 'Single_tracking_analysis')
    var = results.Single_tracking_analysis;
else
    if isfield(results, 'Dual_tracking_analysis')
        var = results.Dual_tracking_analysis;
    else
        disp('Neither Single nor Dual tracking analysis field exists.');
    end
end

table_out=var{1};
for k=2:length(var)
    table_out=vertcat(table_out,var{k});
end