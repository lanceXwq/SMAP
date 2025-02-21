function out=existnotempty(st,varargin)
out=true;
for k=1:length(varargin)
    if ~isfield(st,varargin{k})
        out=false;
        return
    end
    st=st.(varargin{k});
end
if isempyt(st)
    out=false;
end
end