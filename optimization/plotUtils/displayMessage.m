function [  ] = displayMessage( format,vars,varargin )


opt = struct('fid',[],'std',true);
opt = merge_options(opt, varargin{:});

if ~isempty(opt.fid)
    fprintf(opt.fid,format,vars{:});
end
if opt.std
    fprintf(format,vars{:});
end
    



end

