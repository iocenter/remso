function [ state ] = statePopulateScalar(state, varargin )
%
% fill state mock object fields with 'values' option (a scalar).  If the cell option
% is given, modify state only on the cell positions


opt = struct('fields',[],'values',[],'cells',[]);
opt = merge_options(opt, varargin{:});

if isempty(opt.cells)  %% do it on all cells!
    for k = 1:numel(opt.fields)
        switch opt.fields{k}
            case 'pressure'
                state.(opt.fields{k}) = ones(size(state.(opt.fields{k})))*opt.values(k);
            case 's' %% water saturation
                state.(opt.fields{k}) = ones(size(state.(opt.fields{k}),1),1)*opt.values(k);
            otherwise
                error(['Cannot handle state : ', opt.fields{k}]);
        end
    end
else
    n = size(opt.cells,1);
    for k = 1:numel(opt.fields)
        switch opt.fields{k}
            case 'pressure'
                state.(opt.fields{k})(opt.cells) = ones(n,1)*opt.values(k);
            case 's' %% water saturation
                state.(opt.fields{k})(opt.cells) = ones(n,1)*opt.values(k);
            otherwise
                error(['Cannot handle state : ', opt.fields{k}]);
        end
    end
    
    
end


end

