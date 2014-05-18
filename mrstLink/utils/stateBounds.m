function [b] = stateBounds(state,stateBound,varargin)
%
% fill a state mock object with the values at stateBound, if cells is
% provided change only on the given positions
%


opt = struct('cells',[]);
opt = merge_options(opt, varargin{:});


fields = fieldnames(state);


fieldsU = cell(2,1);
values = inf(2,1);

j = 1; 
for k = 1:numel(fields)
    switch fields{k}
        case 'pressure'
            values(j) = stateBound.(fields{k});
            fieldsU{j} =  fields{k};
            j = j+1;
        case 's' %% water saturation
            values(j) = stateBound.(fields{k});
            fieldsU{j} =  fields{k};
            j = j+1;
        case 'flux'
            state = rmfield(state,'flux');
        otherwise
            error(['Cannot handle state : ', fields{k}]);
    end
end

b = statePopulateScalar(state, 'fields',fieldsU,'values',values,'cells',opt.cells);


end
