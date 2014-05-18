function [state] = stateScaling(state,varargin)
%
% fill a state mock object with scaling information according to the
% options
%

opt = struct('pressure',barsa,'s',1);
opt = merge_options(opt, varargin{:});


fields = fieldnames(state);


fieldsU = cell(2,1);
values = ones(2,1);

j = 1; 
for k = 1:numel(fields)
    switch fields{k}
        case 'pressure'
            values(j) = opt.(fields{k});
            fieldsU{j} = fields{k};
            j = j+1;
        case 's' %% water saturation
            values(j) = opt.(fields{k});
            fieldsU{j} =  fields{k};
            j = j+1;
        case 'flux'
            state = rmfield(state,'flux');
        otherwise
            error(['Cannot handle state : ', fields{k}]);
    end
end

state = statePopulateScalar(state, 'fields',fieldsU,'values',values);


end