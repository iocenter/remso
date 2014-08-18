function [ stateMrst,Jac ] = stateVector2stateMrst( stateVector,varargin)
%
%  write a state vector as a mrst state structure
%
%

opt = struct('xScale',[],...
    'partials',false);
opt = merge_options(opt, varargin{:});


% TODO: considering oil water
if ~isempty(opt.xScale)
    stateVector = stateVector.*opt.xScale;
end

    
    nx = numel(stateVector)/2;
    
    stateMrst.pressure = stateVector(1:nx);
    stateMrst.s = [stateVector(nx+1:end),1-stateVector(nx+1:end)];
    
    if opt.partials
        if ~isempty(opt.xScale)
            Jac = bsxfun(@times,speye(2*nx),opt.xScale');
        else
            Jac = speye(2*nx);
        end
    else
        Jac = [];
    end
    
    
    


end