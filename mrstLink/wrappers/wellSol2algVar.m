function [ v ] = wellSol2algVar( wellSol,varargin )
%
%  Extract the algebraic variables from the well.  Scale it according to
%  the scaling! if doScale is true!
%

opt = struct('doScale',false,'vScale',[]);
opt = merge_options(opt, varargin{:});

%nw = numel(wellSol);
% TODO: only for oil-Water systems


if isfield(wellSol,'bhp') % TODO: Confirm with MRST Developers bhp or pressure?
    v = [vertcat(wellSol.qWs);
        vertcat(wellSol.qOs);
        vertcat(wellSol.bhp)];
else
    v = [vertcat(wellSol.qWs);
        vertcat(wellSol.qOs);
        vertcat(wellSol.pressure)];
end

if opt.doScale
    if isempty(opt.vScale)
        [wellSolScale] = wellSolScaling(wellSol);
        [ vScale ] = wellSol2algVar( wellSolScale );
        v = v./vScale;
    else
        v = v./opt.vScale;
    end
end



end

