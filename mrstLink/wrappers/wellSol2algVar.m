function [ v ] = wellSol2algVar( wellSol,varargin )
%
%  Extract the algebraic variables from the well.  Scale it according to
%  the scaling! if doScale is true!
%

opt = struct('doScale',false,'vScale',[]);
opt = merge_options(opt, varargin{:});

v = [vertcat(wellSol.qWs);
    vertcat(wellSol.qOs);
    vertcat(wellSol.bhp)];

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

