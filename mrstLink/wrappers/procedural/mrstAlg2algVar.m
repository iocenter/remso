function [ v ] = mrstAlg2algVar( wellSol,netSol,varargin )
%
%  Extract the algebraic variables from the well. 
%

opt = struct('vScale',[]);
opt = merge_options(opt, varargin{:});

v = [vertcat(wellSol.qWs);
    vertcat(wellSol.qOs);
    vertcat(wellSol.bhp);
    netSol];

if ~isempty(opt.vScale)
    v = v./opt.vScale;
end




end

