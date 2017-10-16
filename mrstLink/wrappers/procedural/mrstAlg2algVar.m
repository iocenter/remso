function [ v, n] = mrstAlg2algVar( wellSol, netSol,varargin )
%
%  Extract the algebraic variables from the well. 
%

opt = struct('vScale',[], 'nScale',[]);
opt = merge_options(opt, varargin{:});

v = [vertcat(wellSol.qWs);
    vertcat(wellSol.qOs);
    vertcat(wellSol.bhp)];


if ~isempty(opt.vScale)
    v = v./opt.vScale;
end

%
%  Extract the algebraic variables from the network. 
%
if ~isempty(netSol)
    n = netSol;
end

end

