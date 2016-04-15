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

n = [netSol.qw;
    netSol.qo;
    netSol.qg;
    netSol.pV];

if ~isempty(opt.nScale)
    n = n./opt.nScale;
end


end

