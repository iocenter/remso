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

E = vertcat(netSol.E);
V = vertcat(netSol.V);

   
n = [vertcat(E.qwE);
    vertcat(E.qoE);
    vertcat(E.qgE);
    vertcat(V.pressure);
    ];


if ~isempty(opt.nScale)
    n = n./opt.nScale;
end


end

