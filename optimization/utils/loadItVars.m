function [u,x,xs,v,vs,simVars] = loadItVars(varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

opt = struct('dir','./','name','itVars','it',0,'r',0);

opt = merge_options(opt, varargin{:});

name = fullfile(opt.dir,opt.name);
if opt.r > 0
    name = [name,'_r',num2str(opt.r)];
end
if opt.it > 0
    name = [name,'_it',num2str(opt.it)];
end


loadVars = load(name,'u','x','xs','v','vs','simVars');
u = loadVars.u;
x = loadVars.x;
xs = loadVars.xs;
v = loadVars.v;
vs = loadVars.vs;
simVars = loadVars.simVars;



end

