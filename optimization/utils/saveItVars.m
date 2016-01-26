function [  ] = saveItVars(u,x,xs,v,vs,simVars,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

opt = struct('dir','./','name','itVars','it',0,'r',0,'keepPreviousIt',false);

opt = merge_options(opt, varargin{:});

name = fullfile(opt.dir,opt.name);
if opt.r > 0
    name = [name,'_r',num2str(opt.r)];
end
if ~opt.keepPreviousIt
   delete([name,'_it*.mat'])
end
if opt.it > 0
    name = [name,'_it',num2str(opt.it)];
end

[pathstr,~,~] = fileparts(name);
if ~exist(pathstr,'dir')
    mkdir(pathstr)
end

save(name,'u','x','xs','v','vs','simVars')


end

