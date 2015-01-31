function [ varargout ] = ssFwdMemory(u,ssFunc,recoverSim,varargin)

% Wrapper function used to keep simualations in memory and avoid the need
% to call functions 2 times

opt = struct('replace',false,'gradFirst',false,'debug',true);
opt = merge_options(opt, varargin{:});


if nargin >= 6
    if opt.debug
        fid = fopen('deleteMe2Break.txt','r');
        if fid == -1
            fid = fopen('deleteMe2Break.txt','w');fclose(fid);
            keyboard;
        else
            fclose(fid);
        end
    end
end



[simVars] = recoverSim(u);

emptySim = false;
if isempty(simVars)
    emptySim = true;
end

if opt.replace || emptySim

    [f,gradf,conv,simVars] = ssFunc(u,'simVars',simVars);
    
    recoverSim(u,simVars); % save simulation

else
    [f,gradf,conv,~] = ssFunc(u,'simVars',simVars);
end

if opt.gradFirst
   varargout = {gradf,f,conv};
else
   varargout = {f,gradf,conv}; 
end






end

