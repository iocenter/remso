function [ fM ] = memorizeLastSimulation(u,simVars,debug)

% Save the results of a simulation.
% Usefull for IPOPT, to avoid calling same simulations

uM = u;
sM = simVars;

    function [lastSim] = memoryFunc(uT,simVarsT)
        
        lastSim = [];
        
        if nargin < 2
            % return last simulation if it is equal
            if all(cellfun(@(x1,x2)all(x1==x2),uT,uM))
                lastSim = sM;
            else
                if debug
                    save lastMemory u
                end
            end
            
        else
            % save this pair
            uM = uT;
            sM = simVarsT;
        end
    end

fM = @memoryFunc;


end