function [ fgx ,gfgx ] = dealSnoptSimulateSS( x,f,udim,debug )

if nargin >= 3
    if debug
        firstTime = false;
        firstTime = 1 == sngetStatus; % calling snopt for the first time don't save x and don't stop
        
        if ~firstTime
            save lastCall x
        end
        fid = fopen('deleteMe2Break.txt','r');
        if fid == -1
            fid = fopen('deleteMe2Break.txt','w');fclose(fid);
            if ~firstTime
                keyboard;
            end
        else
            fclose(fid);
        end
        
        
        fidOut = fopen('snoptLog.txt','a');
        
    end
else
    debug = false;
end

fgx = [];
gfgx = [];

% sngetNeedG has to be implemented in the SNOPT standard matlab interface
needGradients = 0 ~= sngetNeedG;

u = mat2cell(x,udim,1);
t1 = tic;
[fx,gfx,converged,simVars] = f(u,'gradients',needGradients);
tf1 = toc(t1);

if ~all(converged)
    % snsetStatus has to be implemented in the SNOPT standard matlab interface
    snsetStatus( -1 ) ; %% shorten lineSearch
    if debug
        displayMessage('f time: %g.2 s. Failed on %i. Gradients required %i\n',{tf1,find(~converged,1,'first'),needGradients},'fid',fidOut);
        fclose(fidOut);
    end
    return;
end


fgx = full(fx);
gfgx = cell2mat(gfx);


if debug
    displayMessage('Call time: %g.2 s. G computed %i\n',{tf1,needGradients},'fid',fidOut);
    fclose(fidOut);
end




end

