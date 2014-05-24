function [ ] = initPool2014a( varargin )


opt = struct('poolSize',[],'restart',false);
opt = merge_options(opt, varargin{:});

if opt.restart
    if (gcp('nocreate') ~= 0) 
        delete(gcp('nocreate'))
    end
else
    if (gcp('nocreate') ~= 0) 
        return
    end
    
end

if isempty(opt.poolSize)
    
    
    [status,machineName] = system('hostname');
    machineName = machineName(1:end-1);
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        if(strcmp(machineName,'tucker') == 1)
            parpool(12);
        elseif(strcmp(machineName,'apis') == 1)
            parpool(8);
        elseif(strcmp(machineName,'honoris') == 1)
            parpool(8);
        elseif(strcmp(machineName,'beehive') == 1)
            parpool(8);
        elseif(strcmp(machineName,'ITK-D1000') == 1)
            parpool(8);
        elseif(strcmp(machineName,'dantzig') == 1)
            parpool(12);
        else
            machineName = machineName(1:end-1);
            if(strcmp(machineName,'service') == 1)  %% vilje
                parpool(12);
            else
                parpool();
            end
        end
    end
    
else
    parpool(opt.poolSize);
end    
    
end

