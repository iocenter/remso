function [ ] = initPool2013b( varargin )


opt = struct('poolSize',[],'restart',true);
opt = merge_options(opt, varargin{:});

if opt.restart
    if (matlabpool('size') ~= 0) 
        matlabpool close
    end
else
    if (matlabpool('size') ~= 0) 
        return
    end
    
end

if isempty(opt.poolSize)

[status,machineName] = system('hostname');
machineName = machineName(1:end-1);
if (matlabpool('size') == 0)
    if(strcmp(machineName,'tucker') == 1) 
        matlabpool open 6
    elseif(strcmp(machineName,'apis') == 1)
        matlabpool open 8
    elseif(strcmp(machineName,'honoris') == 1)
        matlabpool open 8
    elseif(strcmp(machineName,'beehive') == 1)
        matlabpool open 8
    elseif(strcmp(machineName,'ITK-D1000') == 1)
        matlabpool open 8
    elseif(strcmp(machineName,'dantzig') == 1)
        matlabpool open 8
    else
        machineName = machineName(1:end-1);
        if(strcmp(machineName,'service') == 1)  %% vilje
            matlabpool open 12
        else
            matlabpool open
        end
    end
    
else
	matlabpool('OPEN',opt.poolSize);
end    
    
end

