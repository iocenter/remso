function [  ] = plotSolution( x,u,ss,obj,varargin)

opt = struct('simulate',true,'plotObj',true);
opt = merge_options(opt, varargin{:});

nx = numel(x{1});
nu = numel(u{1});

X = cell2mat(x');
U =  cell2mat(u');
% t = 0:size(U,2)-1;
% U = cell2mat(arrayfun(@(v)[v v],U,'UniformOutput',false));
% t = cell2mat(arrayfun(@(v)[v v],t,'UniformOutput',false));
% U = U(:,2:end-1);
% t = t(:,2:end-1);





if opt.simulate
    [~,~,~,~,xp,v,~] = simulateSystemSS(u,ss,[]);
    XS = cell2mat(xp');
end

figure

for k = 1:nx
    subplot(4+opt.plotObj,1,k)
    plot(X(k,:));
    if opt.simulate
        hold on;
        plot(XS(k,:),'x');
        hold off;
    end
    switch k
        case 1
            title('cost');
        otherwise
            title(['state  ',num2str(k-1)])
    end
    
end

for k = 1:nu
    subplot(4+opt.plotObj,1,nx+k)
    plot(U(k,:));
    title(['control  ',num2str(k)])
end

if opt.simulate && opt.plotObj
    f_stage = zeros(1,numel(xp));
    
    for k = 1:numel(xp)
        cik = callArroba(ss.ci,{k});
        [f_stage(k)] = callArroba(obj{k},{xp{k},u{cik},v{k}});
        
    end
    subplot(5,1,nx+nu+1)
    plot(f_stage)
    title('Objective accumulation')
end

end
