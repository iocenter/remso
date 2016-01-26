function [ lb,ub ] = scheduleBounds( schedules,varargin)
%
% fill a schedule mock object with maximums and minimums according to the
% optional parameters
%

opt     = struct('maxProd',struct('GRAT',inf,'ORAT',inf,'WRAT',inf,'LRAT',inf,'RESV',inf,'BHP',inf),...
    'minProd',struct('GRAT',0,'ORAT',0,  'WRAT',0,  'LRAT',0,  'RESV',0,  'BHP',0),...
    'maxInj',struct('RATE',inf,'RESV',inf,  'BHP',inf),...
    'minInj',struct('RATE',0,'RESV',0,  'BHP',0),...
    'useScheduleLims',false);


opt     = merge_options(opt, varargin{:});


n_sp = numel(schedules);
ub = schedules;
lb = schedules;
if opt.useScheduleLims
    for k = 1:n_sp
        for wI = 1:numel(schedules(k).control.W)
            schedules(k).control.W(wI).val =schedules(k).control.W(wI).lims.(schedules(k).control.W(wI).('type'));
        end
        
    end
end

for k = 1:n_sp
    
    
    [ vals,type,control] = schedule2CellControls(schedules(k));
    
    valsUB = vals;
    valsLB = vals;
    
    for j = 1:numel(vals) 
        for i = 1:numel(vals{j}) 
            switch type{j}{i}
                case {'inj',1}
                    switch control{j}{i}
                        case {'BHP'}
                            if opt.useScheduleLims
                                valsUB{j}(i) = min(opt.maxInj.(control{j}{i}),vals{j}(i));  
                            else
                                valsUB{j}(i) = opt.maxInj.(control{j}{i});
                            end
                            valsLB{j}(i) = opt.minInj.(control{j}{i});
                        case {'RATE'}
                            if opt.useScheduleLims
                                valsUB{j}(i) = min(opt.maxInj.(control{j}{i}),vals{j}(i));
                            else
                                valsUB{j}(i) = opt.maxInj.(control{j}{i});
                            end
                            valsLB{j}(i) = opt.minInj.(control{j}{i});
                        otherwise
                            error('control not considered: add control treatment!')
                    end
                    
                case {'prod',-1}
                    switch control{j}{i}
                        case {'BHP'}
                            if opt.useScheduleLims
                                valsLB{j}(i) = max(opt.minProd.(control{j}{i}),vals{j}(i));
                            else
                                valsLB{j}(i) = opt.minProd.(control{j}{i});
                            end
                            valsUB{j}(i) = opt.maxProd.(control{j}{i});
                        case {'LRAT'}
                            if opt.useScheduleLims
                                valsLB{j}(i) = max(-opt.maxProd.(control{j}{i}),vals{j}(i));
                            else
                                valsLB{j}(i) = -opt.maxProd.(control{j}{i});
                            end
                            valsUB{j}(i) = -opt.minProd.(control{j}{i});
                        case {'ORAT'}
                            if opt.useScheduleLims
                                valsLB{j}(i) = max(-opt.maxProd.(control{j}{i}),vals{j}(i));
                            else
                                valsLB{j}(i) = -opt.maxProd.(control{j}{i});
                            end
                            valsUB{j}(i) = -opt.minProd.(control{j}{i});
                        case {'WRAT'}
                            if opt.useScheduleLims
                                valsLB{j}(i) = max(-opt.maxProd.(control{j}{i}),vals{j}(i));
                            else
                                valsLB{j}(i) = -opt.maxProd.(control{j}{i});
                            end
                            valsUB{j}(i) = -opt.minProd.(control{j}{i});
                        case {'GRAT'}
                            if opt.useScheduleLims
                                valsLB{j}(i) = max(-opt.maxProd.(control{j}{i}),vals{j}(i));
                            else
                                valsLB{j}(i) = -opt.maxProd.(control{j}{i});
                            end
                            valsUB{j}(i) = -opt.minProd.(control{j}{i});                            
                        otherwise
                            error('control not considered: add control treatment!')
                    end
                otherwise
                    error(['Cannot handle well type: ', type{j}]);
            end
        end
    end
    [ ub(k)] = cellControls2Schedule( valsUB,ub(k));
    [ lb(k)] = cellControls2Schedule( valsLB,lb(k));
    
end

