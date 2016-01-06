function [ schedule ] = relaxLimsInSchedule( schedule)


for k = 1:length(schedule.control)
    for w = 1:length(schedule.control(k).W)
        
        lims = struct();
        if schedule.control(k).W(w).sign == 1  %% injector
            lims.rate = inf;
            lims.bhp = inf;
        else  % producer
            lims.orat = -inf;
            lims.wrat = -inf;
            lims.grat = -inf;
            lims.lrat = -inf;
            lims.bhp  = -inf;
        end
        schedule.control(k).W(w).lims = lims;
    end
end


end

