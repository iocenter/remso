function [pvtw,pvdo,swof,sWcon] = collectReservoirInformation(reservoirP)


w = functions(reservoirP.fluid.BW);
pvtw = w.workspace{2}.pvtw;

w = functions(reservoirP.fluid.BO);
pvdo = w.workspace{2}.pvdo{1};

w = functions(reservoirP.fluid.krW);
swof = w.workspace{2}.swof{1};

sWcon = reservoirP.fluid.sWcon;


minMaxInitPressure = [min(reservoirP.state.pressure(:)),max(reservoirP.state.pressure(:))];

minMaxInitsW = [min(reservoirP.state.s(:,1)),max(reservoirP.state.s(:,1))];

wellData = {};
for wI = 1:numel(reservoirP.schedule.control(1).W)
    for cI = 1:numel(reservoirP.schedule.control)
        w = reservoirP.schedule.control(cI).W(wI);
        wellData = [wellData;{w.sign,w.type,w.val,w.lims}];
    end
end

end