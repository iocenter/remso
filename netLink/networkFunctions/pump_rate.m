function [ rate ] = pump_rate(f, base_rate, base_freq)
%pump_rate: calculates the pump flow rate for a given frequency f using the
% knowledge of the flow rate for a known frequency.
%
    rate = base_rate.*(f./base_freq);

end

