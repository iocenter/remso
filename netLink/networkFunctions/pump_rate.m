function [ rate ] = pump_rate(f, base_rate, base_freq)
%pump_rate: calculates the pump flow rate for a given frequency f using the
% knowledge of the flow rate for a known frequency.
%
%     rateScale =3.5;
    rateScale =1.0;
    rate = base_rate.*(f./base_freq);
    rate = rate.*rateScale;
    

end

