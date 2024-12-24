function [nu_n] = negative(ee)
    if ee <= -0.5 && ee >= -1
        nu_n = 1;
    elseif ee <= 0 && ee >= -0.5
        nu_n = -2 * ee;
    else
        nu_n = 0;
    end
end
