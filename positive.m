function [nu_p] = positive(ee)
    if ee >= 0 && ee <= 0.5
        nu_p = 2 * ee;
    elseif ee >= 0.5 && ee <= 1
        nu_p = 1;
    else
        nu_p = 0;
    end
end
