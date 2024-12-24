function [nu_z] = zero(ee)
    if ee <= 0 && ee >= -0.5
        nu_z = 1 + 2 * ee;
    elseif ee >= 0 && ee <= 0.5
        nu_z = 1 - 2 * ee;
    else
        nu_z = 0;
    end
end
