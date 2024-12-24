function [dd] = negative2(nu)
    if nu == 0
        dd = 0;
    else
        A1 = (1 - nu / 2) * nu;
        A2 = (nu^2) / 4;
        d1 = -(nu / 4 + 0.5);
        d2 = -nu / 3;
        dd = (d1 * A1 + d2 * A2) / (A1 + A2);
    end
end
