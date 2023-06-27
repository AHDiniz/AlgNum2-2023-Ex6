function [x, y, t, u] = explicit_ivp(bounds, T, m, n, dt, initial_condition, bound_conditions)

    x = linspace(bounds.a, bounds.b, n);
    y = linspace(bounds.c, bounds.d, m);
    t = linspace(0, T, dt);

    for t_i : t
    end

endfunction