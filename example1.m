function example1()
    
    bounds.a = 0;
    bounds.b = 1;
    bounds.c = 0;
    bounds.d = 1;

    T = 1;
    
    f_func = @(x, y, t) -2 * t * t * ((x-1) * x + (y-1) * y) + 2 * t * ((x-1) * x * (y - 1) * y);
    gamma_func = @(x, y) 0;
    beta_func.x = @(x, y) 0;
    beta_func.y = @(x, y) 0;
    g_func = @(x, y, t) 0;
    h_func = @(x, y, t) 0;
    q_func = @(x, y, t) 0;
    kappa = 1;

    initial_conditions(1).x = 0;
    initial_conditions(1).y = 0;
    initial_conditions(1).value = 0;

    bound_conditions(1).bound = "top";
    bound_conditions(1).condition_type = "value";
    bound_conditions(2).bound = "right";
    bound_conditions(2).condition_type = "value";
    bound_conditions(3).bound = "bottom";
    bound_conditions(3).condition_type = "value";
    bound_conditions(4).bound = "left";
    bound_conditions(4).condition_type = "value";

    solution = @(x, y, t) t * t * x * (x-1) * y * (y-1);

    ns = [10, 100, 1000];
    dt_explicit = 0.1;
    dt_implicit = 0.01;

    errors_explicit = [];
    errors_implicit = [];

    for i = 1 : numel(ns)

        n = ns(i);

        [x_e, y_e, t_e, u_e] = explicit_ivp(bounds, T, n, n, dt_explicit, initial_conditions, bound_conditions, kappa, beta_func, gamma_func, g_func, h_func, q_func, f_func);
        [x_i, y_i, t_i, u_i] = implicit_ivp(bounds, T, n, n, dt_implicit, initial_conditions, bound_conditions, kappa, beta_func, gamma_func, g_func, h_func, q_func, f_func);
        
        sol_values_ex = [];
        for x = x_e
            for y = y_e
                sol_values_ex = [sol_values_ex; u_e(y * n + x) - solution(x, y, T)];
            end
        end

        sol_values_im = [];
        for x = x_i
            for y = y_i
                sol_values_im = [sol_values_im; u_e(y * n + x) - solution(x, y, T)];
            end
        end

        errors_explicit = cat(1, errors_explicit, [sqrt(sum(u_e - sol_values_ex) ^ 2) / sqrt(sum(sol_values_ex) ^ 2)]);
        errors_implicit = cat(1, errors_implicit, [sqrt(sum(u_i - sol_values_im) ^ 2) / sqrt(sum(sol_values_im) ^ 2)]);

    end

    space_length = (bounds.b - bounds.a) * (bounds.d - bounds.c);
    h = arrayfun(@(n) space_length / (n - 1), ns);
    logH = log(h);
    
    logE_explicit = log(errors_explicit);
    [p_explicit] = polyfit(logH, logE_explicit, 1);
    px_explicit = polyval(p_explicit, logH);

    hf = figure();
    plot(logH, logE_explicit, "o", logH, px_explicit);
    title("Taxa de Convergência - Explícito");
    print("out/example1_convergence_explicit.png", "-dpng");
    
    logE_implicit = log(errors_implicit);
    [p_implicit] = polyfit(logH, logE_implicit, 1);
    px_implicit = polyval(p_implicit, logH);

    hf = figure();
    plot(logH, logE_implicit, "o", logH, px_implicit);
    title("Taxa de Convergência - Implícito");
    print("out/example1_convergence_implicit.png", "-dpng");

endfunction