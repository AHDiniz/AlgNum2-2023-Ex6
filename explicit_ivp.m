function [x, y, t, u] = explicit_ivp(bounds, T, m, n, dt, initial_conditions, bound_conditions, kappa, beta_func, gamma_func, g_func, h_func, q_func, f_func)

    N = n * m;

    h = [(bounds.b - bounds.a) / (n - 1), (bounds.d - bounds.c) / (m - 1)];
    x = linspace(bounds.a, bounds.b, n);
    y = linspace(bounds.c, bounds.d, m);
    t = linspace(0, T, T / dt);

    x_value = @(i) x(clamp(idivide(int32(i), int32(n), "fix") + 1, 1, n));
    y_value = @(i) y(clamp(mod(int32(i), int32(n)) + 1, 1, m));

    a_coeff = arrayfun(@(i) dt * gamma_func(x_value(i), y_value(i)) + 2 * kappa * (1 / (h(1) * h(1)) + 1 / (h(2) * h(2))), [1 : N]);
    b_coeff = arrayfun(@(i) dt * (-kappa / (h(1) * h(1))) - (beta_func.x(x_value(i), y_value(i)) / (2 * h(1))), [1 : N]);
    c_coeff = arrayfun(@(i) dt * (-kappa / (h(1) * h(1))) + (beta_func.x(x_value(i), y_value(i)) / (2 * h(1))), [1 : N]);
    d_coeff = arrayfun(@(i) dt * (-kappa / (h(2) * h(2))) - (beta_func.y(x_value(i), y_value(i)) / (2 * h(2))), [1 : N]);
    e_coeff = arrayfun(@(i) dt * (-kappa / (h(2) * h(2))) + (beta_func.y(x_value(i), y_value(i)) / (2 * h(2))), [1 : N]);

    A = zeros(N, N);

    A(1,1) = a_coeff(1);
    A(1,2) = c_coeff(1);
    A(1,1+n) = e_coeff(1);

    for I = 2 : n
        A(I,I-1) = b_coeff(I);
        A(I,I) = a_coeff(I);
        A(I,I+1) = c_coeff(I);
        A(I,I+n) = e_coeff(I);
    end

    for I = n + 1 : ((m-1)*n)
        A(I,I-n) = d_coeff(I);
        A(I,I-1) = b_coeff(I);
        A(I,I) = a_coeff(I);
        A(I,I+1) = c_coeff(I);
        A(I,I+n) = e_coeff(I);
    end

    for I = (((m-1)*n)+1):((m*n)-1)
        A(I,I-n) = d_coeff(I);
        A(I,I-1) = b_coeff(I);
        A(I,I) = a_coeff(I);
        A(I,I+1) = c_coeff(I);
    end

    A(N,N) = a_coeff(N);
    A(N,N-1) = b_coeff(N);
    A(N,N-1) = d_coeff(N);


    # Apply bound conditions:

    left_bound = [1 : n : (m - 1) * n + 1];
    right_bound = [n : n : N];
    top_bound = [(m - 1) * n + 2 : N - 1];
    bottom_bound = [2 : (n - 1)];

    for t_i = t
        
        f = zeros(N);
        for i = 1 : n
            for j = 1 : m
                f(i * n + j) = f_func(x(i), y(j), t_i);
            end
        end

        for bound_index = 1 : numel(bound_conditions)

            condition = bound_conditions(bound_index);
            target_bound = [];

            switch condition.bound
                case "top"
                    switch condition.condition_type
                        case "value"
                            for i = top_bound
                                A(i,:) = zeros(N, 1);
                                A(i,i) = 1;
                                f(i) = g_func(x_value(i), y_value(i), t_i);
                            end
                        case "derivative"
                            for i = top_bound
                                A(i,i) += e_coeff(i);
                                f(i) += e_coeff(i) * (h(2) / kappa) * h_func(x_value(i), y_value(i), t_i);
                                A(i,i+n) = 0;
                            end
                        case "mixed"
                            for i = top_bound
                                A(i,i) += e_coeff(i) * (1 - (h(2) * condition.beta_value) / condition.alpha_value);
                                f(i) -= e_coeff(i) * ((h(2) * q_func(x_value(i), y_value(i), t_i)) / condition.alpha_value) * q_func(x_value(i), y_value(i), t_i);
                                A(i,i+n) = 0;
                            end
                        otherwise
                            return;
                    end
                case "right"
                    switch condition.condition_type
                        case "value"
                            for i = right_bound
                                A(i,:) = zeros(N, 1);
                                A(i,i) = 1;
                                f(i) = g_func(x_value(i), y_value(i), t_i);
                            end
                        case "derivative"
                            for i = right_bound
                                A(i,i) += c_coeff(i);
                                f(i) += c_coeff(i) * (h(1) / kappa) * h_func(x_value(i), y_value(i), t_i);
                                A(i,i+1) = 0;
                            end
                        case "mixed"
                            for i = right_bound
                                A(i,i) += c_coeff(i) * (1 - (h(1) * condition.beta_value) / condition.alpha_value);
                                f(i) -= c_coeff(i) * ((h(1) * q_func(x_value(i), y_value(i), t_i)) / condition.alpha_value) * q_func(x_value(i), y_value(i), t_i);
                                A(i,i+1) = 0;
                            end
                        otherwise
                            return;
                    end
                case "bottom"
                    switch condition.condition_type
                        case "value"
                            for i = bottom_bound
                                A(i,:) = zeros(N, 1);
                                A(i,i) = 1;
                                f(i) = g_func(x_value(i), y_value(i), t_i);
                            end
                        case "derivative"
                            for i = bottom_bound
                                A(i,i) += d_coeff(i);
                                f(i) += d_coeff(i) * (h(2) / kappa) * h_func(x_value(i), y_value(i), t_i);
                                A(i,i-n) = 0;
                            end
                        case "mixed"
                            for i = bottom_bound
                                A(i,i) += d_coeff(i) * (1 - (h(2) * condition.beta_value) / condition.alpha_value);
                                f(i) -= d_coeff(i) * ((h(2) * q_func(x_value(i), y_value(i), t_i)) / condition.alpha_value) * q_func(x_value(i), y_value(i), t_i);
                                A(i,i-n) = 0;
                            end
                        otherwise
                            return;
                    end
                case "left"
                    switch condition.condition_type
                        case "value"
                            for i = left_bound
                                A(i,:) = zeros(N, 1);
                                A(i,i) = 1;
                                f(i) = g_func(x_value(i), y_value(i), t_i);
                            end
                        case "derivative"
                            for i = left_bound
                                A(i,i) += b_coeff(i);
                                f(i) += b_coeff(i) * (h(1) / kappa) * h_func(x_value(i), y_value(i), t_i);
                                A(i,i+1) = 0;
                            end
                        case "mixed"
                            for i = left_bound
                                A(i,i) += b_coeff(i) * (1 - (h(1) * condition.beta_value) / condition.alpha_value);
                                f(i) -= b_coeff(i) * ((h(1) * q_func(x_value(i), y_value(i), t_i)) / condition.alpha_value) * q_func(x_value(i), y_value(i), t_i);
                                A(i,i+1) = 0;
                            end
                        otherwise
                            return;
                    end
                otherwise
                    return;
            end

        end

        # Apply initial conditions:
        prev_u = zeros(N, 1);
        for condition = initial_conditions
            prev_u(clamp(condition.y * n + condition.x, 1, N)) = condition.value;
        end

        u = A * prev_u + dt * f;
        prev_u = u;
    end

endfunction