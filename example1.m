function example1()
    f_func = @(x, y, t) -2 * t * t * ((x-1) * x + (y-1) * y) + 2 * t * ((x-1) * x * (y - 1) * y);
endfunction