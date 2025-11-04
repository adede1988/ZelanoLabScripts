 function f_out = set_slice(f_in, x)
        m = size(x, 1);
        f_out = f_in;
        f_out(1:m, 1, :) = x;
    end