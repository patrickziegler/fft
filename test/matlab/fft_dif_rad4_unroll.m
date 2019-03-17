function x = fft_dif_rad4_unroll(x)

N = length(x);
m = N;

while m > 1
    m = m / 4;
    e0 = exp(-0.5*pi*1i/m);
    e1 = 1;
    for k = 0:m-1
        for i = 1:4*m:N-1
            x0 = x(i + k);
            x1 = x(i + k + m);
            x2 = x(i + k + 2*m);
            x3 = x(i + k + 3*m);
            x(i + k) = x0 + x1 + x2 + x3;
            x(i + k + 2*m) = (x0 - 1i * x1 - x2 + 1i * x3) * e1;
            x(i + k + 1*m) = (x0 - x1 + x2 - x3) * e1^2;
            x(i + k + 3*m) = (x0 + 1i * x1 - x2 - 1i * x3) * e1^3;
        end
        e1 = e1 * e0;
    end
end

x = bitrevorder(x);

end
