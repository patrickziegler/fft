function X = fft_dit_rad2_unroll(x)

X = bitrevorder(x);
n = length(X);
m = 1;

while m < n
    m = m * 2;
    wc = exp(2*pi*1i/m);
    for o = 0:m:n-1
        w = 1;
        for k = 1:m/2
            tmp = X(k + o) * w;
            X(k + o) = X(k + o + m/2) + tmp;
            X(k + o + m/2) = X(k + o + m/2) - tmp;
            w = w * wc;
        end
    end
end

end
