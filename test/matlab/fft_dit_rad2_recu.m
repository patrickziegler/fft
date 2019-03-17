function X = fft_dit_rad2_recu(x)

N = length(x);

if N > 2
    x = reshape(x, 2, N / 2)';
    U = fft_dit_rad2_recu(x(:,1));
    G = fft_dit_rad2_recu(x(:,2));
    k = 0:N/2-1;
    E = exp(-2*pi*1i*k/N);
    X = [G + U .* E, G - U .* E];
else
    X = x(2) + [x(1), -x(1)];
end

end
