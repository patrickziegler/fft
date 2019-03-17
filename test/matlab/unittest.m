function unittest()
addpath('../../build/bindings/matlab');

N = 4^6;

fprintf('Using signal length: %d\n', N);

[t, x, f, X] = make_data(N);

% dlmwrite('data.csv', [t(:), x(:), f(:), X(:)], ...
%     'delimiter', ';', 'newline', 'pc');

funcs = {
    % 'fft_dif_rad2_recu', ...
    % 'fft_dif_rad4_recu', ...
    'fft_dif_rad2_unroll', ...
    'fft_dif_rad4_unroll', ...
    'fft_dif_rad2_c', ...
    'fft_dif_rad4_c'
    };

for i = 1:length(funcs)
    c = funcs{i};
    
    fprintf(['Testing ', c, ': ']);
    fh = str2func(c);
    
    tic;
    Xc = fh(x);
    dt = toc * 1e3;
    
    d = abs(abs(Xc(:)) - abs(X(:)));
    md = max(d);
    
    if md < 1e-3
        fprintf('OK\n -> %.3f ms\n\n', dt);
    else
        fprintf('FAILED\n -> %.3f ms, Max. diff.: %.2f\n\n', dt, md);
        plot_data(t, x, f, X, Xc);
        error('Faulty result detected (see plot)!');
    end
end

plot_data(t, x, f, X);
end

function [t, x, f, X] = make_data(N)
fs = 1e3;
f0 = 20;
f1 = 22;
f2 = 440;

t = ((1 / fs) * (0:(N-1)))';
x = cos(2*pi*f0*t);
x = x + 2 * cos(2*pi*f1*t);
x = x + cos(2*pi*f2*t);
% x = x + 0.2 * rand(size(x));

f = ifftshift((-N/2 : N/2 - 1) * (fs / N));

tic;
X = fft(x);
dt = toc * 1e3;
fprintf(' -> %.3f ms (built-in fft)\n\n', dt);
end

function plot_data(t, x, f, X, Xc)
f = fftshift(f);
X = fftshift(X);

figure();

subplot(2,1,1);
plot(t, x);
xlim([min(t), max(t)]);
ylim([1.1 * min(x), 1.1 * max(x)]);
grid on;

subplot(2,1,2);
semilogy(f, abs(X));
grid on;
zoom on;

if nargin > 4
    hold on;
    % Xc = fftshift(Xc);
    semilogy(f, abs(Xc));
end
end
