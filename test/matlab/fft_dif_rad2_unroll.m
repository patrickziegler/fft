function x = fft_dif_rad2_unroll(x)

N = length(x);
m = N;

while m > 1
    m = m / 2;
    e0 = exp(-pi*1i/m);
    e1 = 1;
    for k = 0:m-1
        for i = 1:2*m:N-1
            buf = x(i + k + m);
            x(i + k + m) = (x(i + k) - buf) * e1;
            x(i + k) = x(i + k) + buf;
        end
        e1 = e1 * e0;
    end
end

x = bitrevorder(x);

end

% Alternative (see below):
% https://gist.github.com/Eight1911/b7cc6b8cc9645822f0d340c3a4b056e1

% function X = fft_dif_rad2_unroll(x)
%
% N = length(x);
% h = N / 2;
% X = zeros(size(x));
%
% s = N;
% sublen = 1;
% while sublen < N
%     s = s / 2;
%     for k = 0:2*s:N-1
%         w = exp(pi * 1i * k / N);
%         for i = 1:s
%             buf = w * x(i + k + s);
%             X(i + k/2) = x(i + k) + buf;
%             X(i + k/2 + h) = x(i + k) - buf;
%         end
%     end
%     x = X;
%     sublen = 2 * sublen;
% end
%
% X = x;
%
% end

% function X = fft_dif_rad2_unroll(x)
%
% N = length(x);
% h = N / 2;
% X = zeros(size(x));
%
% for s = N*2.^-(1:log2(N))
%     for k = 0:2*s:N-1
%         w = exp(pi * 1i * k / N);
%         for i = 1:s
%             buf = w * x(i + k + s);
%             X(i + k/2) = x(i + k) + buf;
%             X(i + k/2 + h) = x(i + k) - buf;
%         end
%     end
%     x = X;
% end
%
% X = x;
%
% end
