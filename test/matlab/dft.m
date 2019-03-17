function X = dft(x)

N = length(x);
k = (0:N-1)' * (0:N-1);
E = exp(-2*pi*1i*k/N);
X = E * x';

end

% function X = dft(x)
%
% N = length(x);
% X = zeros(1,N);
%
% for n = 0:(N-1)
%     for k = 0:(N-1)
%         X(n+1) = X(n+1) + x(k+1)*exp(1i*-2*pi*n*k/N);
%     end
% end
%
% end
