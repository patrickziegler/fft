function X = fft_dif_rad4_recu(x)

N = length(x);

if N > 4
    k = (0:N/4-1)';
    E = exp(-2*pi*1i*k/N);
    x = reshape(x(:), N / 4, 4);
    A = fft_dif_rad4_recu(x(:,1) + x(:,2) + x(:,3) + x(:,4));
    B = fft_dif_rad4_recu((x(:,1) - 1i * x(:,2) - x(:,3) + 1i * x(:,4)) .* E);
    C = fft_dif_rad4_recu((x(:,1) - x(:,2) + x(:,3) - x(:,4)) .* E.^2);
    D = fft_dif_rad4_recu((x(:,1) + 1i * x(:,2) - x(:,3) - 1i * x(:,4)) .* E.^3);
    X = reshape([A B C D]', N, 1);
else
    X = [x(1) + x(2) + x(3) + x(4); x(1) - 1i * x(2) - x(3) + 1i * x(4); x(1) - x(2) + x(3) - x(4); x(1) + 1i * x(2) - x(3) - 1i * x(4)];
end

end

% function X = fft_dif_rad4_recu(x)
%
% N = length(x);
% T = [1 1 1 1; 1 -1i -1 1i; 1 -1 1 -1; 1 1i -1 -1i];
%
% if N > 4
%     k = (0:N/4-1)';
%     E = exp(-2*pi*1i*k/N);
%     x = reshape(x(:), N / 4, 4) * T;
%     A = fft_dif_rad4_recu(x(:,1));
%     B = fft_dif_rad4_recu(x(:,2) .* E);
%     C = fft_dif_rad4_recu(x(:,3) .* E.^2);
%     D = fft_dif_rad4_recu(x(:,4) .* E.^3);
%     X = reshape([A B C D]', N, 1);
% else
%     X = T * x(:);
% end
%
% end
