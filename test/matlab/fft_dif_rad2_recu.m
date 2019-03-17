function X = fft_dif_rad2_recu(x)

N = length(x);

if N > 2
    k = (0:N/2-1)';
    E = exp(-2*pi*1i*k/N);
    x = reshape(x(:), N / 2, 2);
    A = fft_dif_rad2_recu(x(:,1) + x(:,2));
    B = fft_dif_rad2_recu((x(:,1) - x(:,2)) .* E);
    X = reshape([A B]', N, 1);
else
    X = [x(1) + x(2); x(1) - x(2)];
end

end

% function X = fft_dif_rad2_recu(x)
%
% N = length(x);
% T = [1 1; 1 -1];
%
% if N > 2
%     k = (0:N/2-1)';
%     E = exp(-2*pi*1i*k/N);
%     x = reshape(x(:), N / 2, 2) * T;
%     A = fft_dif_rad2_recu(x(:,1));
%     B = fft_dif_rad2_recu(x(:,2) .* E);
%     X = reshape([A B]', N, 1);
% else
%     X = T * x(:);
% end
%
% end
