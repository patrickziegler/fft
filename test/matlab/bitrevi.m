function idx = bitrevi(n)

if mod(log2(n),1) > 0
    error(['Input ' num2str(n) ' is not a power of 2!']);
end

idx = zeros(n, 1);
tmp = 0;

for i = 2:n
    m = n / 2;
    
    while bitand(tmp, m)
        tmp = bitxor(tmp, m);
        m = bitshift(m, -1); % shift right
    end
    
    tmp = bitor(tmp, m);
    idx(i) = tmp;
end

idx = idx + 1;

end
