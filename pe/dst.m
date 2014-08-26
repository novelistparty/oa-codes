function [ y ] = dst( x )
%DST Discrete Sine Transform of x

N = length(x);

y = zeros(size(x));

for k = 1:N
    y(k) = sum(x'.*sin(pi/(N+1)*k.*(1:N))); 
end


end

