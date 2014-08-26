function [ x ] = idst( y )
%DST Discrete Inverse Sine Transform of y

N = length(y);

x = zeros(size(y));

for k = 1:N
    x(k) = 2/(N+1)*sum(y.*sin(pi/(N+1)*k.*(1:N)')); 
end


end