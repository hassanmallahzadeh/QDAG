function fftmatrix(n)
%w = sym('w');
%M = sym('a',n);
N = 2^n;
M = zeros(N,N);
w = exp(2*pi*1i/(N));
for i = 1 : N
    for j = 1 : N
%         el = syms(w^(i*j));
        M(i,j) = w^mod((i-1)*(j-1),N);
    end
end
M = round(M,2);
M = array2table(M)
end