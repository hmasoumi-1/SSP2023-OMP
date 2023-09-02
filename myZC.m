function z = myZC(N)

if mod(N,2)==0  % When the length of sequence is even
    tau = 1;
%     tau = 3;
    k = (0:(N-1)).';
    z = exp(-1j*pi*tau*(k.^2)/N);
else  % When the length of sequence is odd
    tau = 2;
    k = (0:(N-1)).';
    z = exp(-1j*pi*tau*(k.*(k+1))/N);
end

end