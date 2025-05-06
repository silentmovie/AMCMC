function [neighbor_v] = return_neighbor_Hyper_Cube(k, d)

% suppose 0 <= k <= 2^d-1
% k indicates a vertex of the hyper cube
% k has binary reps (........)2
% this function returns all integers (in base of 10)
% whose binary representation just differs 1 digit from
% that of k.

% e.g. if k = 10 (d = 4)
% k = (1010)2
% then all the neighbours are:
% (0010)2, (1110)2, (1000)2, (1011)2
% that are the numbers
% 2, 14, 8, 11

% this function returns a vector of d decimal numbers representing the
% vertices of the hyper cube

neighbor_v = zeros(d, 1);

binary_str_of_decimal_k = dec2bin(k);
len_binary_k = strlength(binary_str_of_decimal_k);


for i = (len_binary_k + 1):d
    neighbor = k + 2^(i-1);
    neighbor_v(i) = neighbor;
end

for i = 1:len_binary_k
    
    binary_digit = str2num(['uint8(',binary_str_of_decimal_k(len_binary_k + 1 - i),')']);

    if binary_digit == 0
        neighbor = k + 2^(i - 1);
    else
        neighbor = k - 2^(i - 1);
    end
    
    neighbor_v(i) = neighbor;

end

end








