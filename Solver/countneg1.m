function out = countneg1(alphat)

chunk_size = 1000;
num_chunks = int64(length(alphat) / chunk_size);
count_neg1 = zeros(num_chunks, 1);

for i = 1:num_chunks
    idx_start = (i-1)*chunk_size + 1;
    idx_end = i*chunk_size;
%     count_neg1(i) = sum(alphat(idx_start:idx_end) == -1);
    count_neg1(i) = sum(alphat(idx_start:idx_end) < 0.01);
end

out = count_neg1;