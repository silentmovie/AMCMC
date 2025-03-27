minnow = -1;
seednow = 1;
for seed = 1:1e8
    [~,~,minEig] = ID_Cn(seed,3);
    if minEig > minnow
        minnow = minEig;
        seednow = seed;
    end
end