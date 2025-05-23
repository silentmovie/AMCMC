% rng(8);
% rhohist(1,:) = ones(1,5)/5;
% inpo = randsample(5,5e5/625, true, rhohist(1,:));

function all_samples = testrng(T, N, k, master_seed)
    % T: outer iterations
    % N: inner parfor iterations
    % k: how many elements to randsample
    % master_seed: base seed for reproducibility

    all_samples = cell(T, N);

    for t = 1:T
        parfor i = 1:N
            % Derive a unique seed for (t, i)
            seed_ti = master_seed + 100000 * t + i;

            % Create a local RNG stream
            stream = RandStream('Threefry', 'Seed', seed_ti);

            % Use randsample with local stream
            all_samples{t, i} = randsample(stream, 1:100, k);
        end
    end
end

        