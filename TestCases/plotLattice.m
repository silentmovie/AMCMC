function out=plotLattice(rho)



% select frames to plot
% rows_to_plot(1:4) = [1, 1000, 2000, 2999];
% rows_to_plot(5:8) = [3000,6000,9500,14000];
rows_to_plot(1:4) = [3000, 8000, 13000, 18000];
rows_to_plot(5:8) = [28000,38000,48000,58000];
% rows_to_plot(9:24) = [22500:8500:150000];
num_plots = length(rows_to_plot);

% Create large enough figure for square subplots
figure('Units', 'pixels', 'Position', [100, 100, 600, 300]);

% Set up 6x4 layout with tight spacing
t = tiledlayout(2, 4, 'TileSpacing', 'compact', 'Padding', 'none');

for idx = 1:num_plots
    i = rows_to_plot(idx);
    data_grid = reshape(rho(i, 1:625), [25, 25]);

    nexttile;
    imagesc(data_grid');
    axis image off;  % 'axis image' ensures square pixels
    colormap('hot');
    title(sprintf('Iter %d', i), 'FontSize', 12);
end

% Optional: add a colorbar
% cb = colorbar;
cb.Layout.Tile = 'east';  % place on the right outside the grid