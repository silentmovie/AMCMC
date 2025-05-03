function out=plotLattice(rho)

% rows_to_plot = 3000:6250:150000;  % Total of 24 plots
% num_plots = length(rows_to_plot);
% 
% figure('Position', [100, 100, 1200, 800]);
% t = tiledlayout(6, 4, 'Padding', 'compact', 'TileSpacing', 'compact');  % remove gaps
% for idx = 1:num_plots
%     i = rows_to_plot(idx);
%     data_grid = reshape(rho(i, 1:625), [25, 25]);
% 
%     subplot(6, 4, idx);  % 4 rows, 6 columns
%     imagesc(data_grid');
%     colormap('hot');     % Customize colormap as needed
%     axis off;            % Hide axes for clarity
%     title(sprintf('%d', i), 'FontSize', 8);
% end
% 
% % Optional: add a colorbar for the last subplot
% % subplot(6, 4, 24);
% % colorbar;
% 
% set(findall(gcf, 'Type', 'axes'), 'FontSize', 18);
% 
% 
% end

rows_to_plot(1:4) = [1, 1000, 2000, 2999];
rows_to_plot(5:8) = [3000,6000,9500,14000];
rows_to_plot(9:24) = [22500:8500:150000];
% rows_to_plot = 6250:6250:150000;  % 24 frames
num_plots = length(rows_to_plot);

% Create large enough figure for square subplots
figure('Units', 'pixels', 'Position', [100, 100, 600, 1500]);

% Set up 6x4 layout with tight spacing
t = tiledlayout(6, 4, 'TileSpacing', 'compact', 'Padding', 'none');

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
cb = colorbar;
cb.Layout.Tile = 'east';  % place on the right outside the grid