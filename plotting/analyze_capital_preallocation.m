function [CapitalStats] = analyze_capital_preallocation(SimulDeterm, params, opts, ModData)
% ANALYZE_CAPITAL_PREALLOCATION Analyze capital allocation from perfect foresight simulation
%
% This function analyzes how capital is allocated across sectors in the
% perfect foresight (deterministic) simulation, showing log deviations
% from the deterministic steady state.
%
% INPUTS:
%   SimulDeterm - Deterministic simulation output (n_vars x T), in log levels
%   params      - Model parameters structure with:
%                 - n_sectors: Number of sectors (37)
%   opts        - Options structure:
%                 - plot_figure: Create bar plot (default: true)
%                 - save_figure: Save figure to file (default: false)
%                 - figures_folder: Folder for saving (default: 'output')
%                 - save_label: Label for saved file (default: '')
%                 - highlighted_sector: Sector to highlight (default: 1, Mining/Oil/Gas)
%   ModData     - Model data with steady state (for endostates_ss)
%
% OUTPUTS:
%   CapitalStats - Structure with:
%                  - avg_k_logdev: Average capital (log deviation) by sector
%                  - std_k_logdev: Std of capital (log deviation) by sector
%                  - sector_labels: Display labels for each sector
%                  - highlighted_sector_idx: Index of highlighted sector
%                  - highlighted_sector_avg: Average capital for highlighted sector

%% Set default options
if nargin < 3
    opts = struct();
end
if nargin < 4
    error('analyze_capital_preallocation:MissingModData', ...
        'ModData required to compute log deviations from steady state');
end

opts = set_default(opts, 'plot_figure', true);
opts = set_default(opts, 'save_figure', false);
opts = set_default(opts, 'figures_folder', 'output');
opts = set_default(opts, 'save_label', '');
opts = set_default(opts, 'highlighted_sector', 1);

n_sectors = params.n_sectors;
idx = get_variable_indices(n_sectors);

% Get SS capital values (log levels)
k_ss = ModData.endostates_ss(1:n_sectors);  % log(K_ss) for each sector

%% Extract capital from simulation
% Capital is in indices [1, n_sectors] (first n elements of state vector)
k_simul = SimulDeterm(idx.k(1):idx.k(2), :);  % n_sectors x T, in log levels

%% Compute log deviations from SS
k_logdev = k_simul - k_ss;  % n_sectors x T, log deviation from SS

%% Compute statistics
avg_k_logdev = mean(k_logdev, 2);  % Average log deviation for each sector
std_k_logdev = std(k_logdev, 0, 2);  % Std of log deviation for each sector

%% Get sector labels
all_labels = SectorLabel(1:n_sectors);
sector_labels = all_labels.display;

%% Build output structure
CapitalStats = struct();
CapitalStats.avg_k_logdev = avg_k_logdev;
CapitalStats.std_k_logdev = std_k_logdev;
CapitalStats.sector_labels = sector_labels;
CapitalStats.highlighted_sector_idx = opts.highlighted_sector;
CapitalStats.highlighted_sector_avg = avg_k_logdev(opts.highlighted_sector);
CapitalStats.highlighted_sector_label = sector_labels{opts.highlighted_sector};

%% Print summary
fprintf('\n  ┌─ Capital Preallocation Analysis (Perfect Foresight) ───────────┐\n');
fprintf('  │  Simulation periods: %d                                        │\n', size(k_simul, 2));
fprintf('  │  Highlighted sector: %d (%s)                         │\n', ...
    opts.highlighted_sector, CapitalStats.highlighted_sector_label);
fprintf('  └────────────────────────────────────────────────────────────────┘\n');

fprintf('\n  Average Capital (log deviation from SS):\n');
fprintf('    %-25s: %+8.5f\n', CapitalStats.highlighted_sector_label, CapitalStats.highlighted_sector_avg);

fprintf('\n  Top 5 sectors by average capital (log deviation):\n');
[sorted_avg, sort_idx] = sort(avg_k_logdev, 'descend');
for i = 1:min(5, n_sectors)
    fprintf('    %2d. %-20s: %+8.5f\n', i, sector_labels{sort_idx(i)}, sorted_avg(i));
end

fprintf('\n  Bottom 5 sectors by average capital (log deviation):\n');
for i = n_sectors:-1:max(1, n_sectors-4)
    fprintf('    %2d. %-20s: %+8.5f\n', n_sectors-i+1, sector_labels{sort_idx(i)}, sorted_avg(i));
end

%% Create bar plot
if opts.plot_figure
    fig = figure('Position', [100, 100, 1400, 600], 'Color', 'w');
    
    % Create bar chart
    bar_colors = repmat([0.3, 0.5, 0.7], n_sectors, 1);
    bar_colors(opts.highlighted_sector, :) = [0.85, 0.33, 0.10];
    
    b = bar(avg_k_logdev, 'FaceColor', 'flat');
    b.CData = bar_colors;
    
    hold on;
    
    % Add error bars for standard deviation
    errorbar(1:n_sectors, avg_k_logdev, std_k_logdev, 'k.', 'LineWidth', 1);
    
    % Add zero line
    plot([0.5, n_sectors+0.5], [0, 0], 'k--', 'LineWidth', 1.5);
    
    hold off;
    
    % Styling
    set(gca, 'XTick', 1:n_sectors);
    set(gca, 'XTickLabel', sector_labels);
    set(gca, 'XTickLabelRotation', 45);
    set(gca, 'FontSize', 10);
    set(gca, 'Box', 'on');
    set(gca, 'TickDir', 'out');
    
    xlabel('Sector', 'FontSize', 12);
    ylabel('Average Capital (log deviation from SS)', 'FontSize', 12);
    title(sprintf('Capital Preallocation Across Sectors (Perfect Foresight Simulation)\nHighlighted: %s (Sector %d)', ...
        CapitalStats.highlighted_sector_label, opts.highlighted_sector), ...
        'FontSize', 14, 'Interpreter', 'none');
    
    % Add grid
    grid on;
    set(gca, 'GridAlpha', 0.3);
    
    % Add legend
    legend_str = sprintf('%s: avg = %+.5f', CapitalStats.highlighted_sector_label, CapitalStats.highlighted_sector_avg);
    text(opts.highlighted_sector, avg_k_logdev(opts.highlighted_sector) + std_k_logdev(opts.highlighted_sector) + 0.002, ...
        sprintf('%+.4f', avg_k_logdev(opts.highlighted_sector)), ...
        'HorizontalAlignment', 'center', 'FontSize', 9, 'FontWeight', 'bold');
    
    % Tight layout
    ax = gca;
    ax.Position(2) = 0.20;
    ax.Position(4) = 0.70;
    
    %% Save figure if requested
    if opts.save_figure
        if ~exist(opts.figures_folder, 'dir')
            mkdir(opts.figures_folder);
        end
        
        if isempty(opts.save_label)
            filename = fullfile(opts.figures_folder, 'capital_preallocation.png');
        else
            filename = fullfile(opts.figures_folder, ['capital_preallocation_', opts.save_label, '.png']);
        end
        
        try
            saveas(gcf, filename);
            fprintf('\n  ✓ Figure saved: %s\n', filename);
        catch ME
            fprintf('\n  ⚠ Could not save figure: %s\n', ME.message);
        end
    end
end

end

