function [client_indices, ranking] = compute_client_rankings(ionet_data, sector_indices, n_sectors)
    n_analyzed = numel(sector_indices);
    client_indices = zeros(n_analyzed, 1);
    ranking = zeros(n_analyzed, n_sectors);

    for i = 1:n_analyzed
        s_idx = sector_indices(i);

        ionet_excl = [ionet_data(s_idx, 1:s_idx-1), ionet_data(s_idx, s_idx+1:end)];
        [~, col_index] = max(ionet_excl);
        if col_index >= s_idx
            col_index = col_index + 1;
        end
        client_indices(i) = col_index;

        shares_vector = ionet_data(s_idx, :);
        [~, sort_idx] = sort(shares_vector, 'descend');
        rank = zeros(1, n_sectors);
        rank(sort_idx) = 1:n_sectors;

        if rank(s_idx) ~= 1
            old_rank = rank(s_idx);
            rank(sort_idx(1)) = old_rank;
            rank(s_idx) = 1;
        end
        ranking(i, :) = rank;
    end
end
