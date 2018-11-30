%% Generate watts strogatz network model for topology data analysis
    % https://jp.mathworks.com/help/matlab/examples/build-watts-strogatz-small-world-graph-model.html

%% Generate network
% Data path
data_dir = 'mat\WS-net';
M=1;   % number of generalizations
K = 16; % average degree
step = 0.01;
N = 128; % number of nodes

for m = 1:M
    for beta = 0.00:step:1.00
        G = WattsStrogatz(N, K/2, beta);
        A = adjacency(G);
        if sum(sum(A)==0) > 0
            error('Error. There is obsolte node N=%s K=%s beta=%s', num2str(N), num2str(K), num2str(beta));
        end
        A = full(A);
        A = uint8(A);

        filebase = strcat(data_dir, '/nodes_', num2str(N), '_K_', num2str(K), '_beta_', num2str(beta, '%1.2f'),  '_index_', num2str(m));
        adjfile = strcat(filebase, '_adj.mat');
        save(adjfile, 'A');
    end
end

        
        