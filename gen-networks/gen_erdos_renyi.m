%% Generate erdos renyi network model for topology data analysis
    % Ref. https://epfl-lts2.github.io/gspbox-html/doc/graphs/gsp_erdos_renyi.html
a = [pwd '/gspbox-0.7.5']
addpath(a)
addpath([a '/graphs'])
addpath([a '/utils'])

%% Generate network
data_dir = 'mat/ER-net';
N = 128; % number of nodes
M = 1;   % number of generations
pstep = 0.001;
pend = 0.1;
pstart = 0.02;

for m = 1:M
    for p = pstart:pstep:pend
        filebase = strcat(data_dir, '/nodes_', num2str(N), '_plink_', num2str(p, '%1.3f'),  '_index_', num2str(m));
        adjfile = strcat(filebase, '_adj.mat');
        
        G = gsp_erdos_renyi(N, p);
        A = G.A;
        A = full(A);
        A = uint8(A);
        save(adjfile, 'A');
    end
end
