%% Generate girvan-newman benchmark
    % Using CommuniyTool_v091: 
    % https://jp.mathworks.com/matlabcentral/fileexchange/45867-community-detection-toolbox

%% Data path
data_dir = 'mat\GN-net'

%% Generate girvan-newman network
    % Size of 128 nodes, 4 communities for 32 nodes per community

N = 32; % number of nodes in each community
K = 4;  % number of communities
M = 1;  % number of generated network for each p_out/p_in rate
C = 16; % average degree

% community strengths rate
beta = (N-1) / (N*(K-1));
for i = 1:100
    r = i / 100;
    z_in  = C * beta / (r + beta);
    z_out = C * r / (r + beta);
    for m = 1:M
        [A,V] = GGGirvanNewman(N, K, z_in, z_out, 0);
        if sum(sum(A)==0) > 0
            error('Error. There is no obsolte node rate=%s index=%s', num2str(r), num2str(m));
        end
        A = uint8(A);
        V = uint8(V);
        
        filebase = strcat(data_dir, '/nodes_', num2str(N), '_group_', num2str(K), '_avgdeg_', num2str(C), '_rate_', num2str(r, '%3.2e'),  '_index_', num2str(m));
        adjfile = strcat(filebase, '_adj.mat');
        save(adjfile, 'A');
    end
end