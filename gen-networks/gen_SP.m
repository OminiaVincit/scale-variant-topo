%% Generate sales-pardo benchmark
    % Reference from:  MSCD_Wav Toolbox available on http://www.gipsa-lab.fr/~nicolas.tremblay/index.php?page=downloads

%% Data path
data_dir = 'mat\SP-net'

%% Generate hierarchical network
M = 1; % number of realizations at each configuration
N = 640; % number of nodes
S3=10;S2=30;S1=120; N=640;

COM3=repmat([1:N/S3]',1,S3)'; COM3=COM3(:);
COM2=repmat([1:N/(S2+S3)]',1,S2+S3)'; COM2=COM2(:);
COM1=repmat([1:N/(S1+S2+S3)]',1,S1+S2+S3)'; COM1=COM1(:);
COM=[COM1,COM2,COM3];

for rho = 0.05:0.05:2.0
    for kbar = 16
        for m = 1:M
            [A]=SP_hbenchmark(COM, rho, kbar);
            if sum(sum(A)==0) > 0
                fprintf('Error. There is obsolte node rho=%s, bar=%s, index=%s', num2str(rho), num2str(bar), num2str(m));
            end
            A = uint8(A);
            
            filebase = strcat(data_dir, '/nodes_', num2str(N), '_rho_', num2str(rho, '%3.2e'), '_kbar_', num2str(kbar, '%3.2e'), '_index_', num2str(m));
            adjfile = strcat(filebase, '_adj.mat');
            save(adjfile, 'A');
        end
    end
end