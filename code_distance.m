% Define generator matrix G and prime p

% R-S code: capset in PG(2,q)

G = [1 1 1 1 1 1 1; 
     0 1 2 3 4 5 6; 
     0 1 4 2 2 4 1];
p = 7;

% 4-general set in PG(3,q)

% G = [1 0 0 0 1 1; 
%      0 1 0 0 1 2; 
%      0 0 1 0 1 3;
%      0 0 0 1 1 4];
% p = 5;

% 5-general set in PG(4,q)

% G = [1 0 0 0 0 1 1 1; 
%      0 1 0 0 0 1 2 6; 
%      0 0 1 0 0 1 3 5;
%      0 0 0 1 0 1 4 3;
%      0 0 0 0 1 1 5 4];
% p = 7;


% Compute all possible codewords by multiplying G with all possible vectors over the field
n = size(G, 1); % Number of rows in G
disp(n);
field_elements = 0:p-1;
field_vectors = dec2base(0:p^n-1, p) - '0'; % Generate all possible vectors over the field
disp(field_vectors);
codewords = mod(field_vectors * G, p); % Multiply field vectors by G (modulo p)
disp(codewords)

% Find rows where all elements are zero
nonzero_rows = any(codewords, 2);

% Filter out the rows where all elements are zero
codewords = codewords(nonzero_rows, :);
disp(codewords);

% Compute the weight (number of non-zero elements) for each codeword
weights = sum(codewords ~= 0, 2);

% Find the minimum weight
min_weight = min(weights);

fprintf('Minimum distance: %d\n', min_weight);