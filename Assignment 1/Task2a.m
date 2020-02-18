%% Task 2a)
clear all
close all
load -ascii twitter.mat
load -ascii users.mat

W = spconvert(twitter);
W(6893,6893)=0; % Making W into a square matrix
[N,M] = size(W);

% Checking the out-degrees of all nodes
out_degree = W*ones(N,1);
for i =1:N 
    if out_degree(i) ==0
        W(i,i) = 1;
    end
end % Adding self-loops to all nodes that has out-degree 0 (in this case only node 1).

% Assigning values for beta and my.
Beta = 0.15;
my = ones(N,1);

% Calculating P
w = W*ones(N,1); % This calculates the out-degree of each node.
D = diag(w); % This puts the out-degree in a diagonal matrix.
P = sparse(inv(D)*W); % This the normalized Weight matrix (all rows sum to 1)

z = 0;
L = sparse(eye(N));

tic;
for k=0:100
    z = z + Beta*L*my;
    L = sparse(L*(1-Beta)*P');
end % Calculating the page rank value z.
toc;

% Top 5 central nodes
Y_max = maxk(z,5);
final_users = [users(Y_max(1,2)), users(Y_max(2,2)), users(Y_max(3,2)), users(Y_max(4,2)), users(Y_max(5,2))]; % I've used these user ID's and looked up their usernames in the given web address.
disp('Most central nodes (user IDs) according to PageRank algorithm: ')
disp('1 @gustavnilsson, 2 AVPapadopoulos, 3 @Asienfoset, 4 @Vikingafoset, 5 @bianca_grossi94') % These were looked up online.

% Top 3 least central nodes
Y_min = mink(z,3); % 17 is (one of) the least central node

% Two nodes that are semi-central (will be used in task 2c)
Y_semi_centr = maxk(z,100);
semi_centr_node1 = Y_semi_centr(99,2); % =node 355
semi_centr_node2 = Y_semi_centr(100,2); % =node 482