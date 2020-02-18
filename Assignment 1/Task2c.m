%% Task 2c
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
end % Adding self-loops to node 1 since it had out-degree 0.

w = W*ones(N,1);
D = diag(w);
P = (inv(D)*W);

%%  Node selection.
% I chose the most and least central nodes (with respect to Page Rank).
%stuborn_node1 = 1; % The most central node
%stuborn_node2 = 17; % One of the least central nodes

% I chose the two least central nodes (with respect to Page Rank).
%stuborn_node1 = 17;
%stuborn_node2 = 39;

% I chose the 99th and 100th most central nodes (with respect to Page Rank).
stuborn_node1 = 355;
stuborn_node2 = 482;

%%
% Making a binary/logical vector to be used for extracting the right
% matrices later.
bin_matr = ones(N);
bin_matr(stuborn_node1,stuborn_node1) = 0;
bin_matr(stuborn_node2,stuborn_node2) = 0;
bin_vec = all(bin_matr);
 
% Calculating the matrices.
Q = P(bin_vec,bin_vec);
E = P(bin_vec,~bin_vec);
G = P(~bin_vec,~bin_vec);
F = P(~bin_vec,bin_vec);

% Making a new vector with the correct indecies for the "user" vector (this
% is not used later)
new_users = [users(bin_vec);users(~bin_vec)];

x_underscore = 1/2*ones(length(Q),1); % Here I chose the starting value for all the regualar nodes.
u = [1;0]; % Here I chose the values for the stuborn nodes.
x = [x_underscore;u];

for time_step = 1:500
    %pause(0.2)
    x_underscore = Q*x_underscore + E*u;
    x = [x_underscore;u];
    %histogram(x)
end
disp('Done')
histogram(x)
xlabel('Opinion value');
ylabel('Number of nodes')
title('Opinion distribution after 500 time steps')
set(gca, 'FontSize', 13);