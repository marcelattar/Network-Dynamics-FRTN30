%% Task 2b
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
P = sparse(inv(D)*W);

% From Task 2a) we know that the two most central nodes are 1 and 2. I
% therefore chose these as the stuborn nodes (not necessary).
Q = P(3:N,3:N);
E = P(3:N,1:2);
G = P(1:2,1:2);
F = P(1:2,3:N);

% Making a new vector with the correct indecies for the "user" vector (this
% is not used later)
new_users = [users(3:N);users(1:2)];

x_underscore = 1/2*ones(length(Q),1);
u = [1;0];
x = [x_underscore;u];

% Calculating values from different nodes.
node3 = zeros(500,1);
node1273 = zeros(500,1);
node6600 = zeros(500,1);
for time_step = 1:500
    %pause(0.1)
    x_underscore = Q*x_underscore + E*u;
    x = [x_underscore;u];
    node3(time_step) = x_underscore(3);
    node1273(time_step) = x_underscore(1273);
    node6600(time_step) = x_underscore(6600);
    %histogram(x)
end
%disp('Done')
figure(1)
histogram(x)
xlabel('Opinion value');
ylabel('Number of nodes')
title('Opinion distribution after 500 time steps')
set(gca, 'FontSize', 13);


figure(2)
plot(linspace(1,500,500),node3)
hold on
plot(linspace(1,500,500),node1273)
hold on
plot(linspace(1,500,500),node6600)
xlabel('Time steps');
ylabel('Opinion value')

set(gca, 'FontSize', 13);
legend('Random node 1', 'Random node 2','Random node 3');