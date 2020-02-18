%% Part 1: Always run this part first. 
clear all

% Transistion rate matrix
Lambda = [0 2/5 1/5 0 0;
          0 0 3/4 1/4 0;
          1/2 0 0 1/2 0;
          0 0 1/3 0 2/3;
          0 1/3 0 1/3 0];

[M,N] = size(Lambda);

% Out-degree of each node (also called rate of the distribution)
w = Lambda*ones(M,1);

% Normalized transition rate matrix
D = diag(w);
P = D\Lambda;

w_star = max(w);

% Transition probability matrix
P_bar = zeros(M,N);
for i=1:M
    for j=1:N
       if i==j
           P_bar(i,j) = 0;
       else
           P_bar(i,j) = Lambda(i,j)/w_star;
       end
    end
end
Pii = ones(M,1) - P_bar*ones(M,1);
P_bar = P_bar + diag(Pii);


%% a) the average return time for a particle to leave node a (2) and return
start_pos = 2; % node a
iter = 10^6;

time_vec = zeros(iter,1); % This is keeping track of the time after each step
time_stamp = zeros(iter,1); % This is keeping track of the times we hit node a
next_pos = start_pos;
F = cumsum(P_bar,2);

% state matrix
x = zeros(5,iter); 
x(start_pos,1) = 1;

prev_val = 0;
for i=2:iter
    
    % time between two ticks, t_next
    u = rand; % uniform distribution from 0-1
    r = w_star;
    t_next = -log(u)/r;
    time_vec(i) = time_vec(i-1)+t_next; % Updates the time vector 
   
    next_pos = find(u<F(next_pos,:),1); % Returns the index of the next position (node)
    x(next_pos,i) = 1; % Update the state
    
    % Checking if I'm in the starting position, if so I note the time of
    % reaching it.
    if next_pos==start_pos
        time_stamp(i) = time_vec(i)-prev_val;
        prev_val = time_vec(i);
    end
end

Average_return_time = mean(time_stamp(find(time_stamp))) % Ans: 6.7418

%% b)
% The invariant distribution, pi_bar
[V,~] = eig(P_bar');
pi_bar = V(:,1);
pi_bar = pi_bar/sum(pi_bar); % Normalizing so that pi_bar sums to 1

% Return-time of node a
return_time = 1/(w(2)*pi_bar(2));
 
%% c) average time of moving from node o (1) to node d (5)

start_pos = 1; % node o
goal_pos = 5; % node d
iter = 10^2;

time_vec = zeros(iter,1);
time_stamp = zeros(iter,1);
next_pos = start_pos;
F = cumsum(P_bar,2);


% state matrix
x = zeros(5,iter); 
x(start_pos,1) = 1;

prev_val = 0;
bool = true;
for i=2:iter
    
    % time between two ticks, t_next
    u = rand; % uniform distribution from 0-1
    r = w_star;
    t_next = -log(u)/r;
    %time_vec(i) = time_vec(i-1)+t_next; % Updates the time vector 
    if bool 
        time_vec(i) = time_vec(i-1)+t_next;
        next_pos = find(u<F(next_pos,:),1); % Returns the index of the next position (node)
        x(next_pos,i) = 1;
    else % If I hit the goal_pos in my last run, then I will reset the next_pos to the start_pos.
        time_vec(i) = 0;
        next_pos = start_pos;
        x(next_pos,i) = 1; 
        bool = true;
    end
    
    % Checking if I'm in the goal position, if so I note the time of
    % reaching it.
    if next_pos==goal_pos
        time_stamp(i) = time_vec(i);
        bool = false;
    end
end

Average_return_time = mean(time_stamp(find(time_stamp))) % Ans: 8.7637


%% d) Theoretical hitting time from o (1) to d (5)

S = 5;
R = setdiff(1:5, S); %returns the values in 1:5 that are not in S
P_bar_mod = P_bar(R,R); %removing rows and columns beloging to node 5

B = ones(length(R),1);
A = (diag(ones(length(R),1))-P_bar_mod);
X = linsolve(A,B);

Hitting_time = X./w_star;
Hitting_time(1) % Hitting time from node o to node d.

