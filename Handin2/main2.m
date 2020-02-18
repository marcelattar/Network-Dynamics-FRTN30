%% Part 2: Always run this part first. 
clear all
clc

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

%% a) Particle perspective
start_pos = 2; % node a
iter = 10000;
nbr_of_part = 100;
pre_all = 50; % This is just a number to pre-allocate a vector (have to be quite big)

time_stamp = zeros(iter,1); % This is keeping track of the times we hit node a
F = cumsum(P_bar,2);


for i=1:iter
    
    time_vec = zeros(pre_all,1); % This is keeping track of the time after each step
    next_pos = start_pos;
    
    % state matrix
    x = zeros(5,pre_all); 
    x(start_pos,1) = 1;
    
    particle = 1;
    time_particles = zeros(nbr_of_part,1); % This will record the hitting time for all 100 particles.
    j = 2;
    while particle <= nbr_of_part
        % time between two ticks, t_next
        u = rand; % uniform distribution from 0-1
        r = w_star;
        t_next = -log(u)/r;
        time_vec(j) = time_vec(j-1)+t_next; % Updates the time vector 
   
        next_pos = find(u<F(next_pos,:),1); % Returns the index of the next position (node)
        x(next_pos,j) = 1; % Update the state
    
        % Checking if I'm in the starting position, if so I note the time of
        % reaching it.
        if next_pos==start_pos
            % I note that the particle has reached a.
            time_particles(particle) = time_vec(j);
            particle = particle + 1; 
            
            % Resetting
            time_vec = zeros(pre_all,1); 
            x = zeros(5,pre_all); 
            x(start_pos,1) = 1;
            j = 1; 
        end
        j = j+1;
    end

    % Here I note the shortest amount of time it took for the 100 particles
    % to reach a (i.e. the fastest particle), I then note that time in 
    % "time_stamp".
    min_time = min(time_particles);
    time_stamp(i) = min_time;
end
Average_return_time = mean(time_stamp) % Ans: 1.6024

%% b) Node perspective ATTEMPT NBR 1 (WRONG)
close all

% Setting starting position and nbr of particles
start_pos = 1; % node o
nbr_of_particles = 100;

% Creating a time vector and a variable to keep track on a particles
% position.
next_pos = start_pos;
time_vec = zeros(1,1); % This is keeping track of the time after each step

% A matrix that records how many particles there are at each node for every
% time-step.
node_vec = zeros(5,1);
node_vec(start_pos,1) = nbr_of_particles; 

% This is used to find the next node to move a particle to
F = cumsum(P,2); % Do NOT use P_bar here

% Creating a counter
i = 1;
while time_vec(i)<=60
    i=i+1;
    % Uniform distribution vector
    u_vec = rand(5,1);
    
    % Rate vector
    rate_vec = node_vec(:,i-1).*w;
    
    % Time vector (the time that the first particle from each node moves)
    t_next_vec = -log(u_vec)./rate_vec;
    [t_next,remove] = min(t_next_vec);
    
    time_vec(i) = time_vec(i-1) + t_next;
    
    % Copying the old particle placements
    node_vec(:,i) = node_vec(:,i-1);
    
    % Removing the particle from the previous node
    node_vec(remove,i) = node_vec(remove,i) - 1;
    
    % Finding what node to move the particle to
    add = find(u_vec(remove)<F(remove,:),1); % Returns the index of the next position (node)
    
    % Adding a particle to the node
    node_vec(add,i) = node_vec(add,i) + 1;
end

% Calculating the average number of particles in each node
Avg_nbr_of_particles = mean(node_vec,2)

% Ploting the particles in the 5 different nodes
plot(time_vec,node_vec)
xlabel('Time unit');
ylabel('Number of particles')
lgd = legend('Node o','Node a','Node b','Node c','Node d');
lgd.FontSize = 15;
title('Particle movement')
set(gca, 'FontSize', 13);

%% b) Node perspective, ATTEMPT NBR 2 (CORRECT)
close all

iter = 1000;
last_value_in_nodes = zeros(5,iter);

for l=1:iter

% Setting starting position and nbr of particles
start_pos = 1; % node o
nbr_of_particles = 100;

% Creating a time vector
time_vec = zeros(1,1); % This is keeping track of the time after each step

% Create a vector to keep track on a particles position (each index 
% represents a particle and the value represent what node it's in, 1-5).
particle_vec = ones(nbr_of_particles,1); % they all start in node 1

% Creating a matrix that keeps track on how many particles there are in
% every node for each time step (this will be used to calc. the average).
part_per_node = zeros(5,1);
part_per_node(start_pos) = nbr_of_particles;

% This is used to find the next node to move a particle to
F = cumsum(P_bar,2); % Important to use BAR

% Creating a counter
i = 1;
while time_vec(i)<=60
    i=i+1;
    
    % Time of next step (when first particle moves)
    u = rand;
    rate = nbr_of_particles*w_star;
    t_next = -log(u)./rate;
    
    time_vec(i) = time_vec(i-1) + t_next;
    
    % Deciding what random particle to move.
    remove = randi(nbr_of_particles);
    
    % Looking up where that particle currently is
    old_node_placement = particle_vec(remove);
    
    % Deciding where it should be moved next
    next_pos = find(u<F(old_node_placement,:),1); % Returns the index of the next position (node)
    
    % Updating the particle vector
    particle_vec(remove) = next_pos;
    
    % Calculating how many particles there are per node and saving it for
    % each time step
    for j=1:5
        part_per_node(j,i) = sum(particle_vec==j);
    end
end
last_value_in_nodes(:,l) = part_per_node(:,end);
end
% Calculating the average number of particles in each node
%Avg_nbr_of_particles = mean(part_per_node,2)
Avg_nbr_of_particles = mean(last_value_in_nodes,2)


% Ploting the particles in the 5 different nodes
plot(time_vec,part_per_node)
xlabel('Time unit');
ylabel('Number of particles')
lgd = legend('Node o','Node a','Node b','Node c','Node d');
lgd.FontSize = 15;
title('Particle movement')
set(gca, 'FontSize', 13);


%% 

