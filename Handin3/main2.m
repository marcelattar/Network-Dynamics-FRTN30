%% 2a. Coloring
clear 
clc

W = zeros(10);
W(1,2)=1;
W(2,1)=1; W(2,3)=1; W(3,2)=1; W(3,4)=1; W(4,3)=1; W(4,5)=1;
W(5,4)=1; W(5,6)=1; W(6,5)=1; W(6,7)=1; W(7,6)=1; W(7,8)=1;  
W(8,7)=1; W(8,9)=1; W(9,8)=1; W(9,10)=1; W(10,9)=1;

nbrIter =1000;
a=1; b=1; c=0; d=0; %cost values; cost decreases when neighbor node is choosing the same color
phi = [a d; c b]; %u_i(x_i , y_i) = sum(W_ij*phi(x_i , x_j)
nbrNodes = 10;
actions = [1 2]; %actions; 1 = Red, 2 = Green
nbrActions = length(actions);
x = ones(10,1); %initially every node is Red
potential = zeros(1,nbrIter); %potential for each time
coord = [0 0; 1 0; 2 0; 3 0; 4 0; 5 0;
        6 0; 7 0; 8 0; 9 0; 10 0];
        
% Draw graph 
figure
set(gcf,'color','white')

% Plot the graph and mark the node that the particle is in with green
subplot(211)
gplot(W,coord,'-k');

hold on
for i = 1:nbrNodes
    if x(i, 1) == 1
        scatter(coord(i,1),coord(i,2),200,'markeredgecolor','k','markerfacecolor', 'r');
    else % if node = 2 (=Green)
        scatter(coord(i,1),coord(i,2),200,'markeredgecolor','k','markerfacecolor', 'g');
    end

end
set(gca,'xtick',[],'ytick',[],'xcolor','w','ycolor','w')

 
close all
for k=1:nbrIter % loop during which each node chooses a new action
    eta = k/100;
    
    % Choosing random node
    node = randi(10);
    
    
        prob = zeros(nbrActions,1); %probabilities of actions
        for action = actions %loop through all possible actions
            c = 0; %utility/cost function for node if action is "action"
            for j = find(W(node,:) ~= 0) %all elements in a row ~= 0, => neighbor nodes
                c = c + phi(action,x(j)); % add 1 or 2 to utility if neighbor color equals color and add 0 otherwise
            end
            prob(action) = exp(-eta*c); %non-normalized prob
        end
        prob = prob/sum(prob); % normalized elements to probabilities
        F = cumsum(prob); % cumulative probability vector
        newAction = find((F>rand(1)),1); % choose new action according to probabilies
        x(node) = newAction; %will be used to draw new colors in the nodes
   
    
    %Compute total potential
    for i=1:nbrNodes
        for j = find(W(i,:) ~= 0) %loop through out-neighbors
            potential(k) = potential(k) + (1/2)*phi(x(i),x(j));
        end
    end

end

% Draw nodes
figure(1)
gplot(W,coord,'-k');
hold on
for i = 1:nbrNodes
    if x(i, 1) == 1
        scatter(coord(i,1),coord(i,2),200,'markeredgecolor','k','markerfacecolor', 'r');
    else % if node = 2 (=Green)
        scatter(coord(i,1),coord(i,2),200,'markeredgecolor','k','markerfacecolor', 'g');
    end

end
set(gca,'xtick',[],'ytick',[],'xcolor','w','ycolor','w')


figure(2)
plot(potential,'o')
xlabel('time')
ylabel('potential function \Phi')
set(gca,'FontSize',15)

%% 2b. Wifi
clear
clc
close all

load -ascii wifi.mat
load -ascii coords.mat

W1 = wifi;
Coord1 = coords;

nbrIter =1000;
phi = [2 1 0 0 0 0 0 0;
       1 2 1 0 0 0 0 0;
       0 1 2 1 0 0 0 0;
       0 0 1 2 1 0 0 0;
       0 0 0 1 2 1 0 0;
       0 0 0 0 1 2 1 0;
       0 0 0 0 0 1 2 1;
       0 0 0 0 0 0 1 2];

nbrNodes = 100;
actions = [1 2 3 4 5 6 7 8]; %actions; 1 = Red, 2 = Green

nbrActions = length(actions);
x = ones(100,1); %initially every node is Red
potential = zeros(1,nbrIter); %potential for each time
 

for k=1:nbrIter % loop during which each node chooses a new action
    %eta = k/100;
    eta = 0.1;
    
    % Selecting a node randomly (with uniform distr.)
    node = randi(nbrNodes);
    
        prob = zeros(nbrActions,1); %probabilities of actions
        for action = actions %loop through all possible actions
            c = 0; %utility function for node if action is "action"
            for j = find(W1(node,:) ~= 0) %all elements in a row ~= 0, => neighbor nodes
                c = c + phi(action,x(j)); % add 1 to utility if neighbor color equals color and add 0 otherwise
            end
            prob(action) = exp(-eta*c); %non-normalized prob
        end
        prob = prob/sum(prob); % normalized elements to probabilities
        F = cumsum(prob); % cumulative probability vector
        newAction = find((F>rand(1)),1); % choose new action according to probabilies
        x(node) = newAction; %will be used to draw new colors in the nodes
    
    
    %Compute total potential
    for i=1:nbrNodes
        for j = find(W1(i,:) ~= 0) %loop through out-neighbors
            potential(k) = potential(k) + (1/2)*phi(x(i),x(j));
        end
    end


end

% Draw nodes ----------
figure(1)
gplot(W1,Coord1,'-k');
hold on
for i = 1:nbrNodes
    if x(i, 1) == 1
        scatter(Coord1(i,1),Coord1(i,2),100,'markeredgecolor','k','markerfacecolor', 'r');
    elseif x(i, 1) == 2 
        scatter(Coord1(i,1),Coord1(i,2),100,'markeredgecolor','k','markerfacecolor', 'g');
    elseif x(i, 1) == 3 
        scatter(Coord1(i,1),Coord1(i,2),100,'markeredgecolor','k','markerfacecolor', 'b');
    elseif x(i, 1) == 4 
        scatter(Coord1(i,1),Coord1(i,2),100,'markeredgecolor','k','markerfacecolor', 'y');
    elseif x(i, 1) == 5 
        scatter(Coord1(i,1),Coord1(i,2),100,'markeredgecolor','k','markerfacecolor', 'm');
    elseif x(i, 1) == 6 
        scatter(Coord1(i,1),Coord1(i,2),100,'markeredgecolor','k','markerfacecolor', 'c');
    elseif x(i, 1) == 7 
        scatter(Coord1(i,1),Coord1(i,2),100,'markeredgecolor','k','markerfacecolor', 'w');
    else 
        scatter(Coord1(i,1),Coord1(i,2),100,'markeredgecolor','k','markerfacecolor', 'k');
    end

end
set(gca,'xtick',[],'ytick',[],'xcolor','w','ycolor','w')
% -----------------

figure(2)

plot(potential,'o')
xlabel('time')
ylabel('potential function \Phi')
title('\eta = t/300')
set(gca,'FontSize',15)