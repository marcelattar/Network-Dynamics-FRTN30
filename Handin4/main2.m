clear all
clc
close all

% ------------- Setting the parameters (Change if you want to) ----------------
beta    = 0.3;    % Probability of infection
rho     = 0.7;    % Probability of recovery
weeks   = 15;     % Nbr of weeks 
N       = 100;    % Nbr of iterations
n       = 500;    % Number of nodes
k       = 6;      % Avg degree 
% ------------------------------------

% The weight matrix, generated 
W = random_graph_generator(n,k);

% The different possible states
S = 0; % Susceptible 
I = 1; % Infected
R = 2; % Recovered

% Will keep track of how many gets infected each week
newly_infected = zeros(weeks,N);

% This will keep track on how many people there are in each state for every
% week.
nbr_susceptible         = zeros(weeks,N);
nbr_infected            = zeros(weeks,N);
nbr_recovered           = zeros(weeks,N);

for iter=1:N
    % State matrix
    X = zeros(n,weeks);
    
    % Initialize the infection randomly
    start_infected = 10;
    r = randperm(n,start_infected);
    X(r,1) = I;
    clear randperm
    
    % Adding the first week
    newly_infected(1,iter) = start_infected; % First week 'nbr_infected' (10) got infected
    
    % Adding the initial state (week 1)
    nbr_susceptible(1,iter) = n - start_infected;
    nbr_infected(1,iter)    = start_infected;
    
    for week=1:weeks-1
        susceptible = find((X(:,week)==S));
        infected = find((X(:,week)==I));
        recovered = find((X(:,week)==R));
        
        
        % I want to check what nodes will get infected
        for l=1:length(susceptible)
            clear rand
            % Finding the neighbors of the susceptible node
            neighbors = find(W(susceptible(l),:));
            
            % Checking if the neighbors are infected
            infected_neighbors = neighbors((X(neighbors,week)==I));
            
            % Nbr of infected neighbors
            m = length(infected_neighbors);
            
            % Infection probability of node
            inf_prob = 1-(1-beta)^m;
            
            % Node is infected
            if rand<=inf_prob
                X(susceptible(l),week+1) = I;
                
                % Update how many have been infected this week
                newly_infected(week+1,iter) = newly_infected(week+1,iter) + 1;
            end
        end
        
        % Same as above but for infected
        for j=1:length(infected)
            clear rand
            % Node is recovered
            if rand<=rho
                X(infected(j),week+1) = R;
                % Node is still infected
            else
                X(infected(j),week+1) = I;
            end
        end
        
        % If a node is recovered it will keep being recovered next week
        for p=1:length(recovered)
            X(recovered(p),week+1) = R;
        end
        
        nbr_susceptible(week+1,iter) = length(find((X(:,week+1)==S)));
        nbr_infected(week+1,iter)    = length(find((X(:,week+1)==I)));
        nbr_recovered(week+1,iter)   = length(find((X(:,week+1)==R)));
    end
    
end

avg_susc = [n;mean(nbr_susceptible,2)];
avg_inf  = [0;mean(nbr_infected,2)];
avg_rec  = [0;mean(nbr_recovered,2)];
avg_newly_inf = [0;mean(newly_infected,2)];

% Plotting
FigH1 = figure('Position', get(0, 'Screensize'));
graph1 = plot(linspace(0,weeks,weeks+1),avg_susc,...,
              linspace(0,weeks,weeks+1),avg_inf,...,
              linspace(0,weeks,weeks+1),avg_rec);
xlabel('time [week]')
ylabel('Number of people')
title('Epidemic spread over time')
set(gca,'FontSize',15);
set(graph1,'LineWidth',2);
lgd = legend('Susceptible','Infected','Recovered');
lgd.FontSize = 15;

FigH2 = figure('Position', get(0, 'Screensize'));
graph2 = plot(linspace(0,weeks,weeks+1),avg_newly_inf);
xlabel('time [week]')
ylabel('Number of newly infected people')
title('Newly infected over time')
set(gca,'FontSize',15)
set(graph2,'LineWidth',2);

% Saving the figures
path = '/figures_png/';
filename1 = [path, 'Task2_fig1.png'];
print(FigH1,[pwd filename1],'-dpng','-r100');
filename2 = [path, 'Task2_fig2.png'];
print(FigH2,[pwd filename2],'-dpng','-r100');