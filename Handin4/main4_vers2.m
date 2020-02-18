%% 4. The H1N1 pandemic in Sweden 2009
clear all
clc

n = 934; %nbr of nodes
S = 0;
I = 1;
R = 2;
V = 3;
weeks = 16;

Vacc = [5; 9; 16; 24; 32; 40; 47; 54; 59; 60; 60; 60; 60; 60; 60; 60];
I0 = [1; 1; 3; 5; 9; 17; 32; 32; 17; 5; 2; 1; 0; 0; 0; 0];

kStart = 9; %Initial testing value
betaStart = 0.3; %Initial testing value
pStart = 0.5; %Initial testing value
deltaK=1;
deltaBeta=0.1;
deltaP=0.1;

RMSparameters=zeros(3,1);

while true
    
    RMSold = 100000;
    
    for TestForK=1:3
        if TestForK==1
            k = kStart-deltaK;
        elseif TestForK==2
            k = kStart;
        else
            k = kStart+deltaK;
        end
        
        k0=k+1; % Initial number of nodes in current graph, k0 > c
        W = ones(k0,k0)-diag(ones(k0,1)); % Adjacency matrix, initial
        W = sparse(W); % Transform W into a sparse matrix
        
        for i=(k0+1):n
            if mod(i,2)==0 % If i even number
                c = floor(k/2); % Degree of new node, rounds up to lowest
            else % If i odd number
                c = ceil(k/2); % Degree of new node, rounds up to highest
            end
            w = sum(W,2);
            P = w./sum(w);
            for j=1:c % Choose c neighbours
                % Choose one neighbour with prob. prop. to node degrees
                neighbour = randsample(1:k0,1,true,full(P));
                % Assure not to choose the same neighbour again
                P(neighbour) = 0;
                W(k0+1,neighbour) = 1;
                W(neighbour,k0+1) = 1;
            end
            k0 = k0 + 1; % Number of nodes in current graph
        end
        
        wout = sum(W,2);
        [Val,TestCentralNodes]= maxk(wout,10); %Finding the nodes with highest outdegree
        [Val2,TestLeastCentralNodes]= mink(wout,10); %Finding the nodes with lowest outdegree
        
        for TestForBeta = 1:3
            if TestForBeta == 1
                beta = betaStart-deltaBeta;
            elseif TestForBeta == 2
                beta = betaStart;
            else
                beta = betaStart+deltaBeta;
            end
            
            for TestForP = 1:3
                if TestForP == 1
                    p = pStart-deltaP;
                elseif TestForP == 2
                    p = pStart;
                else
                    p = pStart+deltaP;
                end
                
                
                
                NewlyInfec = zeros(weeks,1);
                areInfected = zeros(weeks,1);
                areSusceptible = zeros(weeks,1);
                areRecovered = zeros(weeks,1);
                GotVacc = zeros(weeks,1);
                areVaccinated=zeros(weeks,1);
                
                niter = 10;
                
                for iter = 1:niter
                    
                    X = zeros(n,1);
                    
                    Infected0 = randperm(n,1); %Test random nodes
                    %Infected0 = TestLeastCentralNodes'; %test for least central nodes
                    %Infected0 = TestCentralNodes'; %test for most central nodes
                    X(Infected0) = I; % 1 randomly chosen node to initially infect
                    Vaccinated0 = randperm(n,round(Vacc(1)*n/100));
                    X(Vaccinated0)=V;
                    
                    %saving initial values from the first week
                    areInfected(1) = areInfected(1) + length(find(X==I)); %initial number of infected nodes
                    areSusceptible(1) = areSusceptible(1) + length(find(X==S));
                    areRecovered(1) = areRecovered(1) + length(find(X==R));
                    areVaccinated(1) = areVaccinated(1) + length(find(X==V));
                    
                    [nodes,value] = find(X~=V);
                    GotVacc(1) = round((Vacc(2)-Vacc(1))*n/100); %Number of people to vaccinate under a week
                    NodesToVacc = randperm(length(find(X~=V)),GotVacc(1)); %nodes vaccinated
                    NodesToVacc = nodes(NodesToVacc);
                    X(NodesToVacc) = V;
                    
                    for week = 2:weeks
                        susceptible = find(X==S); %find nodes that are susceptible
                        infected = (X==I);
                        StoI = 0;
                        ItoR = 0;
                        counter1 = 1;
                        counter2 = 1;
                        for index = 1:length(susceptible) %Test if susceptible nodes get infected
                            neighbours = find(W(susceptible(index),:));
                            m = sum(infected(neighbours)); %number of neighbor nodes that are infected
                            prob = 1-(1-beta)^m;
                            clear rand;
                            rand = rand();
                            if(rand <= prob)
                                StoI(counter1) = susceptible(index)';
                                counter1 = counter1 + 1;
                            end
                        end
                        infected =  find(infected);
                        for index = 1:length(infected) %Test if infected nodes get recovered
                            clear rand;
                            rand = rand();
                            if(rand <= pStart)
                                ItoR(counter2) = infected(index);
                                counter2 = counter2 + 1;
                            end
                        end
                        if ItoR~=0
                            X(ItoR)=R; %Store which nodes that recovers
                        end
                        if StoI~=0
                            X(StoI)=I; %Store which nodes that gets infected
                        end
                        NewlyInfec(week-1) = NewlyInfec(week-1) + length(find(StoI));
                        areInfected(week) = areInfected(week) + length(find(X==I));
                        areSusceptible(week) = areSusceptible(week) + length(find(X==S));
                        areRecovered(week) = areRecovered(week) + length(find(X==R));
                        areVaccinated(week) = areVaccinated(week) + length(find(X==V));
                        if week~=weeks
                            [nodes,value] = find(X~=V);
                            GotVacc(week) = round((Vacc(week+1)-Vacc(week))*n/100); %Number of people to vaccinate under a week
                            NodesToVacc = randperm(length(find(X~=V)),GotVacc(week)); %nodes vaccinated
                            NodesToVacc = nodes(NodesToVacc);
                            X(NodesToVacc) = V;
                        end
                    end %end week
                end
                
                
                NewlyInfec = NewlyInfec/niter;
                areInfected = areInfected/niter;
                areSusceptible = areSusceptible/niter;
                areRecovered = areRecovered/niter;
                areVaccinated = areVaccinated/niter;
                
                
                RMS = norm(NewlyInfec-I0)*sqrt(1/16)
                if RMS<RMSold
                    RMSold=RMS;
                    RMSparameters(1)=TestForK; %store best k value
                    RMSparameters(2)=TestForBeta; %store best Beta value
                    RMSparameters(3)=TestForP; %store best p value
                    NewlyInfBest = NewlyInfec;
                    BestInfected = areInfected;
                    BestSusceptible = areSusceptible;
                    BestRecovered = areRecovered;
                    BestVaccinated = areVaccinated;
                end
            end
        end
    end
    
    if RMSparameters ==2
        disp('Best estimate has been found')
        
        return
    else
        if RMSparameters(1) ~= 2
            if RMSparameters(1) == 1
                kStart = kStart-deltaK;
            else
                kStart = kStart+deltaK;
            end
        end
        if RMSparameters(2) ~= 2
            if RMSparameters(2) == 1
                betaStart = betaStart-deltaBeta;
            else
                betaStart = betaStart+deltaBeta;
            end
        end
        if RMSparameters(3) ~= 2
            if RMSparameters(3) == 1
                pStart = pStart-deltaP;
            else
                pStart = pStart+deltaP;
            end
        end
    end
    
    %plot
    
end
RMSparameters 
     %%   
FigH1 = figure('Position', get(0, 'Screensize'));
hold on
title('Newly infected for the simulation vs real')
xlabel('week')
ylabel('Number of people')
xlim([1 16])
plot(NewlyInfBest,'r','Linewidth',2)
plot(I0,'g', 'LineWidth', 2)
lgd = legend('simulation','Real');
set(gca,'FontSize',15);
lgd.FontSize = 15;

FigH2 = figure('Position', get(0, 'Screensize'));
hold on
title('Health of 500 people over 15 weeks')
xlabel('Week') 
ylabel('Number of people')
xlim([1 16])
plot(BestInfected,'r', 'LineWidth', 2)
plot(BestSusceptible,'g', 'LineWidth', 2)
plot(BestRecovered,'b', 'LineWidth', 2)
plot(BestVaccinated,'y', 'LineWidth', 2)
lgd = legend('Infected','Susceptible','Recovered','Vaccinated');
set(gca,'FontSize',15);
lgd.FontSize = 15;
%%
% Saving the figures
path = '/figures_png/';
filename1 = [path, 'Task4_fig1.png'];
print(FigH1,[pwd filename1],'-dpng','-r100');
filename2 = [path, 'Task4_fig2.png'];
print(FigH2,[pwd filename2],'-dpng','-r100');