clear
close

% This is it desired output
I0 = [1, 1, 3, 5, 9, 17, 32, 32, 17, 5, 2, 1, 0, 0, 0, 0]';

% ------------- Setting the fixed parameters ------------
weeks       = 15;     % Nbr of weeks 
N           = 10;    % Nbr of iterations
n           = 932;    % Number of nodes/people
delta_k     = 1;
delta_beta  = 0.025;
delta_rho   = 0.025;
% ----------------------

    
% ----- These parameter will update -----
k0       = 5;      % Avg degree 
beta0    = 0;    % Probability of infection
rho0     = 1;    % Probability of recovery
% ----------------------

prev_RMSE = 0;

for itr=1:20
    
    % Creating a 3D array with the different configurations of the params.
    RMSE = zeros(3,3);
    RMSE(:,:,2) = zeros(3,3);
    RMSE(:,:,3) = zeros(3,3);
    
    idx_k = 1;
    for k_test = [k0 - delta_k, k0, k0 + delta_k]
        
        idx_beta = 1;
        for beta_test = [beta0 - delta_beta, beta0, beta0 + delta_beta]
            
            idx_rho = 1;
            for rho_test = [rho0 - delta_rho, rho0, rho0 + delta_rho]
                
                
                
                % The weight matrix, generated
                W = random_graph_generator(n,k_test);
                
                % The different possible states
                S = 0; % Susceptible
                I = 1; % Infected
                R = 2; % Recovered
                V = 3; % Vaccinated
                
                % Vaccination vector
                Vacc = [0, 5, 15, 25, 35, 45, 55, 60, 60, 60, 60, 60, 60, 60, 60].*1e-2;
                
                % Will keep track of how many gets infected each week
                newly_infected   = zeros(weeks,N);
                
                % This will keep track on how many people there are in each state for every
                % week.
                nbr_susceptible         = zeros(weeks,N);
                nbr_infected            = zeros(weeks,N);
                nbr_recovered           = zeros(weeks,N);
                nbr_vaccinated          = zeros(weeks,N);
                
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
                        not_vaccinated = find(~(X(:,week)==V));
                        vaccinated = find((X(:,week)==V));
                        
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
                            inf_prob = 1-(1-beta_test)^m;
                            
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
                            if rand<=rho_test
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
                        
                        % This is the amount of people I need to vaccinate
                        nbr_to_vacc = int16((Vacc(week+1)-Vacc(week))*n);
                        
                        % Return the indicies of the people we should vaccinate
                        ppl_to_vacc = not_vaccinated(randperm(length(not_vaccinated),nbr_to_vacc));
                        
                        for o=1:nbr_to_vacc
                            X(ppl_to_vacc(o),week+1) = V;
                        end
                        
                        % If a node is vaccinated it will keep being vaccinated next week
                        for q=1:length(vaccinated)
                            X(vaccinated(q),week+1) = V;
                        end
                        
                        nbr_susceptible(week+1,iter) = length(find((X(:,week+1)==S)));
                        nbr_infected(week+1,iter)    = length(find((X(:,week+1)==I)));
                        nbr_recovered(week+1,iter)   = length(find((X(:,week+1)==R)));
                        nbr_vaccinated(week+1,iter)  = length(find((X(:,week+1)==V)));
                    end
                    
                end
                
                avg_susc = [n;mean(nbr_susceptible,2)];
                avg_inf  = [0;mean(nbr_infected,2)];
                avg_rec  = [0;mean(nbr_recovered,2)];
                avg_vac  = [0;mean(nbr_vaccinated,2)];
                
                avg_newly_inf = [0;mean(newly_infected,2)]; % This one is of interest
                I = avg_newly_inf;
                
                RMSE(idx_rho,idx_beta,idx_k) = sqrt((1/15)*sum((I-I0).^2));
                idx_rho = idx_rho + 1;
            end
            idx_beta = idx_beta + 1;
        end
        idx_k = idx_k + 1;
    end
    
    % Looking at what combination of params. create the smallest error
    [v,loc] = min(RMSE(:));
    [ii,jj,ll] = ind2sub(size(RMSE),loc);
    
    
    
    disp(['Iteration: ', int2str(itr)])
    disp(['The smallest RMSE was: ', int2str(v)])
    
    
    if ii==2 & jj==2 & ll==2
        disp(['No update was done to the parameters!'])
        break;
    end
    
    % This is used for the updating of the parameters
    k_temp      = [k0 - delta_k, k0, k0 + delta_k];
    beta_temp   = [beta0 - delta_beta, beta0, beta0 + delta_beta];
    rho_temp    = [rho0 - delta_rho, rho0, rho0 + delta_rho];
    
    % Updating the parameters
    k0       = k_temp(ll);      % Avg degree
    beta0    = beta_temp(jj);   % Probability of infection
    rho0     = rho_temp(ii);    % Probability of recovery
    
    disp(['New k: ', int2str(k0)])
    disp(['New beta: ', int2str(beta0)])
    disp(['New rho: ', int2str(rho0)])
    disp(['-----'])
    
    if (prev_RMSE == v)
        delta_beta = delta_beta/2;
        delta_rho = delta_rho/2;
    end
    prev_RMSE = v;
    
end

disp(['The best k: ', int2str(k0)])
disp(['The best beta: ', int2str(beta0)])
disp(['The best rho: ', int2str(rho0)])

