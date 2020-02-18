clear all
load -ascii flow.mat
load -ascii traffic.mat
load -ascii capacities.mat
load -ascii traveltime.mat

% Declaring the initial variables
B = traffic; % node-link incidence matrix
c = capacities;
l = traveltime; % the minimum travel time 
[M,N] = size(B);

% Here I create vectors that will note the different links between nodes.
out_vec = zeros(1,M); 
in_vec = zeros(1,M);
for col=1:N
   for row=1:M
       if B(row,col)==1
           out_vec(col) = row;
       elseif B(row,col)==-1
           in_vec(col) = row;
       end
   end
end
% This is my direct-graph matrix, I will be able to use this in the
% function "graphshortestpath".
DG = sparse(out_vec,in_vec,l,17,17);
C = [0 c(1) 0 0 0 c(5) 0 0 0 0 0 0 0 0 0 0 0;
    0 0 c(2) 0 0 0 c(10) 0 0 0 0 0 0 0 0 0 0;
    0 0 0 c(3) 0 0 0 c(11) c(12) 0 0 0 0 0 0 0 0;
    0 0 0 0 c(4) 0 0 0 c(13) 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 c(14) 0 0 0;
    0 0 0 0 0 0 c(6) 0 0 c(15) 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 c(7) 0 c(18) 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 c(8) 0 c(19) 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 c(20) c(9) 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 c(16) 0 0 0 c(17) 0 0;
    0 0 0 0 0 0 0 0 0 0 0 c(21) 0 0 c(24) 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 c(22) 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 c(23) 0 0 c(25);
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 c(26);
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 c(27) 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 c(28);
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
C = sparse(C);



% Assigning my different nodes.
node1 = 1; % Start node
node2 = 17; % End node

% a) Here I calculate the distance (time) and the shortest path.
[shortest_time,path,~] = graphshortestpath(DG,node1,node2); 

% b) Calculating the maximum flow between node 1 and node 17.
[Max_flow,~,~] = graphmaxflow(C,node1,node2);

% c) compute the external inflow or outflow at each node.
out_in_flow = B*flow; % '+' = out flow, '-' = in flow

% d) Find the social optimum flow & f)
lambda = zeros(M,1); % External in-flow
lambda(1) = out_in_flow(1);

mu = zeros(M,1); % External out-flow
mu(17) = out_in_flow(1);

% Calculating the social optimum flow, f, for each link
cvx_begin
    variable f(N);
    minimize sum(l.*c.*inv_pos(1-f.*inv_pos(c))-l.*c); % I use in_pos instead of './'
    subject to
        B*f == lambda - mu;
        0 <= f <= c;
cvx_end

%% e) find the Wardrop equilibrium 

% Calculating the Wardrop equilibrium flow, f0, for each link
cvx_begin
    variable f0(N);
    minimize sum(-c.*l.*log(1-f0.*inv_pos(c))); % I use in_pos instead of './'
    subject to
        B*f0 == lambda - mu;
        0 <= f0 <= c;
cvx_end

%% f) Introducing tolls
fd = f;

w = fd.*((traveltime.*capacities)./((capacities-fd).*(capacities-fd)));

cvx_begin
    variable f_f(size(B,2))
    minimize sum(-traveltime.*capacities.*log(capacities-f_f)+traveltime.*capacities.*log(capacities) + f_f.*w)
    subject to
        B*f_f == lambda-mu
        0 <= f_f <= capacities
cvx_end

round(w);
round(f_f);
fdiff=round(f_f-fd);

%% g
clear all
load -ascii traffic.mat
load -ascii traveltime.mat
load -ascii capacities.mat
B=traffic; clear traffic
lambda=zeros(17,1); 
lambda(1)=16806;
mu=zeros(17,1); 
mu(17)=16806;

cvx_begin
    variable f(size(B,2));
    minimize sum(inv_pos((ones(size(B,2),1)-f./capacities)./(traveltime.*capacities))-traveltime.*capacities-traveltime.*f)
    subject to
        B*f == lambda-mu;
        0 <= f <= capacities
cvx_end

fsocial=f;
w = fsocial.*((traveltime.*capacities)./((capacities-fsocial).*(capacities-fsocial)))-traveltime;

cvx_begin
    variable f(size(B,2))
    minimize sum(-traveltime.*capacities.*log(capacities-f)+traveltime.*capacities.*log(capacities) + f.*w)
    subject to
        B*f == lambda-mu;
        0 <= f <= capacities
cvx_end

fdiff=round(f-fsocial);

