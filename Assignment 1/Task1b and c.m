%% b)Sweden
close all
clear all
load('IOdownload.mat')

Wswe = io.swe2000; % By looking at this matrix we see that Wswe is not strongly connected since there are 0's in entire rows and cols
[N,M] = size(Wswe);

Wswe_new = Wswe(~all(~Wswe),:); % Removing all the rows with only 0's
Wswe_new = Wswe_new(:,~all(~Wswe_new)); % Removing all the cols with only 0's
name_new = name(~all(~Wswe));
% Wswe_new is now the largest connected component

% Getting the eigenvalues and eigenvectors of "Wswe_new'"
[V,D] = eig(Wswe_new');
[N_new,M_new] = size(Wswe_new);

eigenvalues = diag(D); % Creates a vector with all the eigenvalues
norm_of_eigenvalues = abs(eigenvalues); % Creates a vector with the absolute eigenvalues

lambda = maxk(norm_of_eigenvalues,1); % lambda_W

k = size(lambda);
eigenvector = V(:,lambda(1,2)); % My corresponding eigenvector to lambda_W
Y_eig = maxk(eigenvector,3); % Chosing the top 3 values in the eigenvector

final_name_swe = [name_new(Y_eig(1,2)), name_new(Y_eig(2,2)), name_new(Y_eig(3,2))]; % Vector with the names

disp('Eigenvector centrality Sweden 2000: ')
disp(final_name_swe)

%% c) Sweden

Beta = 0.15;

% my is a vector of ones
my = ones(N,1); 

% my is a vector of zeros except for '31 Wholesale & retail trade; repairs'
%my = zeros(N,1); 
%my(31) = 1; 

% Calculating the Katz centrality z and looking at the top 3 sectors
z = inv(diag(ones(N,1))-1/lambda(1,1)*(1-Beta)*Wswe')*Beta*my;
Indecies = maxk(z,3);

final_name_c = [name(Indecies(1,2)), name(Indecies(2,2)), name(Indecies(3,2))];

disp('Katz centrality Sweden 2000: ');
disp(final_name_c);

%% b) Indonesia
close all
clear all
load('IOdownload.mat')

Widn = io.idn2000; % No neeed to remove any rows or cols here since it's already connected
[N,M] = size(Widn);

% Getting the eigenvalues and eigenvectors of "Widn'"
[V,D] = eig(Widn');

eigenvalues = diag(D); % Creates a vector with all the eigenvalues
norm_of_eigenvalues = abs(eigenvalues); % Creates a vector with the absolute eigenvalues

lambda = maxk(norm_of_eigenvalues,1); % lambda_W
eigenvector = V(:,lambda(1,2)); % My corresponding eigenvector to lambda_W
Y_eig = maxk(eigenvector,3); % Chosing the top 3 values in the eigenvector

final_name_idn = [name(Y_eig(1,2)), name(Y_eig(2,2)), name(Y_eig(3,2))];

disp('Eigenvector centrality Indonesia 2000: ')
disp(final_name_idn)

%% c) Indonesia

Beta = 0.15;

% my is a vector of ones
my = ones(N,1); 

% my is a vector of zeros except for '31 Wholesale & retail trade; repairs'
%my = zeros(N,1); 
%my(31) = 1; 

% Calculating the Katz centrality z and looking at the top 3 sectors
z = inv(diag(ones(N,1))-1/lambda(1,1)*(1-Beta)*Widn')*Beta*my;
Indecies = maxk(z,3);

final_name_c = [name(Indecies(1,2)), name(Indecies(2,2)), name(Indecies(3,2))];

disp('Katz centrality Indonesia 2000: ');
disp(final_name_c);