%% a) In- and out-degree calculation
close all
clear all
load('IOdownload.mat')

Wswe = io.swe2000;
Widn = io.idn2000;
[N,M] = size(Wswe);

% Sweden
in_w_swe = sum(Wswe); % In-degree by summing the columns
out_w_swe = sum(Wswe'); % Out-degree by summing the rows

In_max_swe = maxk(in_w_swe,3); % Making a vector with the top three biggest in-degree nodes
Out_max_swe = maxk(out_w_swe,3); % Making a vector with the top three biggest out-degree nodes

% The names of the corresponding top 3 nodes
final_name_in_swe = [name(In_max_swe(1,2)), name(In_max_swe(2,2)), name(In_max_swe(3,2))]; 
final_name_out_swe = [name(Out_max_swe(1,2)), name(Out_max_swe(2,2)), name(Out_max_swe(3,2))];

disp('In-degree centrality Sweden 2000: ')
disp(final_name_in_swe)

disp('Out-degree centrality Sweden 2000: ')
disp(final_name_out_swe)
disp('---------------------------------------------------------------------')

% Indonesia

in_w_idn = sum(Widn); % In-degree by summing the columns
out_w_idn = sum(Widn'); % Out-degree by summing the rows

In_max_idn = maxk(in_w_idn,3);
Out_max_idn = maxk(out_w_idn,3);

final_name_in_idn = [name(In_max_idn(1,2)), name(In_max_idn(2,2)), name(In_max_idn(3,2))];
final_name_out_idn = [name(Out_max_idn(1,2)), name(Out_max_idn(2,2)), name(Out_max_idn(3,2))];

disp('In-degree centrality Indonesia 2000: ')
disp(final_name_in_idn)

disp('Out-degree centrality Indonesia 2000: ')
disp(final_name_out_idn)