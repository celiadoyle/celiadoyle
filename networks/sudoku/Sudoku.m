%Cecilia Doyle
%hw 2 - sudoku

%% PROBLEM 1

%%--INITIALIZE--%%
%define dimensions
dim1 = 3; %dimension of the micro-board
%get other important dimensions, A matrix, and b vector using function
[A1, b1, numbers1, num_vars1, num_constr1] = initial(dim1);


%%--1A SOLUTION--%
%initial sudoku entries
M_1a = [8,0,0,0,0,0,0,0,0;
    0,0,3,6,0,0,0,0,0;
    0,7,0,0,9,0,2,0,0;
    0,5,0,0,0,7,0,0,0;
    0,0,0,0,4,5,7,0,0;
    0,0,0,1,0,0,0,3,0;
    0,0,1,0,0,0,0,6,8;
    0,0,8,5,0,0,0,1,0;
    0,9,0,0,0,0,4,0,0];
%use functions to find solution
f_1a = objective(num_vars1, dim1, numbers1, M_1a); %create objective function based on initial matrix
x_1a = runopt(num_constr1, num_vars1, A1, b1, f_1a); %run gurobi optimization
M_1a_sol = solution(dim1, numbers1, x_1a); %return final M based on gurobi output


%%--1B SOLUTION--%%
%initial sudoku entries
M_1b = [8,0,0,0,0,0,0,0,3;
    0,0,3,6,0,0,0,0,0;
    0,7,0,0,9,0,2,0,0;
    0,5,0,0,0,7,0,0,0;
    0,0,0,0,4,5,7,0,0;
    9,0,0,1,0,0,6,3,0;
    0,0,1,0,0,0,0,6,8;
    0,0,8,5,0,0,0,1,0;
    5,9,0,0,8,0,4,0,0];
%use functions to find solution
f_1b = objective(num_vars1, dim1, numbers1, M_1b); %create objective function based on initial matrix
x_1b = runopt(num_constr1, num_vars1, A1, b1, f_1b); %run gurobi optimization
M_1b_sol = solution(dim1, numbers1, x_1b); %return final M based on gurobi output


%% PROBLEM 2

%create list of the error ratios
errors = [];

%compare 1000 times
for a = 1:1000
    %generate random M
    M_rand = ceil(9*rand(9,9));

    %use functions to find solution of M_rand
    f_rand = objective(num_vars1, dim1, numbers1, M_rand); %create objective function based on initial matrix
    x_rand = runopt(num_constr1, num_vars1, A1, b1, f_rand); %run gurobi optimization
    M_rand_sol = solution(dim1, numbers1, x_rand); %return final M based on gurobi output

    %check disagreement between M_rand and M_rand_sol
    diff = 0;
    for i = 1:9 %index through rows
        for j = 1:9 %index through columns
            %compare value at each location to count # times the matrices disagree
            if M_rand(i,j) ~= M_rand_sol(i,j)
                diff = diff + 1;
            end
        end
    end

    ratio = diff / 81;
    errors = [errors, ratio];

end

avg = mean(errors);
stdv = std(errors);


%% PROBLEM 3

%%--INITIALIZE--%%
%define dimensions
dim3 = 4; %dimension of the micro-board
%get other important dimensions, A matrix, and b vector using function
[A3, b3, numbers3, num_vars3, num_constr3] = initial(dim3);


%%--3A SOLUTION--%
%initial sudoku entry: empty
M_3a = zeros(16);
%use functions to find solution
f_3a = objective(num_vars3, dim3, numbers3, M_3a); %create objective function based on initial matrix
x_3a = runopt(num_constr3, num_vars3, A3, b3, f_3a); %run gurobi optimization
M_3a_sol = solution(dim3, numbers3, x_3a); %return final M based on gurobi output

%initial sudoku entry: at location 10, 10 the value is 10
M_3b = zeros(16);
M_3b(10,10) = 10;
%use functions to find solution
f_3b = objective(num_vars3, dim3, numbers3, M_3b); %create objective function based on initial matrix
x_3b = runopt(num_constr3, num_vars3, A3, b3, f_3b); %run gurobi optimization
M_3b_sol = solution(dim3, numbers3, x_3b); %return final M based on gurobi output

%% PROBLEM 4


%%--ADJUST FOR NEW CONSTRAINTS--%%
%edit initialized variables, A, and b from problem 1
%start row count at bottom of current A matrix
row = num_constr1;
%add rows to current A matrix
A4 = [A1; zeros(81*3, num_vars1)];
%add rows to b matrix
b4 = [b1; ones(81*3,1)];
%increase number of constraints
num_constr4 = num_constr1 + (81*3);

%new #1
for i = 1:dim1
    for l = 1:dim1
        for m = 1:numbers1
            row = row + 1;
            for j = 1:dim1
                for k = 1:dim1
                    A4(row, (i - 1)*dim1^5 + (j-1)*dim1^4 + (k-1)*dim1^3 + (l-1)*dim1^2 + (m-1) + 1) = 1;
                end
            end
        end
    end
end

%new #2
for j = 1:dim1
    for k = 1:dim1
        for m = 1:numbers1
            row = row + 1;
            for i = 1:dim1
                for l = 1:dim1
                    A4(row, (i - 1)*dim1^5 + (j-1)*dim1^4 + (k-1)*dim1^3 + (l-1)*dim1^2 + (m-1) + 1) = 1;
                end
            end
        end
    end
end

%new #3
for k = 1:dim1
    for l = 1:dim1
        for m = 1:numbers1
            row = row + 1;
            for i = 1:dim1
                for j = 1:dim1
                    A4(row, (i - 1)*dim1^5 + (j-1)*dim1^4 + (k-1)*dim1^3 + (l-1)*dim1^2 + (m-1) + 1) = 1;
                end
            end
        end
    end
end


M4 = [1,2,3,6,4,5,0,0,0;
    4,5,6,0,0,0,0,0,0;
    7,8,9,0,0,0,0,0,0;
    5,0,0,0,0,0,0,0,0;
    8,0,0,0,0,0,0,0,0;
    2,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0,0];

%use functions to find solution
f4 = objective(num_vars1, dim1, numbers1, M4); %create objective function based on initial matrix
x4 = runopt(num_constr4, num_vars1, A4, b4, f4); %run gurobi optimization
M4_sol = solution(dim1, numbers1, x4); %return final M based on gurobi output


%% PRINT ALL RESULTS TOGETHER
disp("Problem 1a")
disp(M_1a_sol)
disp("Problem 1b")
disp(M_1b_sol)
disp(" ")
disp("Problem 2")
disp("Mean = " + avg + " and Standard Deviation = " + stdv)
disp(" ")
disp("Problem 3a")
disp(M_3a_sol)
disp("Problem 3b")
disp(M_3b_sol)
disp(" ")
disp("Problem 4")
disp(M4_sol)

%% 
sense = repmat('=',num_constr1,1);
vtype = repmat('B',num_vars1,1);
model.A = sparse(A1);
model.obj = f_1a;
model.modelsense = 'max';
model.rhs = b1;
model.sense = sense;
model.vtype = vtype;
result = gurobi(model);
result.status
x = result.x;
x=round(x);

%% INITIALIZE
%function to get important dimensions, matrix A, and vector b based on required size of board
function [A, b, numbers, num_vars, num_constr] = initial(dim)

%%--DEFINE DIMENSIONS--%%`
numbers = dim^2; %digit range, in 3x3 1-9, in 4x4 1-16
num_vars = numbers * dim^4; %number of variables
num_constr = 4 * dim^4; %number of constraints
%initialize A matrix and b vector
A = zeros(num_constr,num_vars); %initialize A w 0s, 324 x 729; rows are contraints, columns are variables
b = ones(num_constr, 1); %rhs of constraints, all 1s


%%--DEFINE CONSTRAINTS--%%
%each row of A is a constraint, fill each row with each constraint
row = 0; %to index through rows of A

%family 1: each entry has 1 digit
for i = 1:dim %macro row
    for j = 1:dim %macro col
        for k = 1:dim %micro row
            for l = 1:dim %micro col
                row = row+1;
                for m = 1:numbers
                    A(row, (i - 1)*dim^5 + (j-1)*dim^4 + (k-1)*dim^3 + (l-1)*dim^2 + (m-1) + 1) = 1;
                end
            end
        end
    end
end

%family 2: each micro board has 1-9
for i = 1:dim %macro row
    for j = 1:dim %macro col
        for m = 1:numbers 
            row = row + 1;
            for k = 1:dim %micro row
                for l = 1:dim %micro col
                    A(row, (i - 1)*dim^5 + (j-1)*dim^4 + (k-1)*dim^3 + (l-1)*dim^2 + (m-1) + 1) = 1;
                end
            end
        end
    end
end

%family 3: each column has digit 1-9 only once
for i = 1:dim %macro row
    for k = 1:dim %micro row
        for m = 1:numbers 
            row = row + 1;
            for j = 1:dim %macro col
                for l = 1:dim %micro col
                    A(row, (i - 1)*dim^5 + (j-1)*dim^4 + (k-1)*dim^3 + (l-1)*dim^2 + (m-1) + 1) = 1;
                end
            end
        end
    end
end

%family 4: each row has digit 1-9 only once
for j = 1:dim %macro col
    for l = 1:dim %micro col
        for m = 1:numbers 
            row = row + 1;
            for i = 1:dim %macro row
                for k = 1:dim %micro row
                    A(row, (i - 1)*dim^5 + (j-1)*dim^4 + (k-1)*dim^3 + (l-1)*dim^2 + (m-1) + 1) = 1;
                end
            end
        end
    end
end

end
%% DEFINE OBJECTIVE FUNCTION
function [f] = objective(num_vars, dim, numbers, M)

%define objective function coefficients
f = zeros(num_vars,1); %initialize all coef as 0
M_row = 0; %row index of M
M_col = 0; %column index of M
%check each entry, assign coefficient of 1 to entries that are initialized
%only the initialized variables matter in objective fcn bc want to minimize change
for i = 1:dim %macro row
    for k = 1:dim %micro row
        M_row = M_row + 1; %count which row it is on for position in M
        M_col = 0; %reset column count each time index to new row
        for j = 1:dim %macro col
            for l = 1:dim %micro col
                M_col = M_col + 1; %count which column it is on for position in M
                 for m = 1:numbers
                     if m == M(M_row, M_col)
                         %set coef to 1 if the variable is used in initial 
                         f((i - 1)*dim^5 + (j-1)*dim^4 + (k-1)*dim^3 + (l-1)*dim^2 + (m-1) + 1) = 1;
                     end
                 end
            end
        end
    end
end

end


%% OPTIMIZATION FUNCTION
% use gurobi to compute optimal solution

function [x] = runopt(num_constr, num_vars, A, b, f)

%run optimization
addpath('/Library/gurobi1000/macos_universal2/matlab')
diary solution

sense = repmat('=',num_constr,1);
vtype = repmat('B',num_vars,1);
model.A = sparse(A);
model.obj = f;
model.modelsense = 'max';
model.rhs = b;
model.sense = sense;
model.vtype = vtype;
result = gurobi(model);
result.status
x = result.x;
x=round(x);

diary off

end

%% GET SOLUTION
function[M_sol] = solution(dim, numbers, x)

M_row = 0; %row index of m
M_col = 0; %column index of m

M_sol = zeros(dim^2); %solution matrix M to be filled with final values

%index through each position variable and add the correct value in each entry
for i = 1:dim %macro row
    for k = 1:dim %micro row
        M_row = M_row + 1; %count which row it is on to get position in M
        M_col = 0; %reset column count each time it indexes to a new row
        for j = 1:dim %macro col
            for l = 1:dim %micro col
                M_col = M_col + 1; %count which column it is on to get position in M
                 for m = 1:numbers %check number 1-9
                     %when x_ijklm = 1 the current entry is equal to m (the current # 1-9)
                     if x((i - 1)*dim^5 + (j-1)*dim^4 + (k-1)*dim^3 + (l-1)*dim^2 + (m-1) + 1) == 1
                         M_sol(M_row, M_col) = m; %adds the value to the current location in the matrix
                     end
                 end
            end
        end
    end
end

end