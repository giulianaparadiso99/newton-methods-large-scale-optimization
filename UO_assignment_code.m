clear
close all
clc
%% ROSENBROCK FUNCTION IN 2D
% parameters
c1 = 1e-4; % used for the Armijo condition
rho = 0.8; %it works better with 0.8 for Modified Newton and Truncated Newton without preconditioning
btmax = 10; %max numer of iterations in the backtracking cycle
gradftol = 1e-8; % tolerance on gradient norm (stopping criterion)
kmax = 5e3; %max number of outer iterations

% starting points
x0=[1.2;1.2];
%x0=[-1.2;1];
%x0=[0;0];
%x0=[-1;-1];
%x0=[0.4;-0.4];

% function handles for f(x), its gradient and its Hessian matrix
fR = @(x) 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;
gradfR = @(x) [-400*x(1)*(x(2)-x(1)^2)-2*(1-x(1)); 200*(x(2)-x(1)^2)];
HfR= @(x) [-400*x(2)+1200*x(1)^2+2, -400*x(1) ; -400*x(1), 200];

%% MODIFIED NEWTON METHOD
disp ('°°°°°°°° START OF MODIFIED NEWTON METHOD FOR ROSENBROCK FUNCTION °°°°°°°°')
tic ;
[xk, fk, k, gradfk_norm, xseq, btseq] =...
     modifiednewton(x0, fR, gradfR, HfR, gradftol, kmax, c1, rho, btmax);
elapsed_time = toc ;
disp ('°°°°°°°° FINISHED °°°°°°°°')
disp ('°°°°°°°° MODIFIED NEWTON METHOD RESULTS °°°°°°°°')
disp ([ 'norm ( gradf(x)):', num2str(gradfk_norm) ])
disp(['N. of iterations: ', num2str(k), '/', num2str(kmax), ';'])
disp(['Elapsed Time: ', num2str(elapsed_time), ' seconds'])

% experimental rate of convergence given the minimum point
x_min=[1;1];
[pmedian, pseq] = exp_conv_rate_min(x_min, k, xseq);

disp ('°°°°°°°°°°°°°°°°°°°°°°°°°')

%% TRUNCATED NEWTON METHOD WITHOUT PRECONDITIONING
disp ('°°°°°°°° START OF TRUNCATED NEWTON METHOD WITHOUT PRECONDITIONING FOR ROSENBROCK FUNCTION °°°°°°°°')
tic ;
[xk , fk , gradfk_norm , k, xseq , btseq ] = ...
     truncated_newton(x0 , fR, gradfR , HfR , kmax , gradftol , c1 , rho , btmax);
elapsed_time = toc ;
disp ('°°°°°°°° FINISHED °°°°°°°°')
disp ('°°°°°°°° TRUNCATED NEWTON METHOD WITHOUT PRECONDITIONING RESULTS °°°°°°°°')
disp ([ 'norm ( gradf(x)):', num2str(gradfk_norm) ])
disp(['N. of iterations: ', num2str(k), '/', num2str(kmax), ';'])
disp(['Elapsed Time: ', num2str(elapsed_time), ' seconds'])

% experimental rate of convergence given the minimum point
x_min=[1;1];
[pmedian, pseq] = exp_conv_rate_min(x_min, k, xseq);

disp ('°°°°°°°°°°°°°°°°°°°°°°°°°')

%% TRUNCATED NEWTON METHOD WITH PRECONDITIONING
rho = 0.5; % preconditioned Truncated Newton method works better with rho=0.5
disp ('°°°°°°°° START OF TRUNCATED NEWTON METHOD WITH PRECONDITIONING FOR ROSENBROCK FUNCTION °°°°°°°°')
tic ;
[xk , fk , gradfk_norm , k, xseq , btseq ] = ...
    truncated_newton_pre(x0 , fR, gradfR , HfR , kmax , gradftol , c1 , rho , btmax );
elapsed_time = toc ;
disp ('°°°°°°°° FINISHED °°°°°°°°')
disp ('°°°°°°°° TRUNCATED NEWTON METHOD WITH PRECONDITIONING RESULTS °°°°°°°°')
disp ([ 'norm ( gradf(x)):', num2str(gradfk_norm) ])
disp(['N. of iterations: ', num2str(k), '/', num2str(kmax), ';'])
disp(['Elapsed Time: ', num2str(elapsed_time), ' seconds'])

% experimental rate of convergence given the minimum point
x_min=[1;1];
[pmedian, pseq] = exp_conv_rate_min(x_min, k, xseq);

disp ('°°°°°°°°°°°°°°°°°°°°°°°°°')

clear
close all
clc

%% PROBLEM 27: PENALTY FUNCTION
% parameters
n=1e3; %dimension of the R^n space
c1 = 1e-4; % used for the Armijo condition
rho = 0.8; %works better with 0.8 rather than 0.5
btmax = 10; %max numer of iterations in the backtracking cycle
gradftol = 1e-8; % tolerance on gradient norm (stopping criterion)
kmax = 5e3; %max number of outer iterations

% starting points
x0=zeros(n,1);
%x0=ones(n,1);
%x0=1/n*ones(n,1);
%x0=linspace(1,n,n)';

x0 = zeros(n,1);
x0(10:10:end) = -1;

% function handles for f(x), its gradient and its Hessian matrix
f27 = @(x) 1/2*(sum( (1/sqrt(100000)*(x(n:1)-1)).^2 )) ...
    + 1/2*sum( (x(1:n).^2 -1/4) ).^2;
gradf27 = @ (x) (1/ sqrt(100000)*( x(1: end) - 1) + ...
    2*x(1: end).*sum( (x(1: end ).^2 - 1/4) ) ) ; 
Hf27_diag = @ (x) diag ((1/ sqrt (100000) + ...
    2* sum( x (1: end) .^2 - 1/4) ) + 4*( x(1: end ) .^2) ) ;
Hf27_nondiag = @ (x) 4 * ( x * x') - diag (4 * x .^2) ;
Hf27 = @ (x) Hf27_diag ( x ) + Hf27_nondiag ( x ) ;

%% MODIFIED NEWTON METHOD
disp ('°°°°°°°° START OF MODIFIED NEWTON METHOD FOR PROBLEM 27 FUNCTION °°°°°°°°')
tic ;
[xk, fk, k, gradfk_norm, xseq, btseq] =...
    modifiednewton(x0, f27, gradf27, Hf27, gradftol, kmax, c1, rho, btmax);
elapsed_time = toc ;
disp ('°°°°°°°° FINISHED °°°°°°°°')
disp ('°°°°°°°° MODIFIED NEWTON METHOD RESULTS °°°°°°°°')
disp ([ 'norm ( gradf(x)):', num2str(gradfk_norm) ])
disp(['N. of iterations: ', num2str(k), '/', num2str(kmax), ';'])
disp(['Elapsed Time: ', num2str(elapsed_time), ' seconds'])

% experimental rate of convergence without knowing the minimum point
[pmedian, pseq] = exp_conv_rate(k, xseq);

disp ('°°°°°°°°°°°°°°°°°°°°°°°°°')

%% TRUNCATED NEWTON METHOD WITHOUT PRECONDITIONING
disp ('°°°°°°°° START OF TRUNCATED NEWTON METHOD WITHOUT PRECONDITIONING FOR PROBLEM 27 FUNCTION °°°°°°°°')
tic ;
[xk , fk , gradfk_norm , k, xseq , btseq ] = ...
     truncated_newton(x0 , f27, gradf27 , Hf27 , kmax , gradftol , c1 , rho , btmax);
elapsed_time = toc ;
disp ('°°°°°°°° FINISHED °°°°°°°°')
disp ('°°°°°°°° TRUNCATED NEWTON METHOD WITHOUT PRECONDITIONING RESULTS °°°°°°°°')
disp ([ 'norm ( gradf(x)):', num2str(gradfk_norm) ])
disp(['N. of iterations: ', num2str(k), '/', num2str(kmax), ';'])
disp(['Elapsed Time: ', num2str(elapsed_time), ' seconds'])

% experimental rate of convergence without knowing the minimum point
[pmedian, pseq] = exp_conv_rate(k, xseq);

disp ('°°°°°°°°°°°°°°°°°°°°°°°°°')

%% TRUNCATED NEWTON METHOD WITH PRECONDITIONING
disp ('°°°°°°°° START OF TRUNCATED NEWTON METHOD WITH PRECONDITIONING FOR PROBLEM 27 FUNCTION °°°°°°°°')
tic ;
[xk , fk , gradfk_norm , k, xseq , btseq ] = ...
     truncated_newton_pre(x0 , f27, gradf27 , Hf27 , kmax , gradftol , c1 , rho , btmax );
elapsed_time = toc ;
disp ('°°°°°°°° FINISHED °°°°°°°°')
disp ('°°°°°°°° TRUNCATED NEWTON METHOD WITH PRECONDITIONING RESULTS °°°°°°°°')
disp ([ 'norm ( gradf(x)):', num2str(gradfk_norm) ])
disp(['N. of iterations: ', num2str(k), '/', num2str(kmax), ';'])
disp(['Elapsed Time: ', num2str(elapsed_time), ' seconds'])

% experimental rate of convergence without knowing the minimum point
[pmedian, pseq] = exp_conv_rate(k, xseq);

disp ('°°°°°°°°°°°°°°°°°°°°°°°°°')

clear
close all
clc

%% PROBLEM 28: VARIABLY DIMENSIONED FUNCTION
% parameters
n=1e3; %dimension of the R^n space
c1 = 1e-4; % used for the Armijo condition
rho = 0.8; %works better with 0.8 rather than 0.5
btmax = 10; %max numer of iterations in the backtracking cycle
gradftol = 1e-8; % tolerance on gradient norm (stopping criterion)
kmax = 5e3; %max number of outer iterations

% starting points
x0=zeros(n,1); 
%x0=ones(n,1);
%x0=1/n*ones(n,1);
% x0=linspace(1,n,n)';

% function handles for f(x), its gradient and its Hessian matrix
f28 = @(x) 1/2*sum((x(1:n)-1).^2) ...
    + 1/2*(sum( (1:n)*(x(1:n)-1) ).^2) ...
    + 1/2*(sum( (1:n)*(x(1:n)-1) ).^4);
gradf28 = @ (x) ( x(1: end ) + (1: n )'.* sum ((1: n )'.*( x (1: end ) -1) ) + ...
    2*(1: n )'.* sum ((1: n )'.*( x (1: end) -1) ) .^3) ;
Hf28 = @ (x) ((1: n )'.*(1: n ) + ...
    6*(1: n )'.*(1: n ) .*( sum ((1: n )'.*( x (1:end ) -1) ) ) .^2 + diag ( ones (n ,1) ) ) ;

%% MODIFIED NEWTON METHOD
disp ('°°°°°°°° START OF MODIFIED NEWTON METHOD FOR PROBLEM 28 FUNCTION °°°°°°°°')
tic ;
[xk, fk, k, gradfk_norm, xseq, btseq] =...
     modifiednewton(x0, f28, gradf28, Hf28, gradftol, kmax, c1, rho, btmax);
elapsed_time = toc ;
disp ('°°°°°°°° FINISHED °°°°°°°°')
disp ('°°°°°°°° MODIFIED NEWTON METHOD RESULTS °°°°°°°°')
disp ([ 'norm ( gradf(x)):', num2str(gradfk_norm) ])
disp(['N. of iterations: ', num2str(k), '/', num2str(kmax), ';'])
disp(['Elapsed Time: ', num2str(elapsed_time), ' seconds'])

% experimental rate of convergence without knowing the minimum point
[pmedian, pseq] = exp_conv_rate(k, xseq);

disp ('°°°°°°°°°°°°°°°°°°°°°°°°°')

%% TRUNCATED NEWTON METHOD WITHOUT PRECONDITIONING
disp ('°°°°°°°° START OF TRUNCATED NEWTON METHOD WITHOUT PRECONDITIONING FOR PROBLEM 28 FUNCTION °°°°°°°°')
tic ;
[xk , fk , gradfk_norm , k, xseq , btseq ] = ...
     truncated_newton(x0 , f28, gradf28 , Hf28 , kmax , gradftol , c1 , rho , btmax);
elapsed_time = toc ;
disp ('°°°°°°°° FINISHED °°°°°°°°')
disp ('°°°°°°°° TRUNCATED NEWTON METHOD WITHOUT PRECONDITIONING RESULTS °°°°°°°°')
disp ([ 'norm ( gradf(x)):', num2str(gradfk_norm) ])
disp(['N. of iterations: ', num2str(k), '/', num2str(kmax), ';'])
disp(['Elapsed Time: ', num2str(elapsed_time), ' seconds'])

% experimental rate of convergence without knowing the minimum point
[pmedian, pseq] = exp_conv_rate(k, xseq);

disp ('°°°°°°°°°°°°°°°°°°°°°°°°°')

%% TRUNCATED NEWTON METHOD WITH PRECONDITIONING
disp ('°°°°°°°° START OF TRUNCATED NEWTON METHOD WITH PRECONDITIONING FOR PROBLEM 28 FUNCTION °°°°°°°°')
tic ;
[xk , fk , gradfk_norm , k, xseq , btseq ] = ...
     truncated_newton_pre(x0 , f28, gradf28 , Hf28 , kmax , gradftol , c1 , rho , btmax );
elapsed_time = toc ;
disp ('°°°°°°°° FINISHED °°°°°°°°')
disp ('°°°°°°°° TRUNCATED NEWTON METHOD WITH PRECONDITIONING RESULTS °°°°°°°°')
disp ([ 'norm ( gradf(x)):', num2str(gradfk_norm) ])
disp(['N. of iterations: ', num2str(k), '/', num2str(kmax), ';'])
disp(['Elapsed Time: ', num2str(elapsed_time), ' seconds'])

% experimental rate of convergence without knowing the minimum point
[pmedian, pseq] = exp_conv_rate(k, xseq);

disp ('°°°°°°°°°°°°°°°°°°°°°°°°°')

clear
close all
clc

%% PROBLEM 5: Generalized Broyden tridiagonal function

% parameters
n=1e3; %dimension of the R^n space
c1 = 1e-4; % used for the Armijo condition
rho = 0.8; %works better with 0.8 rather than 0.5
btmax = 10; %max numer of iterations in the backtracking cycle
gradftol = 1e-8; % tolerance on gradient norm (stopping criterion)
kmax = 5e3; %max number of outer iterations

% starting points
x0=zeros(n,1);
%x0=;
%x0=1/n*ones(n,1);
%x0=0.5*ones(n,1);
%x0=;


% function handles for f(x), its gradient and its Hessian matrix

f_5 = @(x) (abs((3-2*x(1))*x(1) - x(2)+1)).^(7/3) + ... 
sum(abs((3-2*x(2:n-1)).*x(2:n-1) -x(1:n-2)- x(3:n)+1).^(7/3)) + ...
(abs((3-2*x(n))*x(n)-x(n-1)+1)).^(7/3);

gradf5 = @(x) [ ... 
    
    7/3 * (abs((3 - 2*x(1)) * x(1) - x(2) + 1)).^(4/3) .* sign((3 - 2*x(1)) * x(1) - x(2) + 1) * (3 - 4*x(1)) + ...
    7/3* (abs((3 - 2*x(2)) * x(2) - x(1) - x(3) + 1)).^(4/3) .* sign((3 - 2*x(2)) * x(2) - x(1) - x(3) + 1) * (-1); ...
    
    7/3 * (abs((3 - 2*x(1)) * x(1) - x(2) + 1).^(4/3).*sign((3 - 2*x(1)) * x(1) - x(2) +1)) * (-1) + ...
    7/3 * (abs((3 - 2*x(2)) * x(2) - x(1) - x(3) + 1).^(4/3).* sign((3 - 2*x(2)) * x(2) - x(1) - x(3) + 1)) * (3 - 4*x(2)) + ...
    7/3 * (abs((3 - 2*x(3)) * x(3) - x(2) - x(4) + 1).^(4/3).* sign((3 - 2*x(3)) * x(3) - x(2) - x(4) + 1)) * (-1);
    
    7/3 * (abs((3 - 2*x(2:end-3)) .* x(2:end-3) - x(1:end-4) - x(3:end-2) + 1)).^(4/3) .* sign((3 - 2*x(2:end-3)) .* x(2:end-3) - x(1:end-4) - x(3:end-2) + 1) .* (-1) + ...
    7/3 * (abs((3 - 2*x(3:end-2)) .* x(3:end-2) - x(2:end-3) - x(4:end-1) + 1)).^(4/3) .* sign((3 - 2*x(3:end-2)) .* x(3:end-2) - x(2:end-3) - x(4:end-1) + 1) .* (3 - 4*x(3:end-2)) + ...
    7/3 * (abs((3 - 2*x(4:end-1)) .* x(4:end-1) - x(3:end-2) - x(5:end) + 1)).^(4/3) .* sign((3 - 2*x(4:end-1)) .* x(4:end-1) - x(3:end-2) - x(5:end) + 1) .* (-1); ...
    
    7/3 * (abs((3 - 2*x(end-2)) * x(end-2) - x(end-3) - x(end-1) + 1).^(4/3) .* sign((3 - 2*x(end-2)) * x(end-2) - x(end-3) - x(end-1) +1)) * (-1) + ...
    7/3 * (abs((3 - 2*x(end-1)) * x(end-1) - x(end-2) - x(end) + 1).^(4/3) .* sign((3 - 2*x(end-1)) * x(end-1) - x(end-2) - x(end) + 1)) * (3 - 4*x(end-1)) + ...
    7/3 * (abs((3 - 2*x(end)) * x(end) - x(end-1) + 1).^(4/3) .* sign((3 - 2*x(end)) * x(end) - x(end-1) + 1)) * (-1);
    
    7/3 * (abs((3 - 2*x(end-1)) * x(end-1) - x(end-2) - x(end) + 1)).^(4/3) .* sign((3 - 2*x(end-1)) * x(end-1) - x(end-2) - x(end) + 1) * (-1) + ...
    7/3 * (abs((3 - 2*x(end)) * x(end) - x(end-1) + 1)).^(4/3) .* sign((3 - 2*x(end)) * x(end) - x(end-1) + 1) * (3-4*x(end))];

H=zeros(n,n);

d = @(x) [ ... 
    
    28/9 * abs((3 - 2*x(1)) .* x(1) - x(2) + 1).^(1/3) .* (3-4*x(1)).^2 + ...
    7/3 * abs((3 - 2*x(1)) .* x(1) - x(2) +1).^(4/3) .* sign((3 - 2*x(1)) .* x(1) - x(2) +1) * (-4) + ...
    28/9 * abs((3 - 2 * x(2)) .* x(2) - x(1) - x(3) +1).^(1/3); ...
    
    28/9 * abs((3 - 2*x(1)) * x(1) - x(2) +1).^(1/3) + ...
    28/9 * abs((3 - 2*x(2)) * x(2) - x(1) - x(3) + 1).^(1/3) * (3-4*x(2)).^2 + ...
    7/3 * abs((3 - 2*x(2)) * x(2) - x(1) - x(3) + 1).^(4/3) .* sign((3 - 2*x(2)) * x(2) - x(1) - x(3) + 1) * (-4) + ...
    28/9 * abs((3 - 2 * x(3)) * x(3) - x(2) - x(4) +1).^(1/3);
    
    28/9 * abs((3 - 2*x(2:end-3)) .* x(2:end-3) - x(1:end-4) - x(3:n-2) +1).^(1/3) + ...
    28/9 * abs((3 - 2*x(3:end-2)) .* x(3:end-2) - x(2:end-3) - x(4:end-1) + 1).^(1/3) .* (3-4*x(3:end-2)).^2 + ...
    7/3 * abs((3 - 2*x(3:end-2)) .* x(3:end-2) - x(2:end-3) - x(4:end-1) + 1).^(4/3) .* sign((3 - 2*x(3:end-2)) .* x(3:end-2) - x(2:end-3) - x(4:end-1) + 1) * (-4) + ...
    28/9 * abs((3 - 2 * x(4:end-1)) .* x(4:end-1) - x(3:end-2) - x(5:end) +1).^(1/3);

    28/9 * abs((3 - 2*x(end-2)) .* x(end-2) - x(end-3) - x(end-1) +1).^(1/3) + ...
    28/9 * abs((3 - 2*x(end-1)) .* x(end-1) - x(end-2) - x(end) +1).^(1/3) .* (3 - 4*x(end-1)).^2 + ...
    7/3 * abs((3 - 2*x(end-1)) .* x(end-1) - x(end-2) - x(end) +1).^(4/3) .* sign((3 - 2*x(end-1)) .* x(end-1) - x(end-2) - x(end) +1) * (-4) + ...
    28/9 * abs((3 - 2*x(end)) .* x(end) - x(end-1) +1).^(1/3);
   
    28/9 * abs((3 - 2*x(end-1)) .* x(end-1) - x(end-2) - x(end) +1).^(1/3) + ...
    28/9 * abs((3 - 2*x(end)) .* x(end) - x(end-1) +1).^(1/3) .* (3 - 4*x(end)).^2 + ...
    7/3 * abs((3 - 2*x(end)) .* x(end) - x(end-1)+1).^(4/3) .* sign((3 - 2*x(end)) .* x(end) - x(end-1) +1) * (-4)];

d_up = @(x) [
    
    28/9 * abs((3 - 2*x(1)) .* x(1) - x(2) +1).^(1/3) *(-1) * (3 - 4*x(1)) + ...
    28/9 * abs((3 - 2*x(2)) .* x(2) - x(1) - x(3) +1).^(1/3) *(-1) * (3 - 4*x(2));
    
    28/9 * abs((3 - 2*x(2:n-2)) .* x(2:n-2) - x(1:n-3) - x(3:n-1) +1).^(1/3) *(-1) .* (3 - 4*x(2:n-2)) + ...
    28/9 * abs((3 - 2*x(3:n-1)) .* x(3:n-1) - x(2:n-2) - x(4:n) +1).^(1/3) *(-1) .* (3 - 4*x(3:n-1));
   
    28/9 * abs((3 - 2*x(end-1)) .* x(end-1) - x(end-2) - x(end) +1).^(1/3) *(-1) * (3 - 4*x(end-1)) + ...
    28/9 * abs((3 - 2*x(end)) .* x(end) - x(end-1) +1).^(1/3) *(-1) * (3 - 4*x(end))];

d_up2 = @(x) [28/9 * abs((3 - 2*x(2:n-1)) .* x(2:n-1) - x(1:n-2) - x(3:n) +1).^(1/3)];
    

Hf5 = @(x) H + diag(d(x)) + diag(d_up(x),1) + diag(d_up2(x),2) + diag(d_up(x),-1) + diag(d_up2(x),-2);

%% MODIFIED NEWTON METHOD
disp ('°°°°°°°° START OF MODIFIED NEWTON METHOD FOR PROBLEM 5 FUNCTION °°°°°°°°')
tic ;
[xk, fk, k, gradfk_norm, xseq, btseq] =...
    modifiednewton(x0, f_5, gradf5, Hf5, gradftol, kmax, c1, rho, btmax);
elapsed_time = toc ;
disp ('°°°°°°°° FINISHED °°°°°°°°')
disp ('°°°°°°°° MODIFIED NEWTON METHOD RESULTS °°°°°°°°')
disp ([ 'norm ( gradf(x)):', num2str(gradfk_norm) ])
disp(['N. of iterations: ', num2str(k), '/', num2str(kmax), ';'])
disp(['Elapsed Time: ', num2str(elapsed_time), ' seconds'])

% experimental rate of convergence without knowing the minimum point
[pmedian, pseq] = exp_conv_rate(k, xseq);

disp ('°°°°°°°°°°°°°°°°°°°°°°°°°')

%% TRUNCATED NEWTON METHOD WITHOUT PRECONDITIONING
disp ('°°°°°°°° START OF TRUNCATED NEWTON METHOD WITHOUT PRECONDITIONING FOR PROBLEM 5 FUNCTION °°°°°°°°')
tic ;
[xk , fk , gradfk_norm , k, xseq , btseq ] = ...
     truncated_newton(x0 , f_5, gradf5 , Hf5 , kmax , gradftol , c1 , rho , btmax);
elapsed_time = toc ;
disp ('°°°°°°°° FINISHED °°°°°°°°')
disp ('°°°°°°°° TRUNCATED NEWTON METHOD WITHOUT PRECONDITIONING RESULTS °°°°°°°°')
disp ([ 'norm ( gradf(x)):', num2str(gradfk_norm) ])
disp(['N. of iterations: ', num2str(k), '/', num2str(kmax), ';'])
disp(['Elapsed Time: ', num2str(elapsed_time), ' seconds'])

% experimental rate of convergence without knowing the minimum point
[pmedian, pseq] = exp_conv_rate(k, xseq);

disp ('°°°°°°°°°°°°°°°°°°°°°°°°°')

%% TRUNCATED NEWTON METHOD WITH PRECONDITIONING
disp ('°°°°°°°° START OF TRUNCATED NEWTON METHOD WITH PRECONDITIONING FOR PROBLEM 5 FUNCTION °°°°°°°°')
tic ;
[xk , fk , gradfk_norm , k, xseq , btseq ] = ...
     truncated_newton_pre(x0 , f_5, gradf5 , Hf5 , kmax , gradftol , c1 , rho , btmax );
elapsed_time = toc ;
disp ('°°°°°°°° FINISHED °°°°°°°°')
disp ('°°°°°°°° TRUNCATED NEWTON METHOD WITH PRECONDITIONING RESULTS °°°°°°°°')
disp ([ 'norm ( gradf(x)):', num2str(gradfk_norm) ])
disp(['N. of iterations: ', num2str(k), '/', num2str(kmax), ';'])
disp(['Elapsed Time: ', num2str(elapsed_time), ' seconds'])

% experimental rate of convergence without knowing the minimum point
[pmedian, pseq] = exp_conv_rate(k, xseq);

disp ('°°°°°°°°°°°°°°°°°°°°°°°°°')

%% FUNCTIONS IMPLEMENTED TO DEAL WITH THE ASSIGNMENT

%% 1) NEWTON METHOD 

function [xk, fk, k, normgradfk, xseq, btseq] =...
    modifiednewton(x0, f, gradf, Hf, gradftol, kmax, c1, rho, btmax)
% input:
% x0=starting point of the sequence
% f=function to minimize
% gradf=gradient of f
% Hf=hessian matrix of f
% gradftol=tolerance for gradf under which we assume convergence to min
% kmax=max number of iterations for the sequence xk
% c1=parameter of Armijo condition
% rho=reduction parameter for backtracking
% btmax=max number of iterations for backtracking

% output:
% xk=approximation of argmin f
% fk=approximation of min f
% k=number of iterates for convergence
% normgradfk=norm of gradf(xk) at convergence
% xseq=sequence produced by the method
% btseq=sequence produced by backtracking


% Function handle for the armijo condition%
farmijo = @(fk, alpha, gradfk, pk) ...
      fk + c1 * alpha * (-sum((gradfk*pk).^2, 1));

%initialization
xseq = zeros(length(x0), kmax);
btseq = zeros(1, btmax);

xk = x0;
fk = f(xk);
gradfk = gradf(xk);
gradfk=gradfk';
hessfk = Hf(xk);
normgradfk = norm(gradf(xk));
Ek = eye(size(hessfk));
Bk = hessfk;

k = 0;
while k <= kmax && normgradfk >= gradftol
    %build Bk=Hf(xk)+Ek
    beta = min( norm ( gradfk , 'fro') , norm ( hessfk ,'fro') ) ;
    tau = 0;
    for i = 1:size(hessfk, 1)
        if hessfk(i,i) <= 0
            tau = beta/2;
            break
        end
    end 
    flag = 1;
    while flag ~= 0 
        Bk = hessfk + tau*Ek;
        %check if Bk positive definite by trying the Cholesky
        %decomposition
        [~, flag] = chol(Bk);
        %if flag~=0 %if Bk not positive definite, increase tau
            if tau > 0
                tau = 2*tau;
            else
                tau = beta/2;
            end
        %end
    end
    Bk = hessfk + tau*Ek;
  
    %solve the linear system Bk*pk=-gradf(xk)
    pk = Bk \ (-gradfk');

    alpha = 1;
    xnew = xk + alpha * pk;
    fnew = f(xnew);

    %find alpha with backtracking
    bt = 0;
    while bt < btmax && fnew > farmijo(fk, alpha, gradfk, pk) 
        %while Armijo condition is not satisfied, reduce the value of alpha
        alpha = rho * alpha;
        xnew = xk + alpha * pk;
        fnew = f(xnew);
        bt = bt + 1;
    end

    xk = xnew;
    fk = f(xk);
    gradfk = gradf(xk);
    gradfk=gradfk';
    normgradfk = norm(gradfk);
    hessfk = Hf(xk);
    k = k+1;

    % Store current xk in xseq
    xseq(:, k) = xk;
    % Store bt iterations in btseq
    btseq(k) = bt;
end

% remove the extra columns
xseq = xseq (:, 1:k);
btseq = btseq (1:k);
end

%% 2) TRUNCATED NEWTON METHOD WITHOUT PRECONDITIONING 

function [xk , fk , gradfk_norm , k, xseq , btseq ] = ...
    truncated_newton(x0 , f, gradf , Hessf , kmax , tolgrad , c1 , rho , btmax)

% INPUTS :
% x0 = n- dimensional column vector ;
% f = function handle for the function R^n->R;
% gradf = function handle for the gradient of f;
% Hessf = function handle for the Hessian of f;
% kmax = maximum number of outer iterations ;
% tolgrad = value used as stopping criterion for the norm of the
% gradient ;
% c1 = factor of the Armijo condition that must be a scalar in (0 ,1);
% rho = fixed factor used for reducing alpha0 ;
% btmax = maximum number of steps for updating alpha during the
% backtracking strategy .
%
% OUTPUTS :
% xk = the last x computed by the function ;
% fk = the value f(xk);
% gradfk_norm = value of the norm of gradf (xk)
% k = index of the last iteration performed
% xseq = n-by -k matrix where the columns are the xk computed during the
% iterations
% btseq = 1-by -k vector where elements are the number of backtracking
% iterations at each optimization step .


% function handle for the armijo condition
%farmijo = @(fk , alpha , gradfk , pk)fk + c1 * alpha * gradfk' * pk;  %QUI HO MESSO I PUNTI, AURO NO

% Function handle for the armijo condition%
farmijo = @(fk, alpha, gradfk, pk) ...
      fk + c1 * alpha * (-sum((gradfk*pk).^2, 1));

% initializations
xseq = zeros(length(x0), kmax );
btseq = zeros(1, kmax );

xk = x0;
fk = f(xk);
gradfk = gradf (xk);
gradfk=gradfk';
hessfk = Hessf (xk);

gradfk_norm = norm ( gradfk );

k = 0; % counter

%the second one is the stopping criterion
while k < kmax && gradfk_norm > tolgrad
    z = zeros ( length (x0) ,1);
    etak = min (0.5 , ( norm ( gradfk )));
    pk = cg_curvtrun_newt ( hessfk , -gradfk, z, etak ); % cg + curvature condition with no preconditioning

    % Backtracking with armijo condition
    % start with alpha = 1
    alpha = 1;

    % the new candidate xk
    xnew = xk + alpha * pk;
    %the value of f in the candidate new xk
    fnew = f(xnew);

    bt = 0;

    % Backtracking strategy :
    % 2nd condition holds if the Armijo condition is not satisfied
    while bt < btmax && fnew > farmijo(fk, alpha, gradfk, pk) 

        % reduce alpha
        alpha = rho * alpha ;

        % update xnew and fnew w.r.t. the reduced alpha
        xnew = xk + alpha * pk;
        fnew = f( xnew );

        % increase the counter by one
        bt = bt + 1;

    end

    % update the values
    xk = xnew ;
    fk = fnew ;
    gradfk = gradf(xk);
    gradfk=gradfk';
    gradfk_norm = norm( gradfk );
    hessfk = Hessf(xk);

    % increase the step of the outer cycle by one
    k = k + 1;

    % store current xk in xseq
    xseq (:, k) = xk;
    % store bt iterations in btseq
    btseq (k) = bt;
end

% remove the extra columns
xseq = xseq (:, 1:k);
btseq = btseq (1:k);

end

%% 3) CONJUGATE GRADIENT METHOD WITH NEGATIVE CURVATURE CONDITION 

function xk = cg_curvtrun_newt ( A, b, x0 , tol)

% INPUTS :
% A = matrix of the coefficients of the lienar system
% b = vector of constants
% x0 = n- dimensional column vector 
% tol = tolerance on the norm of the relative residual

% OUTPUT :
% xk = the last x computed by the function ;
%b=b';
% initializations
xk = x0;
b=b';
rk = b - A * xk; % residual
pk = rk;
norm_b = norm (b);
relres = norm (rk) / norm_b ; % norm of the relative residual
kmax = 100;
k = 0;
% compute cg iterations until the relative residual gets smaller
% than the tolerance or you reach the max number of iterations

while relres > tol && k < kmax
    if pk' * A * pk <= 0 %if the negative curvature condition is satisfied
        if k == 0 % if it ’s the first iteration , compute the new xk and stop
            zk = A * pk;
            alphak = (rk' * pk )/ (pk' * zk);
            xk = xk + alphak * pk;
        end
        % stop computing xk and return it
        break
    end
    % continue the cycle if not
    zk = A * pk;
    alphak = (rk' * pk )/ (pk' * zk);
    xk = xk + alphak * pk;
    rk = rk - alphak * zk;
    betak = -(rk' * zk) / (pk' * zk); 
    pk = rk + betak * pk;

    relres = norm (rk) / norm_b ;
    
    k = k + 1;
end
end

%% 4) TRUNCATED NEWTON METHOD WITH PRECONDITIONING 

function [xk , fk , gradfk_norm , k, xseq , btseq ] = ...
    truncated_newton_pre(x0 , f, gradf , Hessf , kmax , tolgrad , c1 , rho , btmax )

% function that performs the modified newton method ,
% implementing the backtracking strategy with preconditioning on the cg
% iterations .

% INPUTS :
% x0 = n- dimensional column vector ;
% f = function handle for the function R^n->R;
% gradf = function handle for the gradient of f;
% Hessf = function handle for the Hessian of f;
% kmax = maximum number of outer iterations ;
% tolgrad = value used as stopping criterion for the norm of the
% gradient ;
% c1 = factor of the Armijo condition that must be a scalar in (0 ,1);
% rho = fixed factor used for reducing alpha0 ;
% btmax = maximum number of steps for updating alpha during the
% backtracking strategy .
%
% OUTPUTS :
% xk = the last x computed by the function ;
% fk = the value f(xk);
% gradfk_norm = value of the norm of gradf (xk)
% k = index of the last iteration performed
% xseq = n-by -k matrix where the columns are the xk computed during the
% iterations
% btseq = 1-by -k vector where elements are the number of backtracking
% iterations at each optimization step .

% function handle for the armijo condition
%farmijo = @(fk , alpha , gradfk , pk) fk + c1 * alpha * gradfk' * pk;

% Function handle for the armijo condition%
farmijo = @(fk, alpha, gradfk, pk) ...
      fk + c1 * alpha * (-sum((gradfk*pk).^2, 1));


% initializations
xseq = zeros( length(x0), kmax );
btseq = zeros(1, kmax );

xk = x0;
fk = f(xk);
gradfk = gradf (xk);
gradfk=gradfk';
hessfk = Hessf (xk);

gradfk_norm = norm ( gradfk );

k = 0; % counter

%the second one is the stopping criterion
while k < kmax && gradfk_norm > tolgrad
    z = zeros( length (x0) ,1);
    etak = min (0.5 , gradfk_norm);
    pk = cg_curvtrun_newt_pre( hessfk , -gradfk, z, etak ); %cg + negative curvature condition with preconditioning
    
    % Backtracking with armijo condition
    % start with alpha = 1
    alpha = 1;
    
    % the new candidate xk
    xnew = xk + alpha * pk;
    
    
    
    %the value of f in the candidate new xk
    fnew = f( xnew );
    
    bt = 0;
    
    % Backtracking strategy :
    % 2nd condition holds if the Armijo condition is not satisfied

    while bt < btmax && fnew > farmijo (fk , alpha , gradfk , pk)

        % reduce alpha
        alpha = rho * alpha ;
        
        % update xnew and fnew w.r.t. the reduced alpha
        xnew = xk + alpha * pk;
        fnew = f( xnew );
        
        % increase the counter by one
        bt = bt + 1;
    end

    % update the values
    xk = xnew ;
    fk = fnew ;
    gradfk = gradf (xk);
    gradfk=gradfk';
    gradfk_norm = norm ( gradf(xk) );
    hessfk = Hessf (xk);
    
    % increase the step of the outer cycle by one
    k = k + 1;
    
    % store current xk in xseq
    xseq (:, k) = xk;
    % store bt iterations in btseq
    btseq (k) = bt;
end

% remove the extra columns
xseq = xseq (:, 1:k);
btseq = btseq (1:k);

end

%% 5) PRECONDITIONED CONJUGATE GRADIENT METHOD WITH NEGATIVE CURVATURE CONDITION 

function xk = cg_curvtrun_newt_pre( A, b, x0 , tol)

% INPUTS :
% A = matrix of the coefficients of the lienar system
% b = vector of constants
% x0 = n- dimensional column vector ;.
% tol = tolerance on the norm of the relative residual


% OUTPUT :
% xk = the last x computed by the function ;


% preconditioning operations
%either use the ichol when possible or diag ( diag (A)) as preconditioner
try
   L = ichol (sparse(A)); 
catch exception
   L = diag ( diag(A));
end
%L = diag ( diag (A));

b=b';
L_i = L \eye( size (L ,1) ); % compute the inverse
A = L_i * A; %new matrix
b = L_i * b; % new vector

% initilizations
xk = x0;
rk = b - A * xk;
pk = rk;
norm_b = norm (b);
relres = norm (rk) / norm_b ; % relative residual
kmax = 1000;
k = 0;

% compute pcg iterations until the relative residual gets smaller
% than the tolerance or you reach the max number of iterations

while relres > tol && k < kmax
    if pk' * A * pk <= 0 %if the negative curvature condition is satisfied
        if k == 0 % if its the first iteration , compute the new xk and stop
            zk = A * pk;
            alphak = (rk' * pk )/ (pk' * zk);
            xk = xk + alphak * pk;
        end
        break
        % stop computing xk and return it
    end
    zk = A * pk;
    alphak = (rk' * pk )/ (pk' * zk);
    xk = xk + alphak * pk;
    rk = rk - alphak * zk;
    betak = -(rk' * zk) / (rk' * zk); 
    pk = rk + betak * pk;
    relres = norm (rk) / norm_b ;
    k = k + 1;
end

end

%% 6) EXPERIMENTAL RATE OF CONVERGENCE WITH MINIMUN POINT COMPUTATION 

function [p,pseq] = exp_conv_rate_min(x_min, k, xseq)
% function to compute the experimental rate of convergence of a method,
% given:
% x_min the minimum of the function
% k = number of iterations of the method until convergence
% xseq = sequence of points produced by the iterative method
% The function returns:
% pseq = the sequence of experimental rates computed for 3 consecutive
%   values of the sequence x(k-1), x(k), x(k+1)
% p = the median of the convergence rate vector
if k<=100
    pseq=zeros(k,1);
    for i=2:length(xseq)-1
        ekm=norm(xseq(:,i-1)-x_min);
        ek=norm(xseq(:,i)-x_min);
        ekp=norm(xseq(:,i+1)-x_min);
        pseq(i-1)= log( ekp/ek ) / log( ek/ekm ) ;
    end
    p=median(pseq)
else
    xseq_new=xseq(:,length(xseq)-100:length(xseq));
    pseq=zeros(100,1);
    for i=2:length(xseq_new)-1
        ekm=norm(xseq_new(:,i-1)-x_min);
        ek=norm(xseq_new(:,i)-x_min);
        ekp=norm(xseq_new(:,i+1)-x_min);
        pseq(i-1)= log( ekp/ek) / log( ek/ekm );
    end
    p=median(pseq)
end

end

%% 7) EXPERIMENTAL RATE OF CONVERGENCE WITHOUT MINIMUN POINT COMPUTATION 

function [p,pseq] = exp_conv_rate(k, xseq)
% function to compute the experimental rate of convergence of a method,
% given:
% x_min the minimum of the function
% k = number of iterations of the method until convergence
% xseq = sequence of points produced by the iterative method
% The function returns:
% pseq = the sequence of experimental rates computed for 3 consecutive
%   values of the sequence x(k-1), x(k), x(k+1)
% p = the median of the convergence rate vector
if k<=100
    pseq=zeros(k,1);
    for i=3:size(xseq,2)-1
        ek=norm( xseq(:,i)-xseq(:,i-1) ); % xk-x(k-1)
        ekm=norm( xseq(:, i-1)-xseq(:,i-2) ); %x(k-1)-x(k-2)
        ekp=norm(xseq(:,i+1)-xseq(:,i)); %x(k+1)-x(k)
        pseq(i-2)= log( ekp/ek ) / log( ek/ekm ) ;
    end
    p=mean(pseq)
else
    xseq_new=xseq(:,size(xseq,2)-100:size(xseq,2));
    pseq=zeros(100,1);
    for i=3:size(xseq_new,2)-1
        ek=norm( xseq_new(:,i)-xseq_new(:,i-1) );
        ekm=norm( xseq_new(:,i-1)-xseq_new(:,i-2) );
        ekp=norm( xseq_new(:,i+1)-xseq_new(:,i) );
        pseq(i-2)= log( ekp/ek ) / log( ek/ekm );
    end
    p=mean(pseq)
end

end
