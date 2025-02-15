% choose example
example = 2;                  % 1 or 2
obs_operator = 'C';           % 'C' or 'D'
size_of_obs_operator = 9;
end_point_observation = false; %true or false
t_final = 1;

% set parameters for the discretization
nx = 30;
ny = 30;
nt = 200;

% define weight function
    %Example 1
    %w = @(x,y) (x^2 + y^2 + 1)/10;
    
    %Example 2
    %g = @(x) (x^2)/9 -2*x + 10;
    %w = @(x,y) g(x^2 + y^2);
    
    %No weight
    w = @(x,y) 0.1;
%--------define system--------
if end_point_observation
    system = system_for_endpoint_obs(example,t_final,nx,ny,nt);
    system = system.set_w(w);
else
    system = system_for_trajectory_observation(example,t_final,nx,ny,nt);
    system.nt_Lyap = 50;
    system = system.set_w(w);
    system = system.set_C(obs_operator,size_of_obs_operator);
end
%----------------------------

%set parameters for the cost functional
alpha = 0.00000001;
beta  = 0.00001;
system = system.set_alpha_beta(alpha,beta);


%---------define correct solution and generate observation-------
x0 = [-1.5,-1.5];
%x1 = [1,1];
%x2 = [1.5,-0.5];
q = @(x,y) exp(-150*(x-x0(1)).^2-150*(y-x0(2)).^2) ...
    ;%+ exp(-150*(x-x1(1)).^2-150*(y-x1(2)).^2);%+ exp(-150*(x-x2(1)).^2-150*(y-x2(2)).^2);

if end_point_observation
    q_ar = arrayfun(q,system.xx,system.yy);
    q_ar = L1normalize(q_ar,system);
    obs = system.final_state(q_ar);
else
    obs = system.C_for_time(system.solution_time(arrayfun(q,system.xx,system.yy)));
end
system = system.set_observation(obs);

%---------settings for gradient decent algorithm---------

%choose projection for the gradient decent algortihm
L1_sphere = false;     	                % if true, states are projected 
system = system.set_LoneBd(L1_sphere);  % on the Lone-unit ball


if end_point_observation
    type = 'ep';
else
    type = ['t', obs_operator];
end
filename = strcat("data_", type);
filename = strcat(filename,"_nx_");
filename = strcat(filename,string(nx));
%filename = strcat(filename,"_ny_");
%filename = strcat(filename,string(ny));
%filename = strcat(filename,"_nt_");
%filename = strcat(filename,string(nt));
%filename = strcat(filename,"_ex_");
%filename = strcat(filename,string(example));
%filename = strcat(filename,".mat");
filename="data_testvorzeichen.mat"
save(filename)
filename = strcat(filename," saved \n");
fprintf(filename)

function normq = L1normalize(q,sys)
    n = norm((q(:)),1)*sys.hx*sys.hy;
    normq = q/n;
end

