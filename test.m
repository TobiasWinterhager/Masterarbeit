% choose example
example = 2;                  % 1 or 2
obs_operator = 'C';           % 'C' or 'D'
size_of_obs_operator = 8;
end_point_observation = true; %true or false
t_final = 100;

% set parameters for the discretization
nx = 10;
ny = 10;
nt = 10;

% load data
if end_point_observation
    type = 'ep';
else
    type = ['t', obs_operator];
end
filename = strcat("data_", type);
filename = strcat(filename,"_nx_");
filename = strcat(filename,string(nx));
filename = strcat(filename,"_ny_");
filename = strcat(filename,string(ny));
filename = strcat(filename,"_nt_");
filename = strcat(filename,string(nt));
filename = strcat(filename,"_ex_");
filename = strcat(filename,string(example));
filename = strcat(filename,".mat");
fprintf(filename)
load(filename);

alpha = 0.00000001;
beta  = 0.00001;
system = system.set_alpha_beta(alpha,beta);

% ---------define what to plots here----------------------------
%system.plot_cost = true;
%tiles(q,system)
%plot_eigenvalues(system)
%inv_flow_for_final_times_plots()
%inv_flow_for_nt_plots()
%create_plots()
%test_for_different_alpha_beta_plots()
%test_for_w_plots()
%create_plots_observation()
%optcon_for_obs()
%inv_flow_for_nx_plots()
%test_w_1()
%test_L1bd()
%plot_flow()

%------------------

% functions for plotting
function plt = plot_state_heatmap(name,q,system)
    plt = surf(system.xd,system.yd,q,'LineStyle','none');
    
    xl = xlabel('x');
    xl.Position(1) = system.xb + 0.3;
    xl.Position(2) = system.ya - 0.1;
    yl = ylabel('y');
    yl.Position(1) = system.xa - 0.4;
    yl.Position(2) = system.yb + 0;
    yl.Rotation = 0;
    yl.FontWeight = 'bold';
    xl.FontWeight = 'bold';
    zlabel('q')
    name = string(name);
	ttl = title([name], 'Fontsize', 12);
	ttl.Position(1) = ttl.Position(1) - 1;
	ttl.Position(2) = ttl.Position(2) + 1;
    set(gca, 'xlim', [system.xa system.xb]);
    set(gca, 'ylim', [system.ya system.yb]);
    box off
    axis square
    xO = 0.2;  
	yO = 0.1;
	patch(...
		[system.xb system.xa ; system.xa system.xa ], ...
		[system.ya system.yb ; system.ya system.ya ], ...
		[5000 5000;5000 5000],'k', 'clipping', 'off','LineWidth',2)
	patch(...
		[system.xb system.xa-yO; system.xb system.xa+yO; system.xb+xO system.xa], ...
		[system.ya+yO system.yb; system.ya - yO system.yb; system.ya system.yb+xO],...
		[5000 5000;5000 5000;5000 5000], 'k', 'clipping', 'off')
    view(0,90)
    
end

function plt = plot_observation(name,c,system)
    if isa(system, 'system_for_endpoint_obs')
		plt = surf(system.xd,system.yd,c,'LineStyle','none');
        xl = xlabel('x');
		xl.Position(1) = system.xb + 0.3;
		xl.Position(2) = system.ya - 0.1;
		yl = ylabel('y');
		yl.Position(1) = system.xa - 0.4;
		yl.Position(2) = system.yb + 0;
		yl.Rotation = 0;
		yl.FontWeight = 'bold';
		xl.FontWeight = 'bold';
		zlabel('q');
		class(name);
		name = string(name);
		ttl =title([name], 'Fontsize', 12);
		ttl.Position(1) = ttl.Position(1) - 1;
		ttl.Position(2) = ttl.Position(2) + 1;
		set(gca, 'xlim', [system.xa system.xb]);
		set(gca, 'ylim', [system.ya system.yb]); 
		%plt.EdgeColor = 'none';
		box off;
		axis square;
		xO = 0.2;  
		yO = 0.1;
		patch(...
		[system.xb system.xa ; system.xa system.xa ], ...
		[system.ya system.yb ; system.ya system.ya ], ...
		[5000 5000;5000 5000],'k', 'clipping', 'off','LineWidth',2)
		patch(...
		[system.xb system.xa-yO; system.xb system.xa+yO; system.xb+xO system.xa], ...
		[system.ya+yO system.yb; system.ya - yO system.yb; system.ya system.yb+xO],...
		[5000 5000;5000 5000;5000 5000], 'k', 'clipping', 'off')
		view(0,90)
    else
        data= zeros(system.M,system.nt);
        for i = 1:system.nt
            data(:,i) = c{i};
        end
        plt = plot(data');
        name = string(name);
		title([name; string(' ')], 'Fontsize', 12);
        box off;
    end
end

function tiles(q,system,pngname)
    correct_solution = arrayfun(q,system.xx,system.yy);
    obs              = system.obs;
    if isa(system, 'system_for_endpoint_obs')
        sol_inverse    = system.inverse_flow();
        initial_guess   = @(x,y) 1;
        initial_guess   = arrayfun(initial_guess,system.xx,system.yy);
        initial_guess   = L1normalize(initial_guess,system);
        system.max_iter = 250;
        sol_opt_contr   = system.gradient_descend(initial_guess);
        obs_inverse     = system.final_state(sol_inverse);
        obs_opt_contr   = system.final_state(sol_opt_contr);
    else
        nr_ev = 10;
        sol_inverse     = real(system.solve_with_eigenvectors(nr_ev))
        initial_guess   = @(x,y) 1;
        initial_guess   = arrayfun(initial_guess,system.xx,system.yy);
        system.max_iter = 5000;%500;
        sol_opt_contr   = system.gradient_descend(initial_guess);
        obs_inverse     = system.C_for_time(system.solution_time(sol_inverse));
        obs_opt_contr   = system.C_for_time(system.solution_time(sol_opt_contr));
    end

    fig = figure;
    sb(1) = subplot(3,2,1);
    plot_state_heatmap('correct solution',correct_solution, system);
    grid off

    sb(2) = subplot(3,2,2);
    plot_observation("exact observation", obs,system); 
    grid off

    sb(3) = subplot(3,2,3);
    plot_state_heatmap('solution inv Problem',sol_inverse,system);
    grid off

    sb(4) = subplot(3,2,4);
    plot_observation('observation inv Problem', obs_inverse,system);
    grid off

    sb(5) = subplot(3,2,5);
    plot_state_heatmap('solution optimal control',sol_opt_contr,system);
    grid off

    sb(6) = subplot(3,2,6);
    plot_observation('obs. optimal control',obs_opt_contr,system);
    grid off

    clims  = [clim(sb(1)); clim(sb(3)); clim(sb(5))];
    limits = [min(clims(:,1)) max(clims(:,2))];
    clim(sb(1), limits);
    k=2;
    if isa(system, 'system_for_endpoint_obs')
		clims  = [clim(sb(2)); clim(sb(4)); clim(sb(6))];
		limits = [min(clims(:,1)) max(clims(:,2))];
		clim(sb(2), limits);
		k=1;
    end
    size = 0.25;
    for i = 1:k:6
		sb(i).Position(3) = size;
		sb(i).Position(4) = size;
    end
    ax = axes(fig,'visible','off');
    cb = colorbar(sb(1),'Position',[0.01 sb(5).Position(2)-0.05 0.05 sb(1).Position(2)+size+0.01-sb(5).Position(2)],'AxisLocation', 'in');
    if isa(system, 'system_for_endpoint_obs')
        cb = colorbar(sb(2),'Position',[0.51 sb(5).Position(2)-0.05 0.05 sb(1).Position(2)+size+0.01-sb(5).Position(2)],'AxisLocation', 'in');
        fig.Position = [20 20 600 700];
        for i = [2 4 6]
            sb(i).Position(1) = sb(i).Position(1) + 0.07;
        end
        
        for i = 1:6
            sb(i).Position(1) = sb(i).Position(1) + 0.03;
            sb(i).Position(2) = sb(i).Position(2) - 0.05;
        end
    else
        fig.Position = [20 20 1100 900];
        for i= 2:2:6
            sb(i).Position(3) = 2*size;
            sb(i).Position(4) = size;
        end
        for i = [2 4 6]
            sb(i).Position(1) = sb(i).Position(1) - 0.2;
        end
        %cb.Position(1) = cb.Position(1)-0.16;
        for i = 1:6
            sb(i).Position(1) = sb(i).Position(1) - 0.07;
            %sb(i).Position(2) = sb(i).Position(2) + 0.05;
            cb.Position(1) = cb.Position(1);
        end
    end

     
    for i = [2 4 6]
		sb(i).Position(1) = sb(i).Position(1) + 0.07;
    end
	
	for i = 1:6
		sb(i).Position(1) = sb(i).Position(1) + 0.03;
		sb(i).Position(2) = sb(i).Position(2) - 0.05;
	end

    saveas(fig,pngname)

    if isa(system, 'system_for_endpoint_obs')
        fprintf('Error sol_inverse \n')
        printerrors(sol_inverse,correct_solution,system)
        fprintf('Error obs_inverse \n')
        printerrors(obs_inverse,obs,system)
        fprintf('Error sol_opt_contr \n')
        printerrors(sol_opt_contr,correct_solution,system)
        fprintf('Error obs_opt_contr \n')
        printerrors(obs_opt_contr,obs,system)
    end    
end

% plots the eigenvalues of the observability Gramian
function plot_eigenvalues(system)
    nr_ev = 8;
    rows = nr_ev/2;
    fig = figure;
    clims =[];
    
    for row = 1: rows
        i = 4*(row-1);
        e_i = reshape(system.V(:,i/2 + 1),system.nx,system.ny);
        e_i = e_i/sqrt(system.L2normSquared(e_i));
        obs_e_i = system.C_for_time(system.solution_time(e_i));

        sb(i+1) = subplot(rows,4,i+1);
        plot_state_heatmap(strcat('eigenvalue Nr.', string(i/2 + 1)),e_i,system);
        clims = [clims; clim(sb(i+1))];
        colorbar;

        sb(i+2) = subplot(rows,4,i+2);
        plot_observation(strcat('obs of eigenvalue Nr.', string(i/2 + 1)),obs_e_i,system);

        e_i = reshape(system.V(:,i/2 + 2),system.nx,system.ny);
        e_i = e_i/sqrt(system.L2normSquared(e_i));
        obs_e_i = system.C_for_time(system.solution_time(e_i));

        sb(i+3) = subplot(rows,4,i+3);
        plot_state_heatmap(strcat('eigenvalue Nr.', string(i/2 + 2)),e_i,system);
        clims = [clims; clim(sb(i+3))];
        colorbar;

        sb(i+4) = subplot(rows,4,i+4);
        plot_observation(strcat('obs of eigenvalue Nr.', string(i/2 + 2)),obs_e_i,system);

    end
    clim(sb, [min(clims(:,1)) max(clims(:,2))])
    % Link the color limits of both axes
    linkprop(sb(1:2:2*nr_ev), 'CLim');

end

% normalizes the initial guess
function normq = L1normalize(q,sys)
    n = norm((q(:)),1)*sys.hx*sys.hy;
    normq = q/n;
end
   
function printerrors (a,b,sys)
    q=a-b;
    q=q(:);
    %L1
    e = norm(q,1)*sys.hx*sys.hy;
    fprintf('Error in L1 norm is %d \n',e)
    %L1w
    q=(a-b).*sys.w;
    q=q(:);
    e = norm(q,1)*sys.hx*sys.hy;
    fprintf('Error in L1w norm is %d \n',e)
    %L2
    q=(a-b).^2;
    e = norm(reshape((q),1,[]),1)*sys.hx*sys.hy;
    fprintf('Error in L2 norm is %d \n',e)
    %L2w
    e = norm(reshape((q.*sys.w.^2),1,[]),1)*sys.hx*sys.hy;
    fprintf('Error in L2w norm is %d \n',e)
end

% creates the plots form the general results subsection
function create_plots()
    for i = 1:4    
        filename = strcat("example_", string(i), "_for_general_results.mat");
        load(filename)
        alpha = 0.000000001;
        beta  = 0.00001;
        system = system.set_alpha_beta(alpha,beta);
        pngname = strcat('example_', string(i), '_for_general_results.png');
        tiles(q,system,pngname)
    end
end

function create_plots_observation()
    for i = 1:1
        for j =1:1
            %if j == 1
            %    filename = strcat("data_for_Ex", string(i), "C.mat");
            %    load(filename)
            %    pngname = strcat('example_', string(i), '_for_obsC.png');
            %else
            %    filename = strcat("data_for_Ex", string(i), "D.mat");
            %    load(filename)
            %    pngname = strcat('example_', string(i), '_for_obsD.png');
            %end
            filename = "data_tC_nx_40_ny_40_nt_1000_ex_2.mat"
            pngname = 'dfs.png'
            load(filename)
            alpha = 0.00000001;
            beta  = 0.00001;
            system = system.set_alpha_beta(alpha,beta);
            tiles(q,system,pngname)
        end
    end
end

function inv_flow_for_final_times_plots()
    for i = [1 2 4 8]    
        filename = strcat("inv_flow_test_tfinal", string(i), ".mat");
        load(filename)
        alpha = 0.000000001;
        beta  = 0.0000001;
        system = system.set_alpha_beta(alpha,beta);
        
        correct_solution = arrayfun(q,system.xx,system.yy);
        obs              = system.obs;
        sol_inverse    = system.inverse_flow();
        obs_inverse     = system.final_state(sol_inverse);
        
        fig = figure;
        sb(1) = subplot(2,2,1);
        ttl = strcat('correct solution T=', string(i));
        plot_state_heatmap(ttl,correct_solution, system);
        grid off

        sb(2) = subplot(2,2,2);
        ttl = strcat('exact observation T=', string(i));
        plot_observation(ttl, obs,system); 
        grid off

        sb(3) = subplot(2,2,3);
        ttl = strcat('app. solution T=', string(i)); 
        plot_state_heatmap(ttl,sol_inverse,system);
        grid off

        sb(4) = subplot(2,2,4);
        ttl = strcat('app. observation T=', string(i));
        plot_observation(ttl, obs_inverse,system);
        grid off

        clims  = [clim(sb(1)); clim(sb(3))];
        limits = [min(clims(:,1)) max(clims(:,2))];
        clim(sb(1), limits);
        k=2;
        if isa(system, 'system_for_endpoint_obs')
            clims  = [clim(sb(2)); clim(sb(4))];
            limits = [min(clims(:,1)) max(clims(:,2))];
            clim(sb(2), limits);
            k=1;
        end
        size = 0.25;
        for j = 1:k:4
            sb(j).Position(3) = size;
            sb(j).Position(4) = size;
        end
        ax = axes(fig,'visible','off');
        cb = colorbar(sb(1),'Position',[0.01 sb(3).Position(2)-0.05 0.05 sb(1).Position(2)+size+0.01-sb(3).Position(2)],'AxisLocation', 'in');
        cb = colorbar(sb(2),'Position',[0.51 sb(3).Position(2)-0.05 0.05 sb(1).Position(2)+size+0.01-sb(3).Position(2)],'AxisLocation', 'in');
        fig.Position = [20 20 600 500];
    
         
        for j = [2 4]
            sb(j).Position(1) = sb(j).Position(1) + 0.07;
        end
        
        for j = 1:4
            sb(j).Position(1) = sb(j).Position(1) + 0.03;
            sb(j).Position(2) = sb(j).Position(2) - 0.05;
        end
        nameofpng = strcat('inv_flow_for_final_times_plots', string(i), '.png')
        saveas(fig,nameofpng)

        if isa(system, 'system_for_endpoint_obs')
            fprintf('Error solution T=%d \n',i)
            printerrors(sol_inverse,correct_solution,system)
            fprintf('Error inverse T=%d \n',i)
            printerrors(obs_inverse,obs,system)
        end        
    end
end

% tests inv_flow for different nt
function inv_flow_for_nt_plots()
    fig = figure;
    
    j=10^4;
    filename = strcat("inv_flow_test_nt", string(j), ".mat");
    load(filename)
    alpha = 0.000000001;
    beta  = 0.0000001;
    system = system.set_alpha_beta(alpha,beta);
    sol = arrayfun(q,system.xx,system.yy);
    obs = system.obs;    
    
    sb(1) = subplot(4,2,1);
    plot_state_heatmap('correct solution',sol, system);
    grid off

    sb(2) = subplot(4,2,2);
    plot_observation("exact observation", obs,system); 
    grid off

    j=10^2;
    filename = strcat("inv_flow_test_nt", string(j), ".mat");
    load(filename)
    alpha = 0.000000001;
    beta  = 0.0000001;
    system = system.set_alpha_beta(alpha,beta);
    sol    = system.inverse_flow();
    obs     = system.final_state(sol);

    sb(3) = subplot(4,2,3);
    plot_state_heatmap('solution nt=100',sol,system);
    grid off
    
    sb(4) = subplot(4,2,4);
    plot_observation('observation nt=100', obs,system);
    grid off
    
    j=10^3;
    filename = strcat("inv_flow_test_nt", string(j), ".mat");
    load(filename)
    alpha = 0.000000001;
    beta  = 0.0000001;
    system = system.set_alpha_beta(alpha,beta);
    sol    = system.inverse_flow();
    obs     = system.final_state(sol);

    sb(5) = subplot(4,2,5);
    plot_state_heatmap('solution nt=1000',sol,system);
    grid off
    
    sb(6) = subplot(4,2,6);
    plot_observation('observation nt=1000', obs,system);
    grid off

    j=10^4;
    filename = strcat("inv_flow_test_nt", string(j), ".mat");
    load(filename)
    alpha = 0.000000001;
    beta  = 0.0000001;
    system = system.set_alpha_beta(alpha,beta);
    sol    = system.inverse_flow();
    obs     = system.final_state(sol);

    sb(7) = subplot(4,2,7);
    plot_state_heatmap('solution nt=10000',sol,system);
    grid off

    sb(8) = subplot(4,2,8);
    plot_observation('observation nt=10000', obs,system);
    grid off

    clims  = [clim(sb(1)); clim(sb(3));clim(sb(5)); clim(sb(7))];
    limits = [min(clims(:,1)) max(clims(:,2))];
    clim(sb(1), limits);
    k=2;
    if isa(system, 'system_for_endpoint_obs')
        clims  = [clim(sb(2)); clim(sb(4));clim(sb(6)); clim(sb(8))];
        limits = [min(clims(:,1)) max(clims(:,2))];
        clim(sb(2), limits);
        k=1;
    end
    size = 0.25;
    for i = 1:k:8
        sb(i).Position(3) = size;
        sb(i).Position(4) = size;
    end
    ax = axes(fig,'visible','off');
    cb = colorbar(sb(1),'Position',[0.01 sb(7).Position(2)-0.05 0.05 sb(1).Position(2)+size+0.01-sb(7).Position(2)],'AxisLocation', 'in');
    cb = colorbar(sb(2),'Position',[0.51 sb(7).Position(2)-0.05 0.05 sb(1).Position(2)+size+0.01-sb(7).Position(2)],'AxisLocation', 'in');
    fig.Position = [20 20 600 900];

        
    for i = [2 4 6 8]
        sb(i).Position(1) = sb(i).Position(1) + 0.07;
    end
    
    for i = 1:8
        sb(i).Position(1) = sb(i).Position(1) + 0.03;
        sb(i).Position(2) = sb(i).Position(2) - 0.05;
    end
    saveas(fig,'inv_flow_for_nt_plots.png')       
end

% tests GDA for different alpha and beta
function test_for_different_alpha_beta_plots()
    fig = figure;
    
    filename = "data_for_test.mat";
    load(filename)
    initial_guess   = @(x,y) 1;
    initial_guess   = arrayfun(initial_guess,system.xx,system.yy);
    initial_guess   = L1normalize(initial_guess,system);
    system.max_iter = 250;
        


    alpha = 0.000000001;
    beta  = 0.0000001;
    system = system.set_alpha_beta(alpha,beta);    
    sol = q_ar;
    obs = system.obs;    
    
    sb(1) = subplot(4,2,1);
    plot_state_heatmap('correct solution',sol, system);
    grid off

    sb(2) = subplot(4,2,2);
    plot_observation("exact observation", obs,system); 
    grid off
    
    alpha = 0.000000001;
    beta  = 0.0000001;
    system = system.set_alpha_beta(alpha,beta);
    sol    = system.gradient_descend(initial_guess);
    obs    = system.final_state(sol);

    sb(3) = subplot(4,2,3);
    plot_state_heatmap('solution \beta=10^{-7}' ,sol,system);
    grid off

    sb(4) = subplot(4,2,4);
    plot_observation('observation \beta=10^{-7}', obs,system);
    grid off
    
    alpha = 0.000000001;
    beta  = 0.00001;
    system = system.set_alpha_beta(alpha,beta);
    sol    = system.gradient_descend(initial_guess);
    obs     = system.final_state(sol);

    sb(5) = subplot(4,2,5);
    plot_state_heatmap('solution \beta=10^{-5}' ,sol,system);
    grid off

    sb(6) = subplot(4,2,6);
    plot_observation('observation \beta=10^{-5}', obs,system);
    grid off

    alpha = 0.000000001;
    beta  = 0.001;
    system = system.set_alpha_beta(alpha,beta);
    sol    = system.gradient_descend(initial_guess);
    obs     = system.final_state(sol);

    sb(7) = subplot(4,2,7);
    plot_state_heatmap('solution \beta=10^{-3}' ,sol,system);
    grid off
    
    sb(8) = subplot(4,2,8);
    plot_observation('observation \beta=10^{-3}', obs,system);
    grid off

    clims  = [clim(sb(1)); clim(sb(3));clim(sb(5)); clim(sb(7))];
    limits = [min(clims(:,1)) max(clims(:,2))];
    clim(sb(1), limits);
    k=2;
    if isa(system, 'system_for_endpoint_obs')
        clims  = [clim(sb(2)); clim(sb(4));clim(sb(6)); clim(sb(8))];
        limits = [min(clims(:,1)) max(clims(:,2))];
        clim(sb(2), limits);
        k=1;
    end
    size = 0.25;
    for i = 1:k:8
        sb(i).Position(3) = size;
        sb(i).Position(4) = size;
    end
    ax = axes(fig,'visible','off');
    cb = colorbar(sb(1),'Position',[0.01 sb(7).Position(2)-0.05 0.05 sb(1).Position(2)+size+0.01-sb(7).Position(2)],'AxisLocation', 'in');
    cb = colorbar(sb(2),'Position',[0.51 sb(7).Position(2)-0.05 0.05 sb(1).Position(2)+size+0.01-sb(7).Position(2)],'AxisLocation', 'in');
    fig.Position = [20 20 600 900];

        
    for i = [2 4 6 8]
        sb(i).Position(1) = sb(i).Position(1) + 0.07;
    end
    
    for i = 1:8
        sb(i).Position(1) = sb(i).Position(1) + 0.03;
        sb(i).Position(2) = sb(i).Position(2) - 0.05;
    end
    saveas(fig,'test_for_different_alpha_beta_plots2.png')       
end

% tests GDA for trajectory observation
function optcon_for_obs()
    filename = "data_tC_nx_40_ny_40_nt_1000_ex_2.mat"
    pngname = 'dfs.png'
    load(filename)
    alpha = 0.00000001;
    beta  = 0.00001;
    system = system.set_alpha_beta(alpha,beta);
    nr_ev=10;
    correct_solution = arrayfun(q,system.xx,system.yy);
    obs              = system.obs;
    if isa(system, 'system_for_endpoint_obs')
        sol_inverse    = system.inverse_flow();
        initial_guess   = @(x,y) 1;
        initial_guess   = arrayfun(initial_guess,system.xx,system.yy);
        initial_guess   = L1normalize(initial_guess,system);
        system.max_iter = 250;
        sol_opt_contr   = system.gradient_descend(initial_guess);
        obs_inverse     = system.final_state(sol_inverse);
        obs_opt_contr   = system.final_state(sol_opt_contr);
    else
        %sol_inverse     = real(system.solve_with_eigenvectors(nr_ev))
        initial_guess   = @(x,y) 1;
        initial_guess   = arrayfun(initial_guess,system.xx,system.yy);
        system.max_iter = 50;%500;
        sol_opt_contr   = system.gradient_descend(initial_guess);
        %obs_inverse     = system.C_for_time(system.solution_time(sol_inverse));
        obs_opt_contr   = system.C_for_time(system.solution_time(sol_opt_contr));
    end

    fig = figure;
    sb(1) = subplot(2,2,1);
    plot_state_heatmap('correct solution',correct_solution, system);
    grid off

    sb(2) = subplot(2,2,2);
    plot_observation("exact observation", obs,system); 
    grid off

    sb(3) = subplot(2,2,3);
    plot_state_heatmap('solution optimal control',sol_opt_contr,system);
    grid off

    sb(4) = subplot(2,2,4);
    plot_observation('obs. optimal control',obs_opt_contr,system);
    grid off

    clims  = [clim(sb(1)); clim(sb(3))];
    limits = [min(clims(:,1)) max(clims(:,2))];
    clim(sb(1), limits);
    k=2;
    
    size = 0.25;
    for i = 1:k:4
		sb(i).Position(3) = size;
		sb(i).Position(4) = size;
    end
    ax = axes(fig,'visible','off');
    cb = colorbar(sb(1),'Position',[0.01 sb(3).Position(2)-0.05 0.05 sb(1).Position(2)+size+0.01-sb(3).Position(2)],'AxisLocation', 'in');
    if isa(system, 'system_for_endpoint_obs')
        bla=3;
    else
        fig.Position = [20 20 1100 660];
        for i= 2:2:4
            sb(i).Position(3) = 2*size;
            sb(i).Position(4) = size;
        end
        for i = [2 4]
            sb(i).Position(1) = sb(i).Position(1) - 0.2;
        end
        %cb.Position(1) = cb.Position(1)-0.16;
        for i = 1:4
            sb(i).Position(1) = sb(i).Position(1) - 0.07;
            %sb(i).Position(2) = sb(i).Position(2) + 0.05;
            cb.Position(1) = cb.Position(1);
        end
    end

     
    for i = [2 4]
		sb(i).Position(1) = sb(i).Position(1) + 0.02;
    end
	
	for i = 1:4
		sb(i).Position(1) = sb(i).Position(1) + 0.03;
		sb(i).Position(2) = sb(i).Position(2) - 0.05;
	end

    saveas(fig,pngname)

    if isa(system, 'system_for_endpoint_obs')
        fprintf('Error sol_inverse \n')
        printerrors(sol_inverse,correct_solution,system)
        fprintf('Error obs_inverse \n')
        printerrors(obs_inverse,obs,system)
        fprintf('Error sol_opt_contr \n')
        printerrors(sol_opt_contr,correct_solution,system)
        fprintf('Error obs_opt_contr \n')
        printerrors(obs_opt_contr,obs,system)
    end    
end

% tests inv_flow for different nx
function inv_flow_for_nx_plots()
    fig = figure;
    
    j=80;
    filename = strcat("data_ep_nx_", string(j), ".mat");
    load(filename)
    alpha = 0.000000001;
    beta  = 0.0000001;
    system = system.set_alpha_beta(alpha,beta);
    sol = arrayfun(q,system.xx,system.yy);
    obs = system.obs;    
    
    sb(1) = subplot(4,2,1);
    plot_state_heatmap('correct solution',sol, system);
    grid off

    sb(2) = subplot(4,2,2);
    plot_observation("exact observation", obs,system); 
    grid off

    j=40;
    filename = strcat("data_ep_nx_", string(j), ".mat");
    load(filename)
    alpha = 0.000000001;
    beta  = 0.0000001;
    system = system.set_alpha_beta(alpha,beta);
    sol    = system.inverse_flow();
    obs     = system.final_state(sol);

    sb(3) = subplot(4,2,3);
    plot_state_heatmap('solution nx=40',sol,system);
    grid off
    
    sb(4) = subplot(4,2,4);
    plot_observation('observation nx=40', obs,system);
    grid off
    
    j=60;
    filename = strcat("data_ep_nx_", string(j), ".mat");
    load(filename)
    alpha = 0.000000001;
    beta  = 0.0000001;
    system = system.set_alpha_beta(alpha,beta);
    sol    = system.inverse_flow();
    obs     = system.final_state(sol);

    sb(5) = subplot(4,2,5);
    plot_state_heatmap('solution nx=60',sol,system);
    grid off
    
    sb(6) = subplot(4,2,6);
    plot_observation('observation nx=60', obs,system);
    grid off

    j=80;
    filename = strcat("data_ep_nx_", string(j), ".mat");
    load(filename)
    alpha = 0.000000001;
    beta  = 0.0000001;
    system = system.set_alpha_beta(alpha,beta);
    sol    = system.inverse_flow();
    obs     = system.final_state(sol);

    sb(7) = subplot(4,2,7);
    plot_state_heatmap('solution nx=80',sol,system);
    grid off

    sb(8) = subplot(4,2,8);
    plot_observation('observation nx=80', obs,system);
    grid off

    clims  = [clim(sb(1)); clim(sb(3));clim(sb(5)); clim(sb(7))];
    limits = [min(clims(:,1)) max(clims(:,2))];
    clim(sb(1), limits);
    k=2;
    if isa(system, 'system_for_endpoint_obs')
        clims  = [clim(sb(2)); clim(sb(4));clim(sb(6)); clim(sb(8))];
        limits = [min(clims(:,1)) max(clims(:,2))];
        clim(sb(2), limits);
        k=1;
    end
    size = 0.25;
    for i = 1:k:8
        sb(i).Position(3) = size;
        sb(i).Position(4) = size;
    end
    ax = axes(fig,'visible','off');
    cb = colorbar(sb(1),'Position',[0.01 sb(7).Position(2)-0.05 0.05 sb(1).Position(2)+size+0.01-sb(7).Position(2)],'AxisLocation', 'in');
    cb = colorbar(sb(2),'Position',[0.51 sb(7).Position(2)-0.05 0.05 sb(1).Position(2)+size+0.01-sb(7).Position(2)],'AxisLocation', 'in');
    fig.Position = [20 20 600 900];

        
    for i = [2 4 6 8]
        sb(i).Position(1) = sb(i).Position(1) + 0.07;
    end
    
    for i = 1:8
        sb(i).Position(1) = sb(i).Position(1) + 0.03;
        sb(i).Position(2) = sb(i).Position(2) - 0.05;
    end
    saveas(fig,'inv_flow_for_nx_plots.png')       
end

% tests GDA for constant w=1
function test_w_1()
        
    filename = strcat("data_w=1.mat");
    load(filename)
    initial_guess   = @(x,y) 1;
    initial_guess   = arrayfun(initial_guess,system.xx,system.yy);
    initial_guess   = L1normalize(initial_guess,system);
    system.max_iter = 250;	
    alpha = 0.00000001;
    beta  = 0.00001;
    system = system.set_alpha_beta(alpha,beta);
    
    correct_solution = arrayfun(q,system.xx,system.yy);
    obs              = system.obs;
    sol_inverse    = system.gradient_descend(initial_guess);
    obs_inverse     = system.final_state(sol_inverse);
    
    fig = figure;
    sb(1) = subplot(2,2,1);
    ttl = strcat('correct solution');
    plot_state_heatmap(ttl,correct_solution, system);
    grid off

    sb(2) = subplot(2,2,2);
    ttl = strcat('exact observation');
    plot_observation(ttl, obs,system); 
    grid off

    sb(3) = subplot(2,2,3);
    ttl = strcat('app. solution w=1'); 
    plot_state_heatmap(ttl,sol_inverse,system);
    grid off

    sb(4) = subplot(2,2,4);
    ttl = strcat('app. observation w=1');
    plot_observation(ttl, obs_inverse,system);
    grid off

    clims  = [clim(sb(1)); clim(sb(3))];
    limits = [min(clims(:,1)) max(clims(:,2))];
    clim(sb(1), limits);
    k=2;
    if isa(system, 'system_for_endpoint_obs')
        clims  = [clim(sb(2)); clim(sb(4))];
        limits = [min(clims(:,1)) max(clims(:,2))];
        clim(sb(2), limits);
        k=1;
    end
    size = 0.25;
    for j = 1:k:4
        sb(j).Position(3) = size;
        sb(j).Position(4) = size;
    end
    ax = axes(fig,'visible','off');
    cb = colorbar(sb(1),'Position',[0.01 sb(3).Position(2)-0.05 0.05 sb(1).Position(2)+size+0.01-sb(3).Position(2)],'AxisLocation', 'in');
    cb = colorbar(sb(2),'Position',[0.51 sb(3).Position(2)-0.05 0.05 sb(1).Position(2)+size+0.01-sb(3).Position(2)],'AxisLocation', 'in');
    fig.Position = [20 20 600 500];

        
    for j = [2 4]
        sb(j).Position(1) = sb(j).Position(1) + 0.07;
    end
    
    for j = 1:4
        sb(j).Position(1) = sb(j).Position(1) + 0.03;
        sb(j).Position(2) = sb(j).Position(2) - 0.05;
    end
    nameofpng = "test_for_w.png";
    saveas(fig,nameofpng)

    if isa(system, 'system_for_endpoint_obs')
        fprintf('Error solution T=%d \n',i)
        printerrors(sol_inverse,correct_solution,system)
        fprintf('Error inverse T=%d \n',i)
        printerrors(obs_inverse,obs,system)
    end        
end

% tests the projection on L1 ball
function test_L1bd()
        
    filename = strcat("data_testLone.mat");
    load(filename)
    initial_guess   = @(x,y) 1;
    initial_guess   = arrayfun(initial_guess,system.xx,system.yy);
    initial_guess   = L1normalize(initial_guess,system);
    system.max_iter = 30;	
    alpha = 0.00000001;
    beta  = 0.00001;
    system = system.set_alpha_beta(alpha,beta);
    
    correct_solution = arrayfun(q,system.xx,system.yy);
    obs              = system.obs;
    sol_inverse      = system.gradient_descend(initial_guess);
    obs_inverse      = system.final_state(sol_inverse);
    
    fig = figure;
    sb(1) = subplot(2,2,1);
    ttl = strcat('correct solution');
    plot_state_heatmap(ttl,correct_solution, system);
    grid off

    sb(2) = subplot(2,2,2);
    ttl = strcat('exact observation');
    plot_observation(ttl, obs,system); 
    grid off

    sb(3) = subplot(2,2,3);
    ttl = strcat('app. solution w=1'); 
    plot_state_heatmap(ttl,sol_inverse,system);
    grid off

    sb(4) = subplot(2,2,4);
    ttl = strcat('app. observation w=1');
    plot_observation(ttl, obs_inverse,system);
    grid off

    clims  = [clim(sb(1)); clim(sb(3))];
    limits = [min(clims(:,1)) max(clims(:,2))];
    clim(sb(1), limits);
    k=2;
    if isa(system, 'system_for_endpoint_obs')
        clims  = [clim(sb(2)); clim(sb(4))];
        limits = [min(clims(:,1)) max(clims(:,2))];
        clim(sb(2), limits);
        k=1;
    end
    size = 0.25;
    for j = 1:k:4
        sb(j).Position(3) = size;
        sb(j).Position(4) = size;
    end
    ax = axes(fig,'visible','off');
    cb = colorbar(sb(1),'Position',[0.01 sb(3).Position(2)-0.05 0.05 sb(1).Position(2)+size+0.01-sb(3).Position(2)],'AxisLocation', 'in');
    cb = colorbar(sb(2),'Position',[0.51 sb(3).Position(2)-0.05 0.05 sb(1).Position(2)+size+0.01-sb(3).Position(2)],'AxisLocation', 'in');
    fig.Position = [20 20 600 500];

        
    for j = [2 4]
        sb(j).Position(1) = sb(j).Position(1) + 0.07;
    end
    
    for j = 1:4
        sb(j).Position(1) = sb(j).Position(1) + 0.03;
        sb(j).Position(2) = sb(j).Position(2) - 0.05;
    end
    nameofpng = "test_for_Lone.png";
    saveas(fig,nameofpng)

    if isa(system, 'system_for_endpoint_obs')
        fprintf('Error solution T=%d \n',i)
        printerrors(sol_inverse,correct_solution,system)
        fprintf('Error inverse T=%d \n',i)
        printerrors(obs_inverse,obs,system)
    end        
end

% plots the flow of both example form Chapter 'Numeric Analysis' For the first example we choose alpha = 1 and eta = 3. 
function plot_flow()
    ex1=example1(1,3);
    ex1.plotflow();
    ex2=example2();
    ex2.plotflow();
end