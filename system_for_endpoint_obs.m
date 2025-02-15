classdef system_for_endpoint_obs
    properties
        % Space-time domain
        xa % left boundary of the space domain
        xb % right boundary of the space domain
        ya % lower boundary of the space domain
        yb % upper boundary of the space domain
        tfinal % final time

        % Discretization parameters
        nx    % number of discretization points in x
        ny    % number of discretization points in y
        nt    % number of discretization points in time

        % Definiton of the dynamics
        fx % function handle for the x-component of the dynamics
        fy % function handle for the y-component of _nrthe dynamics

        % Definition of the cost functional
        obs   = []; % the desired observation given as Matrix or function handle
        alpha = 0;  % weight for the L2 norm of the control
        beta  = 0;  % weight for the L1 norm of the control
        J           % function handle for the cost functional

        % Definition of the space of admissible controls
        w = []  % weight function for the L^p_w norm
        LoneBd  % eitehr true or false. If LoneBd==true, the control has to satisfy ||q||_{L^1} <= 1
                % In both cases it has to be nonnegative
        Lambda_one_over_w_squared % necessary value for the projection onto the L^1 sphere

        % discretization
        A  % discretized generator of the predual semigroup
        A_star % discretized generator of the Koopman semigroup
        xd % discretization points in x
        yd % discretization points in y
        xx % grid in x
        yy % grid in y	
        hx % grid spacing in x
        hy % grid spacing in y

        % break condition in gradient descent iteration
        max_iter = 0 % maximum number of iterations
        min_update  % minimum update
 
        % minimal step size 
        step_size_epsilon

        % options for the ODE solver
        options = odeset('RelTol',1e-6);

        % plot options
        plot_cost = false;
        plot_steps = false;   
    end
    methods
        function obj = system_for_endpoint_obs(example,tfinal,nx,ny,nt)
            % Constructor
            if example == 1
                exmpl = example1(1,3);
            elseif example == 2
                exmpl = example2();
            else
                error('Example must be 1 or 2')
            end
            obj.xa = exmpl.xa;
            obj.xb = exmpl.xb;
            obj.ya = exmpl.ya;
            obj.yb = exmpl.yb;
            obj.tfinal = tfinal;
            obj.nx = nx;
            obj.ny = ny;
            obj.nt = nt;
            obj.fx = @(x,y) exmpl.f_x(x,y);
            obj.fy = @(x,y) exmpl.f_y(x,y);
            obj.max_iter = 0;
            [obj.A_star,obj.xd,obj.yd,obj.hx,obj.hy] = no_bc_upwind(obj.xa,obj.xb,obj.ya,obj.yb,nx,ny,obj.fx,obj.fy);
            [obj.xx, obj.yy] = meshgrid(obj.xd, obj.yd);

            obj.step_size_epsilon= 10^(-8);
            obj.min_update = 1e-12*obj.hx*obj.hy;
        end

        function obj = set_alpha_beta(obj, alpha, beta)
            obj.alpha = alpha;
            obj.beta = beta;
            if  ~isempty(obj.obs) & ~isempty(obj.w)
                obj.J = obj.cost_functional(obj.obs, obj.alpha, obj.beta);
            end
        end

        function obj = set_observation(obj, obs)
            if isa(obs, 'function_handle')
                obj.obs = arrayfun(obs, obj.xx,obj.yy); 
            else
                obj.obs = obs;
            end
            if  ~isempty(obj.w)
                obj.J = obj.cost_functional(obj.obs, obj.alpha, obj.beta);
            end
        end

        function obj = set_w(obj,w)
            if isa(w, 'function_handle')
                obj.w = arrayfun(w, obj.xx,obj.yy); 
            else
                obj.w = w;
            end
            if obj.LoneBd
                obj.Lambda_one_over_w_squared = obj.final_state(obj.w.^(-2));
            end
            if  ~isempty(obj.obs)
                obj.J = obj.cost_functional(obj.obs, obj.alpha, obj.beta);
            end
            % define the predual of A_star w.r.t. the L^2_(w^2) inner product
            obj.A = diag(obj.w(:).^(-2))*(obj.A_star')*diag(obj.w(:).^2);
             
        end

        function obj = set_LoneBd(obj, LoneBd)
            if isa(LoneBd, 'logical')
                obj.LoneBd = LoneBd;
            else
                error('LoneBd must be ture or false')
            end
            if LoneBd
                obj.Lambda_one_over_w_squared = obj.final_state(obj.w.^(-2));
            end
        end
        
        function J = cost_functional(obj, obs, alpha, beta)
            % define the cost functional J(q) = 1/2*||y-obs||_{L2}^2 + alpha/2*||q||_{L2}^2 + beta/2*||q||_{L1}
            J = @(y,q) (0.5*obj.L2normSquared(y-obs) + 0.5*alpha*obj.L2normSquared(q) + beta*obj.L1norm(q));
        end

        function prod = prod(~, A,q)
            % compute the product of the linear operator A on a state q
            prod = reshape(A*q(:), size(q));
        end

        function int = integral(obj, q)
            % compute the integral of a state q over Omega
            int = sum(q(:))*obj.hx*obj.hy;
        end

        function n = L2normSquared(obj, q)
            % compute the squared L^2_w norm of a state q
            n = norm(q.*obj.w,'fro')^2*obj.hx*obj.hy;
        end

        function v = inner_product_L2w(obj, q1, q2)
            % compute the L^2_w^2 inner product of two states q1 and q2
            v = sum(q1(:).*q2(:).*obj.w(:).^2)*obj.hx*obj.hy;
        end

        function n = L1norm(obj, q)
            % compute the L^1_w norm of a state q
            n = norm(reshape(q.*obj.w,1,[]),1)*obj.hx*obj.hy;
        end

        function signum = signum(~, q)
            % compute the signum of a state q where we define sgn(0)= -1
            if q > 0
                signum = 1;
            else
                signum = -1;
            end
        end

        function v = variational_minimization(obj, Lambdaq, q, lambda)
            % compute the value of the term of the variational inequality
            v = obj.inner_product_L2w(Lambdaq - obj.obs, Lambdaq) - obj.inner_product_L2w(obj.alpha*q + lambda, q);
        end

        function q = projection_plus(~, q)
            % projects q onto the convex  set of nonnegative states
            q = max(0, q);
        end

        function [lambda, p] = projection_Lone(obj, q)
            % projects q in L^2_{w^2} onto the L^1_w sphere
            k=-2;
            lambda = 2^(-k+1);
            interval=0;
            q=max(0,q);
            int=obj.L1norm(q);
            if int <= 1
                p = max(0,q);
                lambda = 0;
            else
                while k < 20
                    x=max(0, q-lambda*(obj.w.^(-2)));
                    int=obj.L1norm(x);
                    if int > 1
                        if interval==1
                            lambda=lambda + (2^(-k));
                            k = k + 1; 
                        else
                            lambda=2*lambda;
                            k = k - 1;
                        end
                    else
                        if interval==0
                            interval=1;
                        end
                        lambda=lambda - (2^(-k));
                        k=k + 1;
                    end
                end
                p = max(0,q-lambda*(obj.w.^(-2)));
            end
        end

        function y = final_state(obj, q)
            % calculate the end state corresponding to a control q
            % L^2(Omega) --> L^2(Omega), q |---> y(T)
            [~,y] = ode45(@(t, x) (obj.A*x), linspace(0, obj.tfinal,obj.nt), q, obj.options);
            y = reshape(y(end, :),[obj.nx,obj.ny]);
        end

        function p = direction_of_maximal_descent(obj, q, y)
            % calculate the adjoint state p as the solution of the adjoint equation
            [~,p] = ode45(@(t, x) (obj.A_star*x),  linspace(0, obj.tfinal, obj.nt), y-obj.obs, obj.options);
            p = reshape(p(end, :),[obj.nx, obj.ny]) + obj.alpha*q + obj.beta*obj.signum(q).*(obj.w.^(-1));
            p = p.*obj.w.^2;
        end

        function [p, y, new_cost] = gd_step(obj, q, dir, cost_at_q)   
            if obj.LoneBd
                [p, y, new_cost] = obj.gd_step_Lone(q, dir, cost_at_q);
            else
                [p, y, new_cost] = obj.gd_step_nonnegative(q, dir, cost_at_q);
            end
        end

        function [p, y, new_cost] = gd_step_Lone(obj, q, dir, cost_at_q)
            % definition of the first step size
            norm_d=obj.L2normSquared(dir);
            step_size = 0.5*(obj.L2normSquared(q))/norm_d*10;
            % calculate the step size
            Stop = 0;
            ST_q=obj.final_state(q);
            ST_dir=obj.final_state(dir);
        
            while Stop == 0
                [lambda, p] = obj.projection_Lone(q + step_size*dir);
                y = obj.projection_plus(ST_q + step_size*ST_dir - lambda*obj.Lambda_one_over_w_squared);
                lam = obj.beta*obj.signum(p).*(obj.w.^(-1));
                new_cost = obj.variational_minimization(y, p, lam);
                if new_cost > (cost_at_q - 0.2*step_size*norm_d) && step_size > obj.step_size_epsilon
                    step_size = step_size/2;
                else
                    Stop = 1;
                end
            end
        end

        function [p, y, new_var] = gd_step_nonnegative(obj, q, dir, cost_at_q)
            % definition of the first step size
            norm_d=obj.L2normSquared(dir);
            step_size = 10.5*(obj.L2normSquared(q))/norm_d;
            step_size = 10;
            % calculate the step size
            Stop = false;
            Lambda_q=obj.final_state(q);
            Lambda_dir=obj.final_state(dir);

            while ~Stop
                p = obj.projection_plus(q + step_size*dir);
                y = obj.projection_plus(Lambda_q + step_size*Lambda_dir);
                lam = obj.beta*obj.signum(p).*(obj.w.^(-1));
                %new_cost = obj.J(y, p);
                new_var = obj.variational_minimization(y, p, lam);
                if new_var > (cost_at_q - 0.2*step_size*norm_d) & step_size > obj.step_size_epsilon
                    step_size = step_size/2;
                else
                    Stop = true;
                end
            end
        end

        function q = gradient_descend(obj,initial_guess)
            %test if the maximum number of iterations is set
            if obj.max_iter == 0
                fprintf('Please set the maximum number of iterations first')
                return
            end
            
            % Define the initial control q_0 to start the iteration with, the initial end state y_0 and the initial cost J_0 
            if isa(initial_guess, 'function_handle')
                q = arrayfun(initial_guess, obj.xx, obj.yy);
            else
                q = initial_guess;
            end
            y = obj.final_state(q);
            %J_new = obj.J(y, q);
            lam = obj.beta*obj.signum(q).*(obj.w.^(-1));
            var_val = obj.variational_minimization(y, q, lam);
            costlist = [];
            iter = 0;
            update = 2*obj.min_update;

            %%%plotting
            if obj.plot_steps
                fig_grad_dec=figure();
                tlo_grad_dec=tiledlayout(4, 3);
                h_grad_dec=gobjects(1,12);
                tile=1;
                nexttile(tlo_grad_dec)
                h_grad_dec(tile)=surf(obj.xd,obj.yd,q);
                xlabel('x')
                ylabel('y')
                zlabel('q')
                title('Initial guess q_0');
            end

            %%% start the iteration
            while iter < obj.max_iter %& update > obj.min_update
                %save the old cost
                %J_old = J_new;
                var_val_old = var_val;
                % calculate the adjoint state p, the subdifferential lambda and the direction of maximal descent
                dir = -obj.direction_of_maximal_descent(q, y);
                [q, y, var_val] = obj.gd_step(q, dir, var_val_old);
                %J_new = obj.J(y, q);
                update = abs(var_val - var_val_old);
                % update the iteration counter
                iter = iter + 1;
                %%%plotting each control
                % plotting_gap= floor((obj.max_iter-1)/11);
                % if mod(iter, plotting_gap) == 0
                if ismember(iter,[1,5,10,50,100,500,1000,2500,5000,75000,10000]) && obj.plot_steps
                    tile=tile+1;
                    nexttile(tlo_grad_dec)
                    h_grad_dec(tile)=surf(obj.xd,obj.yd,q);
                    xlabel('x')
                    ylabel('y')
                    zlabel('q')
                    title('Result after iteration Nr.: ' + string(iter));
                end
                if obj.plot_cost
                    costlist = [costlist, J_new];
                end
                fprintf('Iteration: %d, Update: %f, New var_val: %g \n ', iter, update, var_val)
            end
            if iter == obj.max_iter
                fprintf('Maximum number of iterations reached')
            end
            fprintf('done')
            % plotting the evolution of the cost 
            if obj.plot_cost
                figure
                plot(costlist)
                xlabel('Iteration')
                ylabel('Cost')
                title('Evolution of the cost')
            end
        end
    
        function q = inverse_flow(obj)
            %discretize the Koopman operator for the inverse flow
            fx_inv = @(x,y) -obj.fx(x,y);
            fy_inv = @(x,y) -obj.fy(x,y);
            [A_star_inv,~,~,~,~] = no_bc_upwind(obj.xa,obj.xb,obj.ya,obj.yb,obj.nx,obj.ny,fx_inv,fy_inv);
            %calculate the corresponding preadjoint
            A_inv = diag(obj.w(:).^(-2))*(A_star_inv')*diag(obj.w(:).^2);
            %calculate the control q corresponding to the observation
            linear_function = @(t, x) (A_inv*x);
            tspan = linspace(0, obj.tfinal, obj.nt);
            [~,q] = ode45(linear_function, tspan, obj.obs, obj.options);
            q = reshape(q(end, :),[obj.nx,obj.ny]);
        end
    end
    
end
