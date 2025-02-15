classdef system_for_trajectory_observation
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
        nt_Lyap = 1000 % number of discretization points in time for the differential Lyapunov equation

        % Definiton of the dynamics
        fx % function handle for the x-component of the dynamics
        fy % function handle for the y-component of the dynamics

        % Definition of the cost functional
        obs   = []; % the desired observation given as Matrix or function handle
        alpha = 0;  % weight for the L2 norm of the control
        beta  = 0;  % weight for the L1 norm of the control
        J           % function handle for the cost functional

        % Definition of the space of admissible controls
        w = [];      % weight function for the L^p_w norm
        LoneBd  % eitehr true or false. If LoneBd==true, the control has to satisfy ||q||_{L^1} <= 1
                % In both cases it has to be nonnegative
        Lambda_one_over_w_squared % necessary value for the projection onto the L^1 sphere

        % discretization
        A       % discretized generator of the predual semigroup
        A_star  % discretized generator of the Koopman semigroup
        xd      % discretization points in x
        yd      % discretization points in y
        xx      % grid in x
        yy      % grid in y	
        hx      % grid spacing in x
        hy      % grid spacing in y

        % parameters for the case of full trajectory observation
        Lambda_star_f   % the dual of the desired observation, calculated once in the beginning

        M               % length of C
        C = [];         % the discetized operator C, mapping q in L^2 to R^M
                        % should be a cell array of length M, containing (discretized) functions in L^2(Omega)
        P = [];         % Matrix reprasantation of the bilinear form <CS_tP, CS_tq>_{L^2(0,T,L^2_{w^2})}
        V = [];         % eigenvectors of P
        D = [];         % eigenvalues of P

        % break condition in gradient descent iteration
        max_iter = 0;                       % maximum number of iterations
        min_update     % minimum update
 
        % minimal step size 
        step_size_epsilon

        % options for the ODE solver
        options = odeset('RelTol',1e-6);

        % plot options
        plot_cost    = false;
        plot_steps   = false;  
    end
    methods
        function obj = system_for_trajectory_observation(example,tfinal,nx,ny,nt)
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
            [obj.A_star,obj.xd,obj.yd,obj.hx,obj.hy] = no_bc_upwind(obj.xa,obj.xb,obj.ya,obj.yb,nx,ny,obj.fx,obj.fy);
            fprintf('Disctetized generator of the Koopman semigroup calculated\n')
            [obj.xx, obj.yy] = meshgrid(obj.xd, obj.yd);
            obj.hx
            obj.step_size_epsilon= 10^(-8);
            obj.min_update = 1e-6*obj.hx*obj.hy;
        end

        function obj = set_alpha_beta(obj, alpha, beta)
            obj.alpha = alpha;
            obj.beta = beta;
            if  ~isempty(obj.obs) & ~isempty(obj.w)
                obj.J = obj.cost_functional(obj.obs, obj.alpha, obj.beta);
            end
        end

        function obj = set_observation(obj, obs)
            obj.obs = obs;
            if  ~isempty(obj.w)
                obj.J = obj.cost_functional(obj.obs, obj.alpha, obj.beta);
            end
            if  ~isempty(obj.C)
                obj.Lambda_star_f = obj.adjoint_solution_time(obj.C_star_for_time(obj.obs));
            end
        end

        function obj = set_w(obj,w)
            if isa(w, 'function_handle')
                obj.w = arrayfun(w, obj.xx,obj.yy); 
            elseif ismethod(class(w), 'size') & size(w,1) == obj.nx & size(w,2) == obj.ny
                obj.w = w;
            else
                error('w must be a function handle or a matrix of size nx times ny')
            end
            if obj.LoneBd & ~isempty(obj.C)
                obj.Lambda_one_over_w_squared = obj.C_for_time(obj.solution_time(obj.w.^(-2)));
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
            if LoneBd & ~isempty(obj.C)
                obj.Lambda_one_over_w_squared = obj.C_for_time(obj.solution_time(obj.w.^(-2)));
            end
        end

        function obj = set_C(obj, choice, N)
            if choice == 'C'
                obj.C = observationC(obj,N).observation_operator();
                obj.M = N^2;
            elseif choice == 'D'
                obj.C = observationD(obj,N).observation_operator();
                obj.M = 2*N;
            else
                error('Choice must be ''C'' or ''D''')
            end
            obj.P = obj.differential_Lyapunov();
            [obj.V,obj.D] = eig(obj.P);
            fprintf('Eigenvalues of the controlability Gramian calculated\n')
        end
        
        function J = cost_functional(obj, obs, alpha, beta)
            % define the cost functional J(q) = 1/2*||y-obs||_{L2}^2 + alpha/2*||q||_{L2}^2 + beta/2*||q||_{L1}
            J = @(y,q) (0.5*obj.L2l2normTimeSquared(obj.lin_comb_time(1,y,-1,obs)) + 0.5*alpha*obj.L2normSquared(q) + beta*obj.L1norm(q));
        end

        function f = lin_comb_time(obj,scalar_f, f, scalarg, g)
            for i = 1: obj.nt
                f{i} = scalar_f*f{i} - scalarg*g{i};
            end
        end

        function prod = prod(~, A,q)
            % compute the product of the linear operator A on a state q
            prod = reshape(A*q(:), size(q));
        end

        function int = integral(obj, q)
            % compute the integral of a state q over Omega
            int = sum(q(:))*obj.hx*obj.hy;
        end

        function int = integralTime(obj, q)
            int = 0.5* obj.integral(q{1}) + 0.5* obj.integral(q{obj.nt});
            for i = 2 : obj.nt-1
                int = int + obj.integral(q{i});
            end
        end

        function n = L2normSquared(obj, q)
            % compute the squared L^2_w norm of a state q
            n = norm(q.*obj.w,'fro')^2*obj.hx*obj.hy;

        end

        function n = L2l2normTimeSquared(obj,q)
            for i = 1: obj.nt
                q{i} = q{i}'*q{i};
            end
            n = obj.integralTime(q);
        end

        function n = L1norm(obj, q)
            % compute the L^1_w norm of a state q
            n = norm(reshape(q.*obj.w,1,[]),1)*obj.hx*obj.hy;
        end
        
        function v = inner_product_L2w(obj, q1, q2)
            % compute the L^2_w^2 inner product of two states q1 and q2
            v = sum(q1(:).*q2(:).*obj.w(:).^2)*obj.hx*obj.hy;
        end

        function signum = signum(~, q)
            % compute the signum of a state q where we define sgn(0)= -1
            if q > 0
                signum = 1;
            else
                signum = -1;
            end
        end

        function q = projection_plus(~, q)
            % projects q onto the convex  set of nonnegative states
            q = max(0, q);
        end

        function u = projection_plus_time(obj, u)
            for i = 1: obj.nt
                u{i} = max(0, u{i});
            end
        end
        
        function q = normalize(obj,q)
			n = sqrt(obj.L2normSquared(q));
			q=q./n
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

        function Sq = solution_time(obj, q)
            % L^2(Omega) --> L^2(0,T,L^2(Omega)), q |---> y
            [~,y] = ode45(@(t, x) (obj.A'*x), linspace(0, obj.tfinal,obj.nt), q, obj.options);
            Sq = cell(1,obj.nt);
            for i = 1: obj.nt
                Sq{i}= reshape(y(i, :),[obj.nx,obj.ny]);
            end
        end

        function y = adjoint_solution_time(obj,q)
            %L^2(0,T,L^2(Omega)) -->  L^2(Omega)
            frc = (obj.nt-1)/obj.tfinal;
            idx = @(t) round(t*frc);
            x_0 = zeros(obj.nx*obj.ny,1);
            rhs = @(t, x) (obj.A_star*x + q{obj.nt - idx(t)}(:));
            [~,y] = ode45(@(t, x) rhs (t,x),  linspace(0, obj.tfinal, obj.nt), x_0, obj.options);
            y = reshape(y(end, :),[obj.nx,obj.ny]);
        end

        function Cq = C_for_state(obj, q)
            % L^2(Omega) --> l^2
            Cq=obj.C*(q(:).*(obj.w(:).^2));
        end
        
        function C_starf = C_star(obj, f)
            % l^2 --> L^2(Omega)
            C_starf = zeros(obj.nx, obj.ny);
            for i = 1: obj.M
                c_i = reshape(obj.C(i,:),[obj.nx,obj.ny]);
                C_starf = C_starf + f(i)*c_i;
            end
        end

        function frCq = C_for_time(obj, q)
            % L^2(0,T,L^2(Omega)) --> L^2(0,T,l^2)	
            frCq = cell(1, obj.nt);
            for j = 1:obj.nt
                frCq{j}=obj.C_for_state(q{j});
            end
        end

        function C_starf = C_star_for_time(obj, f)
            % L^2(0,T,l^2) --> L^2(0,T,L^2(Omega))
            C_starf = cell(1, obj.nt);
            for i = 1 : obj.nt
                C_starf{i} = obj.C_star(f{i});
            end
        end
        
        function v = variational_minimization(obj, q, lambda)
            % compute the value of the term of the variational inequality
            v = obj.inner_product_L2w(obj.prod(obj.P,q) - obj.Lambda_star_f + obj.alpha*q + lambda, q);
        end

        function p = direction_of_maximal_descent(obj, q)
            p = obj.prod(obj.P,q) - obj.Lambda_star_f + obj.alpha*q + obj.beta*obj.signum(q).*(obj.w.^(-1));
        end

        function [p, new_cost] = gd_step(obj, q, dir, cost_at_q)   
            if obj.LoneBd 
                [p, new_cost] = obj.gd_step_Lone_time(q, dir, cost_at_q);
            else
                [p, new_cost] = obj.gd_step_nonnegative_time(q, dir, cost_at_q);
            end
        end

        function [p, new_var] = gd_step_nonnegative_time(obj, q, dir, cost_at_q)
            % definition of the first step size
            norm_d=obj.L2normSquared(dir);
            step_size = 0.5*(obj.L2normSquared(q))/norm_d*10;
            % calculate the step size
            Stop = false;
            %St_q=obj.solution_time(q);
            %St_dir=obj.solution_time(dir);

            while ~Stop 
                p =obj.projection_plus( q + step_size*dir);
                lam = obj.beta*obj.signum(p).*(obj.w.^(-1));
                new_var = obj.variational_minimization(p, lam);
                if new_var > (cost_at_q - 0.2*step_size*norm_d) && step_size > 0
                    step_size = step_size/2;
                else
                    Stop = true;
                end
            end
        end

        function [p, new_cost] = gd_step_Lone_time(obj, q, dir, cost_at_q)
            % definition of the first step size
            norm_d=obj.L2normSquared(dir);
            step_size = 200.5*(obj.L2normSquared(q))/norm_d*10;
            % calculate the step size
            Stop = false;
        
            while ~Stop
                [lambda, p] = obj.projection_Lone(q + step_size*dir);
                lam = obj.beta*obj.signum(p).*(obj.w.^(-1));
                %Stq_min_lammda_one_over_w_sqr = obj.lin_comb_time(1,St_q, -lambda, obj.Lambda_one_over_w_squared);
                %y_before_projection = lin_comb_time(1,Stq_min_lammda_one_over_w_sqr, step_size, St_dir);
                %y = obj.C_for_time(obj.projection_plus_time(y_before_projection));
                new_cost = obj.variational_minimization(p, lam);%obj.J(y, p);
                if new_cost > (cost_at_q - 0.2*step_size*norm_d) %&& step_size > obj.step_size_epsilon
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
            %y= obj.C_for_time(obj.solution_time(q));
            lam = obj.beta*obj.signum(q).*(obj.w.^(-1)); 
            J_new = obj.variational_minimization(q, lam);
            costlist = [J_new];
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
                J_old = J_new;
                % calculate the adjoint state p, the subdifferential lambda and the direction of maximal descent
                dir = -obj.direction_of_maximal_descent(q);
                [q, J_new] = obj.gd_step(q, dir, J_old);
                update = abs(J_new - J_old);
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
                if ismember( iter, 1:1000:obj.max_iter)
					fprintf(' Iteration: %d, Update: %f, New cost: %g\n', iter, update, J_new)
				end	
            end
            if iter == obj.max_iter
                fprintf('Maximum number of iterations reached')
            end
            % plotting the evolution of the cost 
            if obj.plot_cost
                figure
                plot(costlist)
                xlabel('Iteration')
                ylabel('Cost')
                title('Evolution of the cost')
            end
        end

        function P = differential_Lyapunov(obj)
            C_starC = obj.C'*obj.C*diag(obj.w(:).^2);
            P = BDF_Diff_Lyapunov(-obj.A, C_starC, zeros(obj.nx*obj.ny,obj.nx*obj.ny), 0, obj.tfinal, obj.nt_Lyap);
            
        end

        function q = solve_with_eigenvectors(obj, N)
            N=min(obj.nx*obj.ny,N)
            %V=V(:,1:N);
            %fprintf('%d\n',diag(D))
            %P=V*D(:,1:N);
            %q = zeros(obj.nx*obj.ny,1);
            %N = min(N, obj.nx*obj.ny);
            
            %for i = 1:N
            %    if D(i,i) > 0    
            %        e_i = V(:,i);
            %        q = q + (sum(obj.Lambda_star_f(:).*obj.w(:).^2.*e_i)*obj.hx*obj.hy/D(i,i))*e_i;
            %    end
            %end
            %fprintf('size of P= %d,%d \n',size(P))
            %fprintf('size of Lambdaf= %d,%d \n',size(obj.Lambda_star_f))
            %q = P\obj.Lambda_star_f(:);
            %X=D*V;
            %lam=X\obj.Lambda_star_f(:);
            %q = reshape(V*lam, [obj.nx,obj.ny]);
            q = zeros(obj.nx,obj.ny);
            for i = 1:N
                n_i = obj.inner_product_L2w(obj.V(:,i),obj.V(:,i));
                q = q + reshape(obj.inner_product_L2w(obj.V(:,i),obj.Lambda_star_f)*obj.V(:,i)/(obj.D(i,i)*n_i), [obj.nx,obj.ny])
            end
            q = real(q)
            %fprintf('%d \n',q)
        end
    end
end
