classdef example1
    properties
        xa  
        xb 
        ya 
        yb
        alpha
        eta
    end
    methods
        function obj = example1(alpha, eta)
            obj.xa = -3; 
            obj.xb = 1.5;
            obj.ya = -3; 
            obj.yb = 1.5;
            obj.alpha = alpha;
            obj.eta = eta;
        end
        
        function f = f_x(obj,x,y)
            f = y-obj.alpha*x^3;
        end

        function f = f_y(obj,x,y)
            f = -x-obj.eta*y;
        end

        function plotflow(obj)
            fx = @(x,y) obj.f_x(x,y)
            fy = @(x,y) obj.f_y(x,y)
            
            % Setting the discretization parameters for all three dimensions
            nx = 30; ny = 30; % seems fine enough for very decent results ...
            % discretization 
            hx = (obj.xb-obj.xa)/(nx-1); 
            hy = (obj.yb-obj.ya)/(ny-1);
            xd = obj.xa:hx:obj.xb;
            yd = obj.ya:hy:obj.yb;
            [xx, yy] = meshgrid(xd, yd);
            U = arrayfun(fx,xx,yy);
            V = arrayfun(fy,xx,yy);
            % make figure
            fig_grad_dec=figure();
            title('Flow of Example 1')
            tlo_grad_dec=tiledlayout(2, 1);
            h_grad_dec=gobjects(1,2);
            tile=1;
            nexttile(tlo_grad_dec)
            
            quiver(xx,yy,U,V,'b')
            title('Flow field');
            set(gca, 'xlim', [obj.xa obj.xb]);
            set(gca, 'ylim', [obj.ya obj.yb]);
            axis square
            tile=tile+1;
            nexttile(tlo_grad_dec)
            title('Trajectories');
            set(gca, 'xlim', [obj.xa obj.xb]);
            set(gca, 'ylim', [obj.ya obj.yb]);
            axis square
            for x = obj.xa:0.5:obj.xb
                for y = obj.ya:0.5:obj.yb
                    if x==obj.xa || x==obj.xb || y==obj.ya || y==obj.yb 
                        object=streamline(xx,yy,U,V,x,y);
                        object.Color='blue';
                    end
                end
            end
        end
    end
end