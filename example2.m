classdef example2
    properties
        xa  
        xb 
        ya 
        yb 
    end
    methods
        function obj = example2()
            obj.xa = -2; 
            obj.xb = 2;
            obj.ya = -2; 
            obj.yb = 2;
        end
        function phi = angle(obj, r)
            if r > 3
                phi = pi;
            elseif r > 2
                phi = 1.5*pi*(r-2)^2 - pi*(r-2)^3 + 0.5*pi;
            elseif r > 1
                phi = 0.5*pi;
            else
                phi = 1.5*pi*(r^2) - pi*(r^3);
            end
        end

        function f = f_x(obj,x,y)
            r = x^2 + y^2;
            phi = obj.angle(r);
            f = cos(phi)*x - sin(phi)*y;
        end

        function f = f_y(obj,x,y)
            r = x^2 + y^2;
            an = @(r) obj.angle(r);
            phi = arrayfun(an,r);
            f = sin(phi).*x + cos(phi).*y;
        end

        function plotflow(obj)
            fx = @(x,y) obj.f_x(x,y);
            fy = @(x,y) obj.f_y(x,y);
            
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
            axis square
            tile=tile+1;
            nexttile(tlo_grad_dec)
            title('Trajectories');
            axis square
            for x = obj.xa:0.5:obj.xb
                for y = obj.ya:0.5:obj.yb
                    if x==obj.xa || x==obj.xb || y==obj.ya || y==obj.yb || x==0.5 || y==0.5 || x==-0.5 || y==-0.5
                        object=streamline(xx,yy,U,V,x,y);
                        object.Color='blue';
                    end
                end
            end
        end
    end
end