classdef observationC
    properties
        % Properties
        system
        N
        w
    end
    methods
        % Constructor
        function obj = observationC(system, N)
            obj.system = system;
            obj.N = N;
        end
        % Methods
        function c_i = indicator_square(obj, x1, x2, y1, y2)
            function Cxy = c(x,y)  
                if x1 <= x && x < x2 && y1 <= y && y < y2
                    Cxy = 100;
                else
                    Cxy = 0;
                end
            end
            c_i = arrayfun(@(x,y) c(x,y), obj.system.xx, obj.system.yy)./obj.system.w;
        end

        function C = observation_operator(obj)
            C = ones(obj.N^2, obj.system.nx*obj.system.ny);
            xa= obj.system.xa;
            xb= obj.system.xb;
            ya= obj.system.ya;
            yb= obj.system.yb;
            sx = xa: (xb-xa)/obj.N : xb;
            sy = ya: (yb-ya)/obj.N : yb;
            for i = 1:obj.N
                for j = 1:obj.N
                    row =obj.N*(i-1)+j;
                    C_row = obj.indicator_square(sx(i), sx(i+1), sy(j), sy(j+1));
                    C(row,:) = C_row(:)';
                end
            end
            fprintf('Observation operator C_N defined with N=%d \n', obj.N);
        end
    end
end