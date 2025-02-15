classdef observationD
    properties
        % Properties
        system
        N
    end
    methods
        % Constructor
        function obj = observationD(system, N)
            obj.system = system;
            obj.N = N;
        end
        % Methods
        function C = observation_operator(obj)
            C = ones(obj.N*2, obj.system.nx*obj.system.ny);
            for i = 1:obj.N
                sinx = arrayfun(@(x,y) sin(i*pi*x), obj.system.xx, obj.system.yy);
                siny = arrayfun(@(x,y) sin(i*pi*y), obj.system.xx, obj.system.yy);
                C(2*i-1,:)    = sinx(:)';
                C(2*i,:)  = siny(:)';
            end
            fprintf('Observation operator D_N defined with N=%d \n', obj.N);
        end
    end
end
