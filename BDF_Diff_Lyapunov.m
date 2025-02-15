function Y_m = BDF_Diff_Lyapunov(A,C,Y0,t0,tf, N)
    % [t,Y,Y_m] = BDF_Diff_Lyapunov(A,B,C,Y0,t0,tf)
    % The resolution of the differential Lyapunov matrix equations
    %  ( backward differentiation formula method)
    %                                  dY/dt=AY+YA'+C (*)
    %                                    Y(T0) = Y0
    % Input:        A          : size matrix  (m,m).
    %               C          : size matrix  (m,m).
    %               Y0         : size matrix  (m,m).
    %               tf , t0    : temps.
    % Output:    
    %             Y_m          : solution to the inslant Tf of (*)  
    %             Y            : solution to (*) 
    %             t
    %    Author                :  LAKHLIFA SADEK.
    %    Last modification     :  19/06/2018
    % E-mail: lakhlifasdek@gmail.com; sadek.l@ucd.ac.ma
    % ORCID : https://orcid.org/0000-0001-9780-2592
    % 
    fprintf('Solving Differential Ricatti equation with %d time steps.\n',N)
    h=(tf-t0)/N; 
    [ms]=size(A,1);
    Z0=Y0;Z1=Y0;
    TAA=(2/3)*h*A-(1/2)*eye(ms); 
      alpha0=4/3;alpha1=-1/3;
    %t(1)=t0;
    %Y=[Y0(:)];
    for k=1:N
        fprintf('Step %d of %d \n',k, N)
        %t(k+1)=t0+(k)*h;
        CC=alpha0*Z0+alpha1*Z1+h*(2/3)*C;
        Yk=lyap(TAA,CC);
       Z1=Z0;
       Z0=Yk;  
       %Y=[Y Yk(:)];
    end
    Y_m=Yk;
    fprintf('Differential Ricatti equation solved.\n')
    end
