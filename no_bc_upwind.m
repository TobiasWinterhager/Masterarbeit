function [A,xd,yd,hx,hy] = no_bc_upwind(xa,xb,ya,yb,nx,ny,fx,fy)

hx = (xb-xa)/(nx-1); 
hy = (yb-ya)/(ny-1); 
xd = xa:hx:xb;
yd = ya:hy:yb;
ex = ones(nx,1);
ey = ones(ny,1);
Ix = speye(nx);
Iy = speye(ny); 

Grad_x_fw = 1/hx*spdiags([-ex ex],0:1,nx,nx);
Grad_y_fw = 1/hy*spdiags([-ey ey],0:1,ny,ny);
Grad_x_bw = 1/hx*spdiags([-ex ex],-1:0,nx,nx);
Grad_y_bw = 1/hy*spdiags([-ey ey],-1:0,ny,ny);

e1x = Ix(:,1);
enx = Ix(:,nx);
e1y = Iy(:,1);
eny = Iy(:,ny);

O = sparse(nx*ny,nx*ny); 

[X,Y] = meshgrid(xd,yd);
Fx = spdiags(reshape(arrayfun(fx,X',Y')',nx*ny,1),0,nx*ny,nx*ny);
Fy = spdiags(reshape(arrayfun(fy,X',Y')',nx*ny,1),0,nx*ny,nx*ny); 


A =  max(Fx,O)*kron(Grad_x_fw,Iy) + min(Fx,O)*kron(Grad_x_bw,Iy) ...
      + max(Fy,O)*kron(Ix,Grad_y_fw) + min(Fy,O)*kron(Ix,Grad_y_bw);
