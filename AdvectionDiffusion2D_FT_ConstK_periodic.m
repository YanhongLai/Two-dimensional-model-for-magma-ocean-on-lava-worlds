  
    function [Tt] = AdvectionDiffusion2D_FT_ConstK_periodic( nx, nz, dt, dx, dz, x, z, kh, kz, Tx, Tbot, Twest, Teast, ut1, wt1, Tt1,tau)
    
    A = sparse ( nx*nz, nx*nz );
    f = zeros ( nx*nz, 1 );
      
    num = 1; Number=zeros(nz,nx);
    for i=1:nz
      for j=1:nx
          Number(i,j)=num;
          num=num+1;
      end
    end
   
     for i=2:nz-1
         for j=1:nx
             
              ii = Number(i,j);
              
%               ii2 = Number(i+1,j);
%               ii1 = Number(i-1,j);              
%               sz2 = kz(ii2)*dt/(z(i+1)-z(i-1))/(z(i+2)-z(i));   sz1 = kz(ii1)*dt/(z(i+1)-z(i-1))/(z(i)-z(i-2));    sz = sz2+sz1;
%               cz = wt1(ii,1)*dt/(z(i+1)-z(i-1));
%               sx2 = kh(ii+1)*dt/4/dx/dx;  sx1 = kh(ii-1)*dt/4/dx/dx;   sx = sx2+sx1;
%               cx = ut1(ii,1)*dt/2/dx;  
              
              sx = kh(ii)*dt/dx/dx;   sz = kz(ii)*dt/dz/dz;
              cx = ut1(ii)*dt/2/dx;   cz = wt1(ii)*dt/2/dz;
                 
              A(ii,ii) = 1+2*sx+2*sz;
              A(ii, Number(i-1,j)) = -cz-sz;
              A(ii, Number(i+1,j)) =  cz-sz;
              
              if j==1
                  A(ii, Number(i,nx-1)) = -cx-sx;
                  A(ii, Number(i,j+1)) =  cx-sx;
              elseif j==nx
                  A(ii, Number(i,j-1)) = -cx-sx;
                  A(ii, Number(i,2)) =  cx-sx;
              else
                  A(ii, Number(i,j-1)) = -cx-sx;
                  A(ii, Number(i,j+1)) =  cx-sx;
              end
                      
              f(ii,1) = Tt1(ii); 
         end
     end
     
     %fixed boundary conditions
       for j=1:nx 
         ii=Number(1,j);
         A(ii,ii)=1;
%          f(ii,1)=Tx(j);
        f(ii,1)=Tt1(ii)-(Tt1(ii)-Tx(j))/tau*dt;
       end         

     for j=1:nx
         ii=Number(nz,j);
         A(ii,ii)=1;
         f(ii,1)=Tbot;
%          f(ii,1)=Tt1(ii);
     end
        
    T2 = A \ f;
    Tt=T2(:);
    
    end