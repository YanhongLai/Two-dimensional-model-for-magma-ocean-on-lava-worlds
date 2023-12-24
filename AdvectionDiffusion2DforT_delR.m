  
    function [Tt] = AdvectionDiffusion2DforT_delR( nx, nz, dt, dx, delRc, kh, kz, Tx, ut1, wt1, Tt1,tau)
    
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
              
              sx = kh(ii)*dt/dx/dx;   sz = kz(ii)*dt/delRc(i-1)/delRc(i-1);
              cx = ut1(ii)*dt/2/dx;   cz = wt1(ii)*dt/2/delRc(i-1);
                 
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
%          f(ii,1)=Tbot;
         f(ii,1)=Tt1(ii);
     end
        
    T2 = A \ f;
    Tt=T2(:);
    
    end