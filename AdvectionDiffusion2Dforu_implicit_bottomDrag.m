  
    function [Tt] = AdvectionDiffusion2Dforu_implicit_bottomDrag( nx, nz, dt, dx, dz, x, z, kh, kz, ut1, wt1, pres1t, rho1t, taux, bottomDragLinear)
    
    A = sparse ( nx*nz, nx*nz );
    f = zeros ( nx*nz, 1 );
    Gt1 = zeros(nx*nz, 1);  Gt2 = zeros(nx*nz, 1);
    
    nxz = nx*nz;
    num = 1; Number=zeros(nz,nx);
    for i=1:nz
      for j=1:nx
          Number(i,j)=num;
          num=num+1;
      end
    end
    
    %surface wind forcing
    rhos=zeros(nx,1); rhos(:)=rho1t(1:nx);
    dz_s=z(1)-z(2);
    Fu=taux./rhos;  Fu=Fu./dz_s;
    
%     ubot=zeros(nx,1); ubot(1:end)=ut1(nxz-nx+1:nxz,1);
%     Fb=-bottomDragLinear.*ubot/dz_s;
     
     for i=1:nz
%      for i=1:nz-1
         for j=1:nx
              ii = Number(i,j);
              
              if (i>2) && (i<nz-1)
                   ii2 = Number(i+1,j);
                   ii1 = Number(i-1,j);
                   
                   sz2 = kz(ii2)*dt/(z(i+1)-z(i-1))/(z(i+2)-z(i));   
                   sz1 = kz(ii1)*dt/(z(i+1)-z(i-1))/(z(i)-z(i-2));    
                   sz = sz2+sz1; 
                   
                   if j==1
                       ii0=Number(i,nx-1);
                       sx2 = kh(ii+1)*dt/4/dx/dx;  sx1 = kh(ii0)*dt/4/dx/dx;
                       A(ii, Number(i,nx-2)) = -sx1;  A(ii, Number(i,j+2)) = -sx2;
                   elseif j==2
                       sx2 = kh(ii+1)*dt/4/dx/dx;  sx1 = kh(ii-1)*dt/4/dx/dx;
                       A(ii, Number(i,nx-1)) = -sx1;  A(ii, Number(i,j+2)) = -sx2;
                   elseif j==nx-1
                       sx2 = kh(ii+1)*dt/4/dx/dx;  sx1 = kh(ii-1)*dt/4/dx/dx;
                       A(ii, Number(i,j-2)) = -sx1; A(ii, Number(i,2)) = -sx2;
                   elseif j==nx
                       ii0=Number(i,2);
                       sx2 = kh(ii0)*dt/4/dx/dx;  sx1 = kh(ii-1)*dt/4/dx/dx;
                       A(ii, Number(i,j-2)) = -sx1; A(ii, Number(i,3)) = -sx2;
                   else
                       sx2 = kh(ii+1)*dt/4/dx/dx;  sx1 = kh(ii-1)*dt/4/dx/dx;   
                       A(ii, Number(i,j-2)) = -sx1;  A(ii, Number(i,j+2)) = -sx2;
                   end
                   
                    sx = sx2+sx1;
                    A(ii,ii) = 1+sx+sz;
                    A(ii, Number(i-2,j)) = -sz1;
                    A(ii, Number(i+2,j)) = -sz2;
              else
                  if i==1
                      sx = 2*kh(ii)*dt/dx/dx;     sz = kz(ii)*dt/dz/dz;
                      A(ii,ii) = 1+sx+sz;
                      A(ii, Number(i+1,j)) = -kz(ii)*dt/dz/dz;
                  elseif i==nz
                      sx = 2*kh(ii)*dt/dx/dx;     sz = kz(ii)*dt/dz/dz;
                      A(ii,ii) = 1+sx+sz+bottomDragLinear*dt/dz_s;
                      A(ii, Number(i-1,j)) = -kz(ii)*dt/dz/dz;
                  else
                    sx = 2*kh(ii)*dt/dx/dx;     sz = 2*kz(ii)*dt/dz/dz;
                    A(ii,ii) = 1+sx+sz;
                    A(ii, Number(i-1,j)) = -kz(ii)*dt/dz/dz;
                    A(ii, Number(i+1,j)) = -kz(ii)*dt/dz/dz;
                  end
                  
                    if j==1
                        A(ii, Number(i,nx-1)) = -kh(ii)*dt/dx/dx;
                        A(ii, Number(i,j+1)) = -kh(ii)*dt/dx/dx;
                    elseif j==nx
                        A(ii, Number(i,j-1)) = -kh(ii)*dt/dx/dx;
                        A(ii, Number(i,2)) = -kh(ii)*dt/dx/dx;
                    else
                        A(ii, Number(i,j-1)) = -kh(ii)*dt/dx/dx;
                        A(ii, Number(i,j+1)) = -kh(ii)*dt/dx/dx;
                    end
              end
              
              %%%%%%%%%%%%%% other terms
              if j==1
                  ii0=Number(i,nx-1);
                  Gt1(ii,1) = -ut1(ii)*(ut1(ii+1)-ut1(ii0))/2/dx -(pres1t(ii+1)-pres1t(ii0))/rho1t(ii)/2/dx;
              elseif j==nx
                  ii0=Number(i,2);
                  Gt1(ii,1) = -ut1(ii)*(ut1(ii0)-ut1(ii-1))/2/dx -(pres1t(ii0)-pres1t(ii-1))/rho1t(ii)/2/dx;
              else    
                  Gt1(ii,1) = -ut1(ii)*(ut1(ii+1)-ut1(ii-1))/2/dx -(pres1t(ii+1)-pres1t(ii-1))/rho1t(ii)/2/dx;
              end
              
              if i==1
                  ii2 = Number(i+1,j);
                  Gt2(ii,1) =  -wt1(ii)*(ut1(ii2)-ut1(ii))/dz + Fu(j);   %surface wind forcing
              elseif i==nz
                  ii1 = Number(i-1,j);
                  Gt2(ii,1) =  -wt1(ii)*(ut1(ii)-ut1(ii1))/dz;   %bottom linear drag
              else
                  ii2 = Number(i+1,j);
                  ii1 = Number(i-1,j);
                  Gt2(ii,1) =  -wt1(ii)*(ut1(ii2)-ut1(ii1))/2/dz;
              end
              
              f(ii,1) = ut1(ii) + dt*Gt1(ii,1) + dt*Gt2(ii,1);
         end
     end
     
%       for j=1:nx
%          ii=Number(nz,j);
%          A(ii,ii)=1;
%          f(ii,1)=ut1(ii);
%       end
     
    T2 = A \ f;
    Tt=T2(:);
    
%     for j=1:nx
%         ii=Number(nz,j);
%         Tt(ii)=0;   %no-slip bottom boundary
%     end
    
    
    end