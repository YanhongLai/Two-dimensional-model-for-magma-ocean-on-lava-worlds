function [T2]=conversion_theta2temp_delR(nx, nz, delRc, alpha, cp, g, theta)
    
  num = 1; Number=zeros(nz,nx);
  for i=1:nz
      for j=1:nx
          Number(i,j)=num;
          num=num+1;
      end
  end

    T2=zeros(nz*nx,1);
    
    for i=1:nz
        for j=1:nx
            ii=Number(i,j);
            if i==1
                T2(ii)=T2(ii);
            else
                for n=2:i
                    ii2=Number(n-1,j);
                    T2(ii)=T2(ii)-theta(ii2)*alpha*g/cp*delRc(n-1);
                end
            end
            T2(ii)=theta(ii)+T2(ii);
        end
    end
         
    
    end