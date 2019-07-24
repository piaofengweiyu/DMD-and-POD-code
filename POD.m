function [phi,energy,Time] = POD(data)
   
  r2 = size(data,2);  
    C=data'*data/r2;
    
    [beta, lmd] = svd(C,'econ');
    phi = data*beta;
    for i=1:r2
        phi(:,i) = phi(:,i)/norm(phi(:,i),2);
        
    end
    
    Time= data'*phi;
    energy=zeros(r2,1);
    for i=1:r2
        energy(i,1) = norm(Time(:,i),2);
        
    end
    
    
    
end