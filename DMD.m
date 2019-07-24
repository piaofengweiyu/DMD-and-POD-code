function [Phi,mu,lambda,P,f,time,delat,Xhat] = DMD(Xori,r,dt)                           
%Phi 模态
% mu 
%f 频率
%Xhat 重构
%time 时间系数
%P 模态能量系数
% lambda 
                                                                                
%% DMD
X = Xori(:, 1:end-1);
Y= Xori(:, 2:end);

[U, S, V] = svd(X, 'econ');
r = min(r, size(U,2));   % r选择可以根据能量大小，选择能量多的几个
% diagS = diag(S);



U_r = U(:, 1:r); % truncate to rank-r
S_r = S(1:r, 1:r);
V_r = V(:, 1:r);

Atilde = U_r' * Y * V_r / S_r;

Ahat = (S_r^(-1/2)) * Atilde * (S_r^(1/2));
[What, D] = eig(Ahat);
W_r = S_r^(1/2) * What;

Phi = Y*V_r/S_r*W_r;

lambda = diag(D);  % discrete-time eigenvalues   
mu = log(lambda)/dt; % continuous-time eigenvalues  

%% 重构
x0 = X(:,1);

b = Phi\x0; % initial conditions


mm1 = size(X, 2); % mm1 = m - 1                                                
time= zeros(r, mm1);                                                  
t = (0:mm1-1)*dt; % time vector                                                 
for iter = 1:mm1                                                                
    time(:,iter) = (b.*exp(mu*t(iter)));                            
end         

% for iter = 1:mm1                                                                
%     time(:,iter) = (b.*lambda.^iter);                            
% end  


Xhat = Phi * time;
%% delat 

% for iter = 1:mm1                                                                
    delat = sqrt(norm((abs(imag(lambda(:)))-1),2)/r);    
    
% end   


%% f, the frequencies of the modes in cycles/sec

f = abs(imag(mu(:))/2/pi); % frequency in cycles/sec

%% P, the power of the modes
P = (diag(Phi'*Phi)); % roughly scales like the fft spectrum

