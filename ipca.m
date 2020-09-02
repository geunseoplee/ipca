function [U, D, V] = ipca(B, b, newb, U, D, V)
% Geunseop Lee. 2020.09.01
% example [U, D] = ipca(B, b, newb, U, D); or [U, D, V] = ipca(B, b, newb,
% U, D, V];
% U, D, V : svd of zero mean matrix
% b : mean vector of previous data matrix
% newb : mean vector after updating data matrix
% B : newly added data appended to data matrix

    p = size(B,2);
    k = size(D,1);
    S = [B - b * ones(1,p), b - newb];
    UtS = U' * S;
    
    [Q, R] = qr(S - U * UtS,0);
    U1 = [U,Q];
    if nargout > 2
         n = size(V,1);
         [Z, W] = qr([zeros(n,p), ones(n,1)],0);
    else
         W = zeros(p+1, p+1);
         W(:,end) = 1;
    end
     T1 = [eye(p,p), ones(p,1)];
    
     Lambda = [D, UtS * W', UtS * T1'; zeros(size(U1,2) - tr, tr), R * W', R * T1'];
     
     [Un, Dn, Vn] = svd(Lambda,'econ');
     
      D = Dn(1:k, 1:k);
      U = U1 * Un(:,1:k);
      
      if nargout > 2
          V1 = [V, Z, zeros(size(V,1), size(T1,1));zeros(size(T1,1), ...
          size(V,2)), zeros(size(T1,1), size(Z,2)), eye(size(T1,1), size(T1,1))];
          V = V1 * Vn(:,1:k);
      end
end
