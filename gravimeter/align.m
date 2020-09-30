function [A] = align(V)
%A!ALIGN(V) returns a constat matrix A which is the real alignment of the
%INVERSE of the complex input matrix V
%from Mohit slides

    if (nargin ==0) || (nargin > 1)
        disp('usage: mat_inv_real = align(mat)')
        return
    end

    D = pinv(real(V'*V));
    A = D*real(V'*diag(exp(1i * angle(diag(V*D*V.'))/2)));


end
