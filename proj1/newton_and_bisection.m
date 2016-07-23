%newton method
function [fun, dfun, x, out] = NewtonMethod(Fun, dFun,x0, params)
[FunFcn,msg] = fcnchk(Fun,0);
if ~isempty(msg)
    error('InvalidFUN',msg);
end
[dFunFcn,msg] = fcnchk(dFun,0);
if ~isempty(msg)
    error('InvalidFUN',msg);
end
out.flg = 1;
x(1)= x0;

N = params.MaxIt;
tol = params.tol;
out.x = [];
out.f =[];

for k = 1:N
    fun(k)  = FunFcn(x(k));
    dfun(k) = dFunFcn(x(k));
    out.x = [out.x;x(k)];
    out.f =[out.f;fun(k)];
    if (abs(fun(k)) < tol)
       out.flg = 0;
       out.it = k;
       return;
    end
    if (dfun(k) == 0)
       out.it = k;
       return;
    end
    x(k+1) = x(k) - fun(k)/dfun(k);
end

       

