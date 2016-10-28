function dxdt=f1(t,A,B,C,Gain,U,X_new,y)
dxdt = A*X_new + B*U + Gain*(y-C*X_new);