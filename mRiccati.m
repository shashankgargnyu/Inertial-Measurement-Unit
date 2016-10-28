function dXdt = mRiccati(t, X, A, B, Q)
X = reshape(X, size(A)); %Convert from "n^2"-by-1 to "n"-by-"n"
dXdt = A*X + X*A' - X*B*X + Q; %Determine derivative
dXdt = dXdt(:); %Convert from "n"-by-"n" to "n^2"-by-1