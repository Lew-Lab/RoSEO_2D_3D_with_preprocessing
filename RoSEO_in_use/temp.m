syms x F

y = ifft(F*fft(x));
gradient(y,x)

%%
syms x
eqn = (-16*x^2+2*x-1/16)^2/(x-1/16)^2==1;
S = solve(eqn)