% Integral of Signum function approximation

function y = int_sgn_approx(x)
    y = x*atan(100*x)-(ln(10000*x^2+1)/200);
end