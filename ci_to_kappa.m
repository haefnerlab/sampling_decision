function kappa = ci_to_kappa(ci)
os = linspace(0,pi,1001);
os = os(1:end-1);
ks = linspace(0, 6);
for ik=length(ks):-1:1
    p1 = exp(ks(ik)*cos(os*2)); p1 = p1 ./ sum(p1);
    p2 = exp(-ks(ik)*cos(os*2)); p2 = p2 ./ sum(p2);
    ratio = p1 ./ (p1 + p2);
    emp_ci(ik) = dot(p1, ratio);
end

ks = [0 ks 7];
emp_ci = [.5 emp_ci 1];

kappa = interp1(emp_ci, ks, ci);
end