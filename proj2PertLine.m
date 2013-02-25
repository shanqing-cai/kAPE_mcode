function lambda = proj2PertLine(F1pert, F2pert, F1, F2)
k = F2pert / F1pert;
lambda = (F1 + k * F2) / sqrt(1 + k ^ 2);
return

