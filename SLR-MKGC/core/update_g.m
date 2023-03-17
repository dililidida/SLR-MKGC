function [ g,fval,exitflag ] = update_g(L, H, Lrank)
view_num = size(Lrank,2);
a = zeros(view_num,1);
for i = 1:view_num
    a(i) = -trace(L*Lrank{i});
end
A = [];
b = [];
Aeq = ones(1,view_num);
beq = 1;
lb = zeros(view_num,1);
ub = ones(view_num,1);
g = ones(view_num, 1)/view_num;
options.Display = 'off';
%options.Display = 'on';
[g,fval,exitflag] = quadprog(H,a,A,b,Aeq,beq,lb,ub,g,options);
g = g./sum(g);
end

