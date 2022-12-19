
%%
plant = tf(1, [1, 15, 50, 0])


%% Designing a PD controller

syms kp td n s a out;

vars = [kp n td];
vals = [200 2 1];

a = kp * (1 + (td * s)/((td/n) * s + 1))

[num, den] = numden(a);
sol = [num; den];

zeros = [];
poles = [];
nogain = 1;

for i = 1:2
    vns = solve(sol(i)==0, s);
    for j = 1:size(vns, 2)
        if i == 1
            nogain = nogain * (s - vns(j));
            zeros(j) = subs(vns(j), vars, vals);
        else
            nogain = nogain / (s - vns(j));
            poles(j) = subs(vns(j), vars, vals);
        end
    end
end

cmul = simplify((num / den) / nogain)

symbolicController = cmul * nogain;
pretty(nogain)
pretty(symbolicController)

controller = zpk(zeros, poles, double(subs(cmul, vars, vals)))

step(feedback(series(controller, plant), 1));
