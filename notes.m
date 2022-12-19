%% Sattellite W Panels
%% Generalized Coordinates

syms q_1(t) q_2(t) t; %%q_1_prime(t) q_2_prime(t) t;
syms E_k(t) E_p(t) D_d(t);

%% Energy Equations

syms kp bp m j0 lp

E_k(t) = .5 * j0 * diff(q_1, t)^2 + m * (diff(q_2, t) + lp * diff(q_1, t))^2
E_p(t) = kp * q_2^2
D_d(t) = bp * diff(q_2, t)^2


%% Euler-Lagrange

L_1 = E_k - E_p

M_1 = [
    subs(diff(E_k(t), diff(q_1(t), t)), [q_2(t), diff(q_1(t), t)], [0, 1]) subs(diff(E_k(t), diff(q_1(t), t)), [q_1(t), diff(q_2(t), t)], [0, 1])
    subs(diff(E_k(t), diff(q_2(t), t)), [q_2(t), diff(q_1(t), t)], [0, 1]) subs(diff(E_k(t), diff(q_2(t), t)), [q_1(t), diff(q_2(t), t)], [0, 1])
]

K_1 = [
    subs(diff(E_p(t), q_1(t)), [q_2(t), q_1(t)], [0 1]) subs(diff(E_p(t), q_1(t)), [q_2(t), q_1(t)], [1 0]);
    subs(diff(E_p(t), q_2(t)), [q_2(t), q_1(t)], [0 1]) subs(diff(E_p(t), q_2(t)), [q_2(t), q_1(t)], [1 0]);
]

C_1 = [
    subs(diff(D_d(t), diff(q_1(t), t)), [diff(q_2(t), t), diff(q_1(t), t)], [0, 1]) subs(diff(D_d(t), diff(q_1(t), t)), [diff(q_2(t), t), diff(q_1(t), t)], [1, 0]);
    subs(diff(D_d(t), diff(q_2(t), t)), [diff(q_2(t), t), diff(q_1(t), t)], [0, 1]) subs(diff(D_d(t), diff(q_2(t), t)), [diff(q_2(t), t), diff(q_1(t), t)], [1, 0]);
]

OUT = [
   "Theta(T)";     % System 1
   "Yp(T)"      % System 2
];

IN = [
   "T_in(t)",     % Input 1
   "NONE"      % Input 2
];

R_1 = [
    1  -1  % System 1 [IN1, IN2]
    0  1; % System 2 [IN1, IN2]
]
M_inv = inv(M_1);

A = [-M_inv * C_1 -M_inv*K_1; [1 0; 0 1] [0 0; 0 0]]
B = [M_inv * R_1; [0 0; 0 0]]
C = [[0 0; 0 0] [1 0; 0 1]]
I = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]

%% Transfer Functions

syms s;

P = C*inv(s*I-A)*B
pretty(simplify(P(1,1)))
pretty(simplify(P(2,1)))

%% Question 4

P_n = subs(P, [m, lp, kp, bp, j0], [20, 2, 320, .5, 580])

% P_n column is determined by input
% P_n row is determined by system output.

for system = 1:2
    for input = 1:2
        disp(OUT(system) + " / " + IN(input) + " P(" + system + ", " + input + ")");
        [num, den] = numden(P_n(system, input));
           
        int = tf(sym2poly(num), sym2poly(den));
        pretty(simplify(P(system, input)))
    
        assignin('base', ['q' num2str(system) 'i' num2str(input)], int);
        
    end
end


