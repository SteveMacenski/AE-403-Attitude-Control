syms q0 q0dot q1 q1dot q2 q2dot q3 q3dot
q = [q1;q2;q3];

qdot = [q1dot;q2dot;q3dot];


Q = [-q1, q0, q3, -q2;
     -q2, -q3, q0, q1;
     -q3, q2, -q1, q0;
     q0, q1, q2, q3]
 
 simplify(det(Q))