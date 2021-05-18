function [output] = modelEulerM2_v3(stateVector,parameterValues, digraph, mode, par1, par2)
% Model file written without using for loops
% Written using Euler's method do reduce the use of ODE15s
% (ex) g1_1 represent gene1(target) in Cell #1
%      g2_30 represent gene2(target) in Cell #30
%  M2 condition and each genes have different degradation values

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        STATES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

g1 = stateVector(1:3:96);
g2 = stateVector(2:3:96);
g3 = stateVector(3:3:96);
% g1_1 = stateVector(1);
% g2_1 = stateVector(2);
% g3_1 = stateVector(3);
% g1_2 = stateVector(4);
% g2_2 = stateVector(5);
% g3_2 = stateVector(6);
% g1_3 = stateVector(7);
% g2_3 = stateVector(8);
% g3_3 = stateVector(9);
% g1_4 = stateVector(10);
% g2_4 = stateVector(11);
% g3_4 = stateVector(12);
% g1_5 = stateVector(13);
% g2_5 = stateVector(14);
% g3_5 = stateVector(15);
% g1_6 = stateVector(16);
% g2_6 = stateVector(17);
% g3_6 = stateVector(18);
% g1_7 = stateVector(19);
% g2_7 = stateVector(20);
% g3_7 = stateVector(21);
% g1_8 = stateVector(22);
% g2_8 = stateVector(23);
% g3_8 = stateVector(24);
% g1_9 = stateVector(25);
% g2_9 = stateVector(26);
% g3_9 = stateVector(27);
% g1_10 = stateVector(28);
% g2_10 = stateVector(29);
% g3_10 = stateVector(30);
% g1_11 = stateVector(31);
% g2_11 = stateVector(32);
% g3_11 = stateVector(33);
% g1_12 = stateVector(34);
% g2_12 = stateVector(35);
% g3_12 = stateVector(36);
% g1_13 = stateVector(37);
% g2_13 = stateVector(38);
% g3_13 = stateVector(39);
% g1_14 = stateVector(40);
% g2_14 = stateVector(41);
% g3_14 = stateVector(42);
% g1_15 = stateVector(43);
% g2_15 = stateVector(44);
% g3_15 = stateVector(45);
% g1_16 = stateVector(46);
% g2_16 = stateVector(47);
% g3_16 = stateVector(48);
% g1_17 = stateVector(49);
% g2_17 = stateVector(50);
% g3_17 = stateVector(51);
% g1_18 = stateVector(52);
% g2_18 = stateVector(53);
% g3_18 = stateVector(54);
% g1_19 = stateVector(55);
% g2_19 = stateVector(56);
% g3_19 = stateVector(57);
% g1_20 = stateVector(58);
% g2_20 = stateVector(59);
% g3_20 = stateVector(60);
% g1_21 = stateVector(61);
% g2_21 = stateVector(62);
% g3_21 = stateVector(63);
% g1_22 = stateVector(64);
% g2_22 = stateVector(65);
% g3_22 = stateVector(66);
% g1_23 = stateVector(67);
% g2_23 = stateVector(68);
% g3_23 = stateVector(69);
% g1_24 = stateVector(70);
% g2_24 = stateVector(71);
% g3_24 = stateVector(72);
% g1_25 = stateVector(73);
% g2_25 = stateVector(74);
% g3_25 = stateVector(75);
% g1_26 = stateVector(76);
% g2_26 = stateVector(77);
% g3_26 = stateVector(78);
% g1_27 = stateVector(79);
% g2_27 = stateVector(80);
% g3_27 = stateVector(81);
% g1_28 = stateVector(82);
% g2_28 = stateVector(83);
% g3_28 = stateVector(84);
% g1_29 = stateVector(85);
% g2_29 = stateVector(86);
% g3_29 = stateVector(87);
% g1_30 = stateVector(88);
% g2_30 = stateVector(89);
% g3_30 = stateVector(90);
% g1_31 = stateVector(91);
% g2_31 = stateVector(92);
% g3_31 = stateVector(93);
% g1_32 = stateVector(94);
% g2_32 = stateVector(95);
% g3_32 = stateVector(96);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s11 = parameterValues(1);
s12 = parameterValues(2);
s13 = parameterValues(3);
s21 = parameterValues(4);
s22 = parameterValues(5);
s23 = parameterValues(6);
s31 = parameterValues(7);
s32 = parameterValues(8);
s33 = parameterValues(9);
dRg1 = parameterValues(10); %degradation rate
dRg2 = parameterValues(11);
dRg3 = parameterValues(12);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       DIFFERENTIAL EQUATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% VERSION 1 of phi

if mode == 1
    phi = @(a)1./(1+exp(par1-par2.*a));
else
    
    phi = @(a)1./(1+ par1 * a);
end
% phi = @(a)1./(1+exp(param1-param2.*a));

% VERSION 2 of phi
% param3 = 1
% phi = @(a)1./(1+ param3 * a);


idx1 = findedge(digraph,[1,2,3],[1,1,1]); % index of all incoming edges for gene 1
idx2 = findedge(digraph,[1,2,3],[2,2,2]); % index of all incoming edges for gene 2
idx3 = findedge(digraph,[1,2,3],[3,3,3]); % index of all incoming edges for gene 3
a1 = digraph.Edges.Weight(idx1(1));
a2 = digraph.Edges.Weight(idx1(2));
a3 = digraph.Edges.Weight(idx1(3));
b1 = digraph.Edges.Weight(idx2(1));
b2 = digraph.Edges.Weight(idx2(2));
b3 = digraph.Edges.Weight(idx2(3));
c1 = digraph.Edges.Weight(idx3(1));
c2 = digraph.Edges.Weight(idx3(2));
c3 = digraph.Edges.Weight(idx3(3));


M_p2 = @(x)1./(1+exp((-3.9.*(55.5-x)+186)./8));

if findedge(digraph, 4, 1)== 10 % Find Morphogen connection
    M1 = 1;
    M2 = 0;
    M3 = 0;
elseif findedge(digraph, 4, 2)== 10
    M1 = 0;
    M2 = 1;
    M3 = 0;        
else
    M1 = 0;
    M2 = 0;
    M3 = 1;
end

stepsize = 0.05;
simtime = 10; % 10 minutes simulation time
numiter = simtime/stepsize;

G1 = g1;
G2 = g2;
G3 = g3;
% This is to add the stopping condition
G1_start = zeros(len(G1))
G2_start = zeros(len(G2))
G3_start = zeros(len(G3))

for n = 1:numiter-1  
    
    g1_param = a1.*s11.*G1 + a2.*s21.*G2 + a3.*s31.*G3 +M_p2(1:32).*M1;
    g1_dot = phi(g1_param)-dRg1.*G1;
    G1 = G1+ stepsize*g1_dot;
    
    g2_param = b1.*s12.*G1 + b2.*s22.*G2 + b3.*s32.*G3 +M_p2(1:32).*M2;
    g2_dot = phi(g2_param)-dRg2.*G2;
    G2 = G2+ stepsize*g2_dot;
    
    g3_param = c1.*s13.*G1 + c2.*s23.*G2 + c3.*s33.*G3 +M_p2(1:32).*M3;
    g3_dot = phi(g3_param)-dRg3.*G3;
    G3 = G3+ stepsize*g3_dot;
    if mod(n,50) == 0 
    % this is to add the stopping condition
        diff = sum(abs(G1-G1_start) + abs(G2-G2_start) + abs(G3-G3_start))/3/32;
        if diff<2
            break
        G1_start = G1
        G2_start = G2
        G3_start = G3
%     g1_1_dot = phi(a1*s11*g1_1(n) + a2*s21*g2_1(n) + a3*s31*g3_1(n) + M_p2(1)*M1)-dRg1*g1_1(n);
%     g2_1_dot = phi(b1*s12*g1_1(n) + b2*s22*g2_1(n) + b3*s32*g3_1(n) + M_p2(1)*M2)-dRg2*g2_1(n);
%     g3_1_dot = phi(c1*s13*g1_1(n) + c2*s23*g2_1(n) + c3*s33*g3_1(n) + M_p2(1)*M3)-dRg3*g3_1(n);
%     g1_2_dot = phi(a1*s11*g1_2(n) + a2*s21*g2_2(n) + a3*s31*g3_2(n) + M_p2(2)*M1)-dRg1*g1_2(n);
%     g2_2_dot = phi(b1*s12*g1_2(n) + b2*s22*g2_2(n) + b3*s32*g3_2(n) + M_p2(2)*M2)-dRg2*g2_2(n);
%     g3_2_dot = phi(c1*s13*g1_2(n) + c2*s23*g2_2(n) + c3*s33*g3_2(n) + M_p2(2)*M3)-dRg3*g3_2(n);
%     g1_3_dot = phi(a1*s11*g1_3(n) + a2*s21*g2_3(n) + a3*s31*g3_3(n) + M_p2(3)*M1)-dRg1*g1_3(n);
%     g2_3_dot = phi(b1*s12*g1_3(n) + b2*s22*g2_3(n) + b3*s32*g3_3(n) + M_p2(3)*M2)-dRg2*g2_3(n);
%     g3_3_dot = phi(c1*s13*g1_3(n) + c2*s23*g2_3(n) + c3*s33*g3_3(n) + M_p2(3)*M3)-dRg3*g3_3(n);
%     g1_4_dot = phi(a1*s11*g1_4(n) + a2*s21*g2_4(n) + a3*s31*g3_4(n) + M_p2(4)*M1)-dRg1*g1_4(n);
%     g2_4_dot = phi(b1*s12*g1_4(n) + b2*s22*g2_4(n) + b3*s32*g3_4(n) + M_p2(4)*M2)-dRg2*g2_4(n);
%     g3_4_dot = phi(c1*s13*g1_4(n) + c2*s23*g2_4(n) + c3*s33*g3_4(n) + M_p2(4)*M3)-dRg3*g3_4(n);
%     g1_5_dot = phi(a1*s11*g1_5(n) + a2*s21*g2_5(n) + a3*s31*g3_5(n) + M_p2(5)*M1)-dRg1*g1_5(n);
%     g2_5_dot = phi(b1*s12*g1_5(n) + b2*s22*g2_5(n) + b3*s32*g3_5(n) + M_p2(5)*M2)-dRg2*g2_5(n);
%     g3_5_dot = phi(c1*s13*g1_5(n) + c2*s23*g2_5(n) + c3*s33*g3_5(n) + M_p2(5)*M3)-dRg3*g3_5(n);
%     g1_6_dot = phi(a1*s11*g1_6(n) + a2*s21*g2_6(n) + a3*s31*g3_6(n) + M_p2(6)*M1)-dRg1*g1_6(n);
%     g2_6_dot = phi(b1*s12*g1_6(n) + b2*s22*g2_6(n) + b3*s32*g3_6(n) + M_p2(6)*M2)-dRg2*g2_6(n);
%     g3_6_dot = phi(c1*s13*g1_6(n) + c2*s23*g2_6(n) + c3*s33*g3_6(n) + M_p2(6)*M3)-dRg3*g3_6(n);
%     g1_7_dot = phi(a1*s11*g1_7(n) + a2*s21*g2_7(n) + a3*s31*g3_7(n) + M_p2(7)*M1)-dRg1*g1_7(n);
%     g2_7_dot = phi(b1*s12*g1_7(n) + b2*s22*g2_7(n) + b3*s32*g3_7(n) + M_p2(7)*M2)-dRg2*g2_7(n);
%     g3_7_dot = phi(c1*s13*g1_7(n) + c2*s23*g2_7(n) + c3*s33*g3_7(n) + M_p2(7)*M3)-dRg3*g3_7(n);
%     g1_8_dot = phi(a1*s11*g1_8(n) + a2*s21*g2_8(n) + a3*s31*g3_8(n) + M_p2(8)*M1)-dRg1*g1_8(n);
%     g2_8_dot = phi(b1*s12*g1_8(n) + b2*s22*g2_8(n) + b3*s32*g3_8(n) + M_p2(8)*M2)-dRg2*g2_8(n);
%     g3_8_dot = phi(c1*s13*g1_8(n) + c2*s23*g2_8(n) + c3*s33*g3_8(n) + M_p2(8)*M3)-dRg3*g3_8(n);
%     g1_9_dot = phi(a1*s11*g1_9(n) + a2*s21*g2_9(n) + a3*s31*g3_9(n) + M_p2(9)*M1)-dRg1*g1_9(n);
%     g2_9_dot = phi(b1*s12*g1_9(n) + b2*s22*g2_9(n) + b3*s32*g3_9(n) + M_p2(9)*M2)-dRg2*g2_9(n);
%     g3_9_dot = phi(c1*s13*g1_9(n) + c2*s23*g2_9(n) + c3*s33*g3_9(n) + M_p2(9)*M3)-dRg3*g3_9(n);
%     g1_10_dot = phi(a1*s11*g1_10(n) + a2*s21*g2_10(n) + a3*s31*g3_10(n) + M_p2(10)*M1)-dRg1*g1_10(n);
%     g2_10_dot = phi(b1*s12*g1_10(n) + b2*s22*g2_10(n) + b3*s32*g3_10(n) + M_p2(10)*M2)-dRg2*g2_10(n);
%     g3_10_dot = phi(c1*s13*g1_10(n) + c2*s23*g2_10(n) + c3*s33*g3_10(n) + M_p2(10)*M3)-dRg3*g3_10(n);
%     g1_11_dot = phi(a1*s11*g1_11(n) + a2*s21*g2_11(n) + a3*s31*g3_11(n) + M_p2(11)*M1)-dRg1*g1_11(n);
%     g2_11_dot = phi(b1*s12*g1_11(n) + b2*s22*g2_11(n) + b3*s32*g3_11(n) + M_p2(11)*M2)-dRg2*g2_11(n);
%     g3_11_dot = phi(c1*s13*g1_11(n) + c2*s23*g2_11(n) + c3*s33*g3_11(n) + M_p2(11)*M3)-dRg3*g3_11(n);
%     g1_12_dot = phi(a1*s11*g1_12(n) + a2*s21*g2_12(n) + a3*s31*g3_12(n) + M_p2(12)*M1)-dRg1*g1_12(n);
%     g2_12_dot = phi(b1*s12*g1_12(n) + b2*s22*g2_12(n) + b3*s32*g3_12(n) + M_p2(12)*M2)-dRg2*g2_12(n);
%     g3_12_dot = phi(c1*s13*g1_12(n) + c2*s23*g2_12(n) + c3*s33*g3_12(n) + M_p2(12)*M3)-dRg3*g3_12(n);
%     g1_13_dot = phi(a1*s11*g1_13(n) + a2*s21*g2_13(n) + a3*s31*g3_13(n) + M_p2(13)*M1)-dRg1*g1_13(n);
%     g2_13_dot = phi(b1*s12*g1_13(n) + b2*s22*g2_13(n) + b3*s32*g3_13(n) + M_p2(13)*M2)-dRg2*g2_13(n);
%     g3_13_dot = phi(c1*s13*g1_13(n) + c2*s23*g2_13(n) + c3*s33*g3_13(n) + M_p2(13)*M3)-dRg3*g3_13(n);
%     g1_14_dot = phi(a1*s11*g1_14(n) + a2*s21*g2_14(n) + a3*s31*g3_14(n) + M_p2(14)*M1)-dRg1*g1_14(n);
%     g2_14_dot = phi(b1*s12*g1_14(n) + b2*s22*g2_14(n) + b3*s32*g3_14(n) + M_p2(14)*M2)-dRg2*g2_14(n);
%     g3_14_dot = phi(c1*s13*g1_14(n) + c2*s23*g2_14(n) + c3*s33*g3_14(n) + M_p2(14)*M3)-dRg3*g3_14(n);
%     g1_15_dot = phi(a1*s11*g1_15(n) + a2*s21*g2_15(n) + a3*s31*g3_15(n) + M_p2(15)*M1)-dRg1*g1_15(n);
%     g2_15_dot = phi(b1*s12*g1_15(n) + b2*s22*g2_15(n) + b3*s32*g3_15(n) + M_p2(15)*M2)-dRg2*g2_15(n);
%     g3_15_dot = phi(c1*s13*g1_15(n) + c2*s23*g2_15(n) + c3*s33*g3_15(n) + M_p2(15)*M3)-dRg3*g3_15(n);
%     g1_16_dot = phi(a1*s11*g1_16(n) + a2*s21*g2_16(n) + a3*s31*g3_16(n) + M_p2(16)*M1)-dRg1*g1_16(n);
%     g2_16_dot = phi(b1*s12*g1_16(n) + b2*s22*g2_16(n) + b3*s32*g3_16(n) + M_p2(16)*M2)-dRg2*g2_16(n);
%     g3_16_dot = phi(c1*s13*g1_16(n) + c2*s23*g2_16(n) + c3*s33*g3_16(n) + M_p2(16)*M3)-dRg3*g3_16(n);
%     g1_17_dot = phi(a1*s11*g1_17(n) + a2*s21*g2_17(n) + a3*s31*g3_17(n) + M_p2(17)*M1)-dRg1*g1_17(n);
%     g2_17_dot = phi(b1*s12*g1_17(n) + b2*s22*g2_17(n) + b3*s32*g3_17(n) + M_p2(17)*M2)-dRg2*g2_17(n);
%     g3_17_dot = phi(c1*s13*g1_17(n) + c2*s23*g2_17(n) + c3*s33*g3_17(n) + M_p2(17)*M3)-dRg3*g3_17(n);
%     g1_18_dot = phi(a1*s11*g1_18(n) + a2*s21*g2_18(n) + a3*s31*g3_18(n) + M_p2(18)*M1)-dRg1*g1_18(n);
%     g2_18_dot = phi(b1*s12*g1_18(n) + b2*s22*g2_18(n) + b3*s32*g3_18(n) + M_p2(18)*M2)-dRg2*g2_18(n);
%     g3_18_dot = phi(c1*s13*g1_18(n) + c2*s23*g2_18(n) + c3*s33*g3_18(n) + M_p2(18)*M3)-dRg3*g3_18(n);
%     g1_19_dot = phi(a1*s11*g1_19(n) + a2*s21*g2_19(n) + a3*s31*g3_19(n) + M_p2(19)*M1)-dRg1*g1_19(n);
%     g2_19_dot = phi(b1*s12*g1_19(n) + b2*s22*g2_19(n) + b3*s32*g3_19(n) + M_p2(19)*M2)-dRg2*g2_19(n);
%     g3_19_dot = phi(c1*s13*g1_19(n) + c2*s23*g2_19(n) + c3*s33*g3_19(n) + M_p2(19)*M3)-dRg3*g3_19(n);
%     g1_20_dot = phi(a1*s11*g1_20(n) + a2*s21*g2_20(n) + a3*s31*g3_20(n) + M_p2(20)*M1)-dRg1*g1_20(n);
%     g2_20_dot = phi(b1*s12*g1_20(n) + b2*s22*g2_20(n) + b3*s32*g3_20(n) + M_p2(20)*M2)-dRg2*g2_20(n);
%     g3_20_dot = phi(c1*s13*g1_20(n) + c2*s23*g2_20(n) + c3*s33*g3_20(n) + M_p2(20)*M3)-dRg3*g3_20(n);
%     g1_21_dot = phi(a1*s11*g1_21(n) + a2*s21*g2_21(n) + a3*s31*g3_21(n) + M_p2(21)*M1)-dRg1*g1_21(n);
%     g2_21_dot = phi(b1*s12*g1_21(n) + b2*s22*g2_21(n) + b3*s32*g3_21(n) + M_p2(21)*M2)-dRg2*g2_21(n);
%     g3_21_dot = phi(c1*s13*g1_21(n) + c2*s23*g2_21(n) + c3*s33*g3_21(n) + M_p2(21)*M3)-dRg3*g3_21(n);
%     g1_22_dot = phi(a1*s11*g1_22(n) + a2*s21*g2_22(n) + a3*s31*g3_22(n) + M_p2(22)*M1)-dRg1*g1_22(n);
%     g2_22_dot = phi(b1*s12*g1_22(n) + b2*s22*g2_22(n) + b3*s32*g3_22(n) + M_p2(22)*M2)-dRg2*g2_22(n);
%     g3_22_dot = phi(c1*s13*g1_22(n) + c2*s23*g2_22(n) + c3*s33*g3_22(n) + M_p2(22)*M3)-dRg3*g3_22(n);
%     g1_23_dot = phi(a1*s11*g1_23(n) + a2*s21*g2_23(n) + a3*s31*g3_23(n) + M_p2(23)*M1)-dRg1*g1_23(n);
%     g2_23_dot = phi(b1*s12*g1_23(n) + b2*s22*g2_23(n) + b3*s32*g3_23(n) + M_p2(23)*M2)-dRg2*g2_23(n);
%     g3_23_dot = phi(c1*s13*g1_23(n) + c2*s23*g2_23(n) + c3*s33*g3_23(n) + M_p2(23)*M3)-dRg3*g3_23(n);
%     g1_24_dot = phi(a1*s11*g1_24(n) + a2*s21*g2_24(n) + a3*s31*g3_24(n) + M_p2(24)*M1)-dRg1*g1_24(n);
%     g2_24_dot = phi(b1*s12*g1_24(n) + b2*s22*g2_24(n) + b3*s32*g3_24(n) + M_p2(24)*M2)-dRg2*g2_24(n);
%     g3_24_dot = phi(c1*s13*g1_24(n) + c2*s23*g2_24(n) + c3*s33*g3_24(n) + M_p2(24)*M3)-dRg3*g3_24(n);
%     g1_25_dot = phi(a1*s11*g1_25(n) + a2*s21*g2_25(n) + a3*s31*g3_25(n) + M_p2(25)*M1)-dRg1*g1_25(n);
%     g2_25_dot = phi(b1*s12*g1_25(n) + b2*s22*g2_25(n) + b3*s32*g3_25(n) + M_p2(25)*M2)-dRg2*g2_25(n);
%     g3_25_dot = phi(c1*s13*g1_25(n) + c2*s23*g2_25(n) + c3*s33*g3_25(n) + M_p2(25)*M3)-dRg3*g3_25(n);
%     g1_26_dot = phi(a1*s11*g1_26(n) + a2*s21*g2_26(n) + a3*s31*g3_26(n) + M_p2(26)*M1)-dRg1*g1_26(n);
%     g2_26_dot = phi(b1*s12*g1_26(n) + b2*s22*g2_26(n) + b3*s32*g3_26(n) + M_p2(26)*M2)-dRg2*g2_26(n);
%     g3_26_dot = phi(c1*s13*g1_26(n) + c2*s23*g2_26(n) + c3*s33*g3_26(n) + M_p2(26)*M3)-dRg3*g3_26(n);
%     g1_27_dot = phi(a1*s11*g1_27(n) + a2*s21*g2_27(n) + a3*s31*g3_27(n) + M_p2(27)*M1)-dRg1*g1_27(n);
%     g2_27_dot = phi(b1*s12*g1_27(n) + b2*s22*g2_27(n) + b3*s32*g3_27(n) + M_p2(27)*M2)-dRg2*g2_27(n);
%     g3_27_dot = phi(c1*s13*g1_27(n) + c2*s23*g2_27(n) + c3*s33*g3_27(n) + M_p2(27)*M3)-dRg3*g3_27(n);
%     g1_28_dot = phi(a1*s11*g1_28(n) + a2*s21*g2_28(n) + a3*s31*g3_28(n) + M_p2(28)*M1)-dRg1*g1_28(n);
%     g2_28_dot = phi(b1*s12*g1_28(n) + b2*s22*g2_28(n) + b3*s32*g3_28(n) + M_p2(28)*M2)-dRg2*g2_28(n);
%     g3_28_dot = phi(c1*s13*g1_28(n) + c2*s23*g2_28(n) + c3*s33*g3_28(n) + M_p2(28)*M3)-dRg3*g3_28(n);
%     g1_29_dot = phi(a1*s11*g1_29(n) + a2*s21*g2_29(n) + a3*s31*g3_29(n) + M_p2(29)*M1)-dRg1*g1_29(n);
%     g2_29_dot = phi(b1*s12*g1_29(n) + b2*s22*g2_29(n) + b3*s32*g3_29(n) + M_p2(29)*M2)-dRg2*g2_29(n);
%     g3_29_dot = phi(c1*s13*g1_29(n) + c2*s23*g2_29(n) + c3*s33*g3_29(n) + M_p2(29)*M3)-dRg3*g3_29(n);
%     g1_30_dot = phi(a1*s11*g1_30(n) + a2*s21*g2_30(n) + a3*s31*g3_30(n) + M_p2(30)*M1)-dRg1*g1_30(n);
%     g2_30_dot = phi(b1*s12*g1_30(n) + b2*s22*g2_30(n) + b3*s32*g3_30(n) + M_p2(30)*M2)-dRg2*g2_30(n);
%     g3_30_dot = phi(c1*s13*g1_30(n) + c2*s23*g2_30(n) + c3*s33*g3_30(n) + M_p2(30)*M3)-dRg3*g3_30(n);
%     g1_31_dot = phi(a1*s11*g1_31(n) + a2*s21*g2_31(n) + a3*s31*g3_31(n) + M_p2(31)*M1)-dRg1*g1_31(n);
%     g2_31_dot = phi(b1*s12*g1_31(n) + b2*s22*g2_31(n) + b3*s32*g3_31(n) + M_p2(31)*M2)-dRg2*g2_31(n);
%     g3_31_dot = phi(c1*s13*g1_31(n) + c2*s23*g2_31(n) + c3*s33*g3_31(n) + M_p2(31)*M3)-dRg3*g3_31(n);
%     g1_32_dot = phi(a1*s11*g1_32(n) + a2*s21*g2_32(n) + a3*s31*g3_32(n) + M_p2(32)*M1)-dRg1*g1_32(n);
%     g2_32_dot = phi(b1*s12*g1_32(n) + b2*s22*g2_32(n) + b3*s32*g3_32(n) + M_p2(32)*M2)-dRg2*g2_32(n);
%     g3_32_dot = phi(c1*s13*g1_32(n) + c2*s23*g2_32(n) + c3*s33*g3_32(n) + M_p2(32)*M3)-dRg3*g3_32(n);  
%   
% 
%     g1_1(n+1) = g1_1(n) + stepsize*g1_1_dot;
%     g2_1(n+1) = g2_1(n) + stepsize*g2_1_dot;
%     g3_1(n+1) = g3_1(n) + stepsize*g3_1_dot;
%     g1_2(n+1) = g1_2(n) + stepsize*g1_2_dot;
%     g2_2(n+1) = g2_2(n) + stepsize*g2_2_dot;
%     g3_2(n+1) = g3_2(n) + stepsize*g3_2_dot;
%     g1_3(n+1) = g1_3(n) + stepsize*g1_3_dot;
%     g2_3(n+1) = g2_3(n) + stepsize*g2_3_dot;
%     g3_3(n+1) = g3_3(n) + stepsize*g3_3_dot;
%     g1_4(n+1) = g1_4(n) + stepsize*g1_4_dot;
%     g2_4(n+1) = g2_4(n) + stepsize*g2_4_dot;
%     g3_4(n+1) = g3_4(n) + stepsize*g3_4_dot;
%     g1_5(n+1) = g1_5(n) + stepsize*g1_5_dot;
%     g2_5(n+1) = g2_5(n) + stepsize*g2_5_dot;
%     g3_5(n+1) = g3_5(n) + stepsize*g3_5_dot;
%     g1_6(n+1) = g1_6(n) + stepsize*g1_6_dot;
%     g2_6(n+1) = g2_6(n) + stepsize*g2_6_dot;
%     g3_6(n+1) = g3_6(n) + stepsize*g3_6_dot;
%     g1_7(n+1) = g1_7(n) + stepsize*g1_7_dot;
%     g2_7(n+1) = g2_7(n) + stepsize*g2_7_dot;
%     g3_7(n+1) = g3_7(n) + stepsize*g3_7_dot;
%     g1_8(n+1) = g1_8(n) + stepsize*g1_8_dot;
%     g2_8(n+1) = g2_8(n) + stepsize*g2_8_dot;
%     g3_8(n+1) = g3_8(n) + stepsize*g3_8_dot;
%     g1_9(n+1) = g1_9(n) + stepsize*g1_9_dot;
%     g2_9(n+1) = g2_9(n) + stepsize*g2_9_dot;
%     g3_9(n+1) = g3_9(n) + stepsize*g3_9_dot;
%     g1_10(n+1) = g1_10(n) + stepsize*g1_10_dot;
%     g2_10(n+1) = g2_10(n) + stepsize*g2_10_dot;
%     g3_10(n+1) = g3_10(n) + stepsize*g3_10_dot;
%     g1_11(n+1) = g1_11(n) + stepsize*g1_11_dot;
%     g2_11(n+1) = g2_11(n) + stepsize*g2_11_dot;
%     g3_11(n+1) = g3_11(n) + stepsize*g3_11_dot;
%     g1_12(n+1) = g1_12(n) + stepsize*g1_12_dot;
%     g2_12(n+1) = g2_12(n) + stepsize*g2_12_dot;
%     g3_12(n+1) = g3_12(n) + stepsize*g3_12_dot;
%     g1_13(n+1) = g1_13(n) + stepsize*g1_13_dot;
%     g2_13(n+1) = g2_13(n) + stepsize*g2_13_dot;
%     g3_13(n+1) = g3_13(n) + stepsize*g3_13_dot;
%     g1_14(n+1) = g1_14(n) + stepsize*g1_14_dot;
%     g2_14(n+1) = g2_14(n) + stepsize*g2_14_dot;
%     g3_14(n+1) = g3_14(n) + stepsize*g3_14_dot;
%     g1_15(n+1) = g1_15(n) + stepsize*g1_15_dot;
%     g2_15(n+1) = g2_15(n) + stepsize*g2_15_dot;
%     g3_15(n+1) = g3_15(n) + stepsize*g3_15_dot;
%     g1_16(n+1) = g1_16(n) + stepsize*g1_16_dot;
%     g2_16(n+1) = g2_16(n) + stepsize*g2_16_dot;
%     g3_16(n+1) = g3_16(n) + stepsize*g3_16_dot;
%     g1_17(n+1) = g1_17(n) + stepsize*g1_17_dot;
%     g2_17(n+1) = g2_17(n) + stepsize*g2_17_dot;
%     g3_17(n+1) = g3_17(n) + stepsize*g3_17_dot;
%     g1_18(n+1) = g1_18(n) + stepsize*g1_18_dot;
%     g2_18(n+1) = g2_18(n) + stepsize*g2_18_dot;
%     g3_18(n+1) = g3_18(n) + stepsize*g3_18_dot;
%     g1_19(n+1) = g1_19(n) + stepsize*g1_19_dot;
%     g2_19(n+1) = g2_19(n) + stepsize*g2_19_dot;
%     g3_19(n+1) = g3_19(n) + stepsize*g3_19_dot;
%     g1_20(n+1) = g1_20(n) + stepsize*g1_20_dot;
%     g2_20(n+1) = g2_20(n) + stepsize*g2_20_dot;
%     g3_20(n+1) = g3_20(n) + stepsize*g3_20_dot;
%     g1_21(n+1) = g1_21(n) + stepsize*g1_21_dot;
%     g2_21(n+1) = g2_21(n) + stepsize*g2_21_dot;
%     g3_21(n+1) = g3_21(n) + stepsize*g3_21_dot;
%     g1_22(n+1) = g1_22(n) + stepsize*g1_22_dot;
%     g2_22(n+1) = g2_22(n) + stepsize*g2_22_dot;
%     g3_22(n+1) = g3_22(n) + stepsize*g3_22_dot;
%     g1_23(n+1) = g1_23(n) + stepsize*g1_23_dot;
%     g2_23(n+1) = g2_23(n) + stepsize*g2_23_dot;
%     g3_23(n+1) = g3_23(n) + stepsize*g3_23_dot;
%     g1_24(n+1) = g1_24(n) + stepsize*g1_24_dot;
%     g2_24(n+1) = g2_24(n) + stepsize*g2_24_dot;
%     g3_24(n+1) = g3_24(n) + stepsize*g3_24_dot;
%     g1_25(n+1) = g1_25(n) + stepsize*g1_25_dot;
%     g2_25(n+1) = g2_25(n) + stepsize*g2_25_dot;
%     g3_25(n+1) = g3_25(n) + stepsize*g3_25_dot;
%     g1_26(n+1) = g1_26(n) + stepsize*g1_26_dot;
%     g2_26(n+1) = g2_26(n) + stepsize*g2_26_dot;
%     g3_26(n+1) = g3_26(n) + stepsize*g3_26_dot;
%     g1_27(n+1) = g1_27(n) + stepsize*g1_27_dot;
%     g2_27(n+1) = g2_27(n) + stepsize*g2_27_dot;
%     g3_27(n+1) = g3_27(n) + stepsize*g3_27_dot;
%     g1_28(n+1) = g1_28(n) + stepsize*g1_28_dot;
%     g2_28(n+1) = g2_28(n) + stepsize*g2_28_dot;
%     g3_28(n+1) = g3_28(n) + stepsize*g3_28_dot;
%     g1_29(n+1) = g1_29(n) + stepsize*g1_29_dot;
%     g2_29(n+1) = g2_29(n) + stepsize*g2_29_dot;
%     g3_29(n+1) = g3_29(n) + stepsize*g3_29_dot;
%     g1_30(n+1) = g1_30(n) + stepsize*g1_30_dot;
%     g2_30(n+1) = g2_30(n) + stepsize*g2_30_dot;
%     g3_30(n+1) = g3_30(n) + stepsize*g3_30_dot;
%     g1_31(n+1) = g1_31(n) + stepsize*g1_31_dot;
%     g2_31(n+1) = g2_31(n) + stepsize*g2_31_dot;
%     g3_31(n+1) = g3_31(n) + stepsize*g3_31_dot;
%     g1_32(n+1) = g1_32(n) + stepsize*g1_32_dot;
%     g2_32(n+1) = g2_32(n) + stepsize*g2_32_dot;
%     g3_32(n+1) = g3_32(n) + stepsize*g3_32_dot;

end
output = [G1;G2;G3];
output = output(:);
% output = [g1_1(numiter) g2_1(numiter) g3_1(numiter) g1_2(numiter) g2_2(numiter) ...
%     g3_2(numiter) g1_3(numiter) g2_3(numiter) g3_3(numiter) g1_4(numiter) g2_4(numiter) ...
%     g3_4(numiter) g1_5(numiter) g2_5(numiter) g3_5(numiter) g1_6(numiter) g2_6(numiter) g3_6(numiter)...
%     g1_7(numiter) g2_7(numiter) g3_7(numiter) g1_8(numiter) g2_8(numiter) g3_8(numiter) g1_9(numiter)...
%     g2_9(numiter) g3_9(numiter) g1_10(numiter) g2_10(numiter) g3_10(numiter) g1_11(numiter) g2_11(numiter)...
%     g3_11(numiter) g1_12(numiter) g2_12(numiter) g3_12(numiter) g1_13(numiter) g2_13(numiter) g3_13(numiter)...
%     g1_14(numiter) g2_14(numiter) g3_14(numiter) g1_15(numiter) g2_15(numiter) g3_15(numiter) g1_16(numiter)...
%     g2_16(numiter) g3_16(numiter) g1_17(numiter) g2_17(numiter) g3_17(numiter) g1_18(numiter) g2_18(numiter)...
%     g3_18(numiter) g1_19(numiter) g2_19(numiter) g3_19(numiter) g1_20(numiter) g2_20(numiter) g3_20(numiter) ...
%     g1_21(numiter) g2_21(numiter) g3_21(numiter) g1_22(numiter) g2_22(numiter) g3_22(numiter) g1_23(numiter)...
%     g2_23(numiter) g3_23(numiter) g1_24(numiter) g2_24(numiter) g3_24(numiter) g1_25(numiter) g2_25(numiter)...
%     g3_25(numiter) g1_26(numiter) g2_26(numiter) g3_26(numiter) g1_27(numiter) g2_27(numiter) g3_27(numiter)...
%     g1_28(numiter) g2_28(numiter) g3_28(numiter) g1_29(numiter) g2_29(numiter) g3_29(numiter) g1_30(numiter)...
%     g2_30(numiter) g3_30(numiter) g1_31(numiter) g2_31(numiter) g3_31(numiter) g1_32(numiter) g2_32(numiter) g3_32(numiter)];

end

