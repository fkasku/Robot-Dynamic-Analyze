%25.01.21
%workout for make up exam
%based on Cenk Karaman
%minor differencess

clear all;close all;clc;
syms the_i a_i alp_i d_i T(alp_i,a_i,d_i,the_i) h1 l1 l2 l3 l4 l5 l6 l7
syms t1 t2 t3 t4 t5 t6 t7 d1 d2 d3 d4 d5 d6 d7 d_d1 d_d2 d_d3 d_d4 d_d5 d_d6 dd_d1 dd_d2 dd_d3 dd_d4 dd_d5 dd_d6
syms th1 th2 th3 th4 th5 th6 d_th1 d_th2 d_th3 d_th4 d_th5 d_th6 dd_th1 dd_th2 dd_th3 dd_th4 dd_th5 dd_th6
syms r1 r12 r13 r21 r22 r23 r31 r32 r33 px py pz M
syms a b c d e f a1 a2 a3 a4 a5 a6 
syms m1 m2 m3 m4 fx fy fz nx ny nz g

 
  
T(alp_i,a_i,d_i,the_i) = [cos(the_i)  -sin(the_i)    0   a_i;
    sin(the_i)*cos(alp_i)   cos(the_i)*cos(alp_i)   -sin(alp_i)   -sin(alp_i)*d_i;
    sin(the_i)*sin(alp_i)   cos(the_i)*sin(alp_i)   cos(alp_i)    cos(alp_i)*d_i;
    0                           0                         0           1];   
    
D_H = [pi/2 0 d1 0 ; 
         -pi/2 0 d2 0 ;
           -pi/2  0  l3  th3;
            pi/2  0   d4  0];


dof = length(D_H(:,1));
TT = cell(dof,1);
TTS = TT;
T_inv = cell(dof,1);
TP = eye(4);
TP_st = cell(dof,1);
T_tr = cell (dof,1);
for i = 1:dof
X = sprintf('transformation matrix %d - %d',i-1,i);
disp(X)
TT{i} = T(D_H(i,1),D_H(i,2),D_H(i,3),D_H(i,4));
TT{i}
TTS{i} = simplify(TT{i});
TP = TP * TT{i};
TP_st{i} = TTS{i}(1:3,1:3).';  % .' operates transpoze of system

T_tr{i} = TTS{i}(1:3,1:3).';   
T_inv{i} = [T_tr{i},-T_tr{i}*TT{i}(1:3,4);0 0 0 1];
end
disp('forward kinematics')
simplify(TP)

w0 = [0 0 0].';
dw0 = [0 0 0].';
v0 = [0 0 0].';
dv0 = [0 g 0].';
Z =    [0 0 1].';

Dp = {d_d1 d_d2 d_th3 d_d4};
DDp = {dd_d1 dd_d2 dd_th3  dd_d4 };
Pc=[ 0   0  -l1/2 ;
      0   0   -l2/2 ;
       0   0  -l3/2 ;
         0   0   -l4/2 ];
M=[m1 m2 m3 m4];

Type = {'P','P','R','P'};
for i=1:dof             %i :1 : robot eksen sayýsý kadar
    syms(['Ixx',num2str(i)])
    syms(['Iyy',num2str(i)])
    syms(['Izz',num2str(i)])
end

for i=1:dof
    Ic{i} = [eval(['Ixx',num2str(i)]) 0     0;...
             0  eval(['Iyy',num2str(i)])    0;...
            0  0   eval(['Izz',num2str(i)])];

end


S = @(x)[0 -x(3) x(2); x(3) 0 -x(1); -x(2) x(1) 0];
%dýþadönük matrislerin yerleþimi
for i = 0:dof-1
    RT{i+1} = TT{i+1}(1:3,1:3).';
    P{i+1} = TT{i+1}(1:3,4);
    if (Type{i+1}=='R')
    if i == 0
         w{i+1} = RT{i+1}*w0 + Dp{i+1}*Z;
         dw{i+1} = RT{i+1}*dw0 + DDp{i+1}*Z + cross(RT{i+1}*w0,(Dp{i+1}*Z));   
         dv{i+1}= RT{i+1}*(cross(dw0,P{i+1}) + cross( w0, cross(w0,P{i+1})) + dv0);
         dvc{i+1} = cross(dw{i+1},Pc(i+1,:).') + cross( w{i+1}, (cross(w{i+1},Pc(i+1, :).'))) + dv{i+1};
         F{i+1} =  M(i+1)*dvc{i+1};
         N{i+1} = Ic{i+1}*dw{i+1} + cross(w{i+1},Ic{i+1}*w{i+1});
         
    else
         w{i+1} = RT{i+1}*w{i} + Dp{i+1}*Z;
         dw{i+1} = RT{i+1}*dw{i} + DDp{i+1}*Z + cross(RT{i+1}*w{i},(Dp{i+1}*Z));
         dv{i+1}= RT{i+1}*(cross(dw{i},P{i+1}) + cross( w{i}, cross(w{i},P{i+1})) + dv{i});
         dvc{i+1} = cross(dw{i+1},Pc(i+1,:).') + cross( w{i+1}, (cross(w{i+1},Pc(i+1, :).'))) + dv{i+1}; 
         F{i+1} =  M(i+1)*dvc{i+1};
         N{i+1} = Ic{i+1}*dw{i+1} + cross(w{i+1},Ic{i+1}*w{i+1});
             
    end
    else 
        
         if i == 0
            w{i+1} = RT{i+1}*w0;
            dw{i+1} = RT{i+1}*dw0;
            dv{i+1} = RT{i+1}*(cross(dw0,P{i+1}) + cross(w0,cross(w0,P{i+1}))  + dv0) + 2*cross(w{i+1},Dp{i+1}*Z)+DDp{i+1}*Z;
            dvc{i+1} = cross(dw{i+1},Pc(i+1,:).')+ cross(w{i+1},(cross(w{i+1},Pc(i+1,:).'))) + dv{i+1};
            F{i+1} = M(i+1)*dvc{i+1};
            N{i+1} = Ic{i+1}*dw{i+1} + cross(w{i+1},(Ic{i+1}*w{i+1}));
        else
            w{i+1} = RT{i+1}*w{i};
            dw{i+1} = RT{i+1}*dw{i};
            dv{i+1} = RT{i+1}*(cross(dw{i},P{i+1}) + cross(w{i},cross(w{i},P{i+1})) + dv{i})+2*cross(w{i+1},Dp{i+1}*Z)+DDp{i+1}*Z;
            dvc{i+1} = cross(dw{i+1},Pc(i+1,:).')+ cross(w{i+1},(cross(w{i+1},Pc(i+1,:).'))) + dv{i+1};
             F{i+1} = M(i+1)*dvc{i+1};  
            N{i+1} = Ic{i+1}*dw{i+1} + cross(w{i+1},(Ic{i+1}*w{i+1}));
        end
    end
end



syms fx fy fz nx ny nz
MyRot= eye(3);
%Myf = [fx;fy;fz];
% Myn = [nx;ny;nz];
Myf=[0;0;0]; %manipulator can move without starting forces   = ff{i+1}
Myn = [0;0;0]; %manipulator can move without starting torque = nn{i+1}

Myp = [0;0;0];
for i = dof:-1:1
    R{i}=TT{i}(1:3,1:3);
    if i == dof
        ff{i}=MyRot*Myf + F{i};
        nn{i} = N{i} + MyRot*Myn + cross(Pc(i,:).',F{i})+cross(Myp,MyRot*Myf);
    else
        ff{i} =R{i+1}*ff{i+1} + F{i};
        nn{i} = N{i} + R{i+1}*nn{i+1} + cross(Pc(i,:).',F{i}) + cross(P{i+1},(R{i+1}*ff{i+1}));
        
    end 
end

%rev joint için hesaplama
for i=1:dof
     if  (Type{i} == 'R')
    Tau{i}=nn{i}.'*Z;    
     end
  if  (Type{i} == 'P')
    Tau{i}=ff{i}.'*Z;    
  end
  
end
    


