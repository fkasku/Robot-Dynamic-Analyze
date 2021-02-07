%30.01.21
%based on Cenk Karaman
%major differencess

clear all;close all;clc;
syms the_i a_i alp_i d_i T(alp_i,a_i,d_i,the_i) h1 l1 l2 l3 l4 l5 l6 l7
syms t1 t2 t3 t4 t5 t6 t7 d1 d2 d3 d4 d5 d6 d7 d_d1 d_d2 d_d3 d_d4 d_d5 d_d6 dd_d1 dd_d2 dd_d3 dd_d4 dd_d5 dd_d6
syms th1 th2 th3 th4 th5 th6 d_th1 d_th2 d_th3 d_th4 d_th5 d_th6 dd_th1 dd_th2 dd_th3 dd_th4 dd_th5 dd_th6
syms r1 r12 r13 r21 r22 r23 r31 r32 r33 px py pz M
syms a b c d e f a1 a2 a3 a4 a5 a6  null
syms m1 m2 m3 m4 fx fy fz nx ny nz g

 
  
T(alp_i,a_i,d_i,the_i) = [cos(the_i)  -sin(the_i)    0   a_i;
    sin(the_i)*cos(alp_i)   cos(the_i)*cos(alp_i)   -sin(alp_i)   -sin(alp_i)*d_i;
    sin(the_i)*sin(alp_i)   cos(the_i)*sin(alp_i)   cos(alp_i)    cos(alp_i)*d_i;
    0                           0                         0           1];   
    

%DH Table
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
tr{i}=TP;
%transformasyon(:,:,i)=TP;  % T01 T02 T03 .... T0N

TP_st{i} = TTS{i}(1:3,1:3).';  % .' operates transpoze of system

T_tr{i} = TTS{i}(1:3,1:3).';    
T_inv{i} = [T_tr{i},-T_tr{i}*TT{i}(1:3,4);0 0 0 1];
end
disp('forward kinematics')
simplify(TP)





%inertia 
for i=1:dof           %i :1 : robot eksen sayýsý kadar
    syms(['Ixx',num2str(i)])
    syms(['Iyy',num2str(i)])
    syms(['Izz',num2str(i)])
end

for i=1:dof
    Ic{i} = [eval(['Ixx',num2str(i)]) 0     0;...
             0  eval(['Iyy',num2str(i)])    0;...
            0  0   eval(['Izz',num2str(i)])];

end


%inertia to base
for i=1:dof
    I{i}=tr{i}(1:3,1:3)*Ic{i}*tr{i}(1:3,1:3).';
end

%dh{i} mass center of joints by their center

dh{1}=[0 0 -l1/2 1].';
dh{2}=[0 0 -l2/2 1].';
dh{3}=[0 0 -l3/2 1].';
dh{4}=[0 0 -l4/2 1].';% if no values put 0 last digit of matrix


%position to base
for i=1:dof
    h{i}=tr{i}*dh{i};
end


% joint variables
m={m1, m2 , m3, m4};  % putt 'null' if no value belongs to last digit  {m1,m2, null}
q={ d1, d2, th3 , d4};  % putt 'null' if no value belongs to last digit  {m1,m2, null}


%derivatives of joint variables
d_q = [d_d1 ; d_d2; d_th3; d_d4];
dd_q = [dd_d1; dd_d2; dd_th3 ; dd_d4];  % put '0' if d_qn=0

%Matrix A of J
for i=1:dof
    
    for j=1:dof
        
  Ax{i,j}= [diff(h{i}(1:3),q{j})];
  
    end
end

%A1= [A{1,1} A{1,2}  A{1,3}]
%A2= [ A{2,1} A{2,2} A{2,3}]

A=cell(dof,1);  % Bu sayede A1 A2 A3 gibi jakobiyen matris elemanlarý hesaplatýlýr
                  % bu iþlemde 3 adet eklem deðiþkeni mevcut idi
                   %dolayýsý ile 3 'e kadar matrise kaydedildi.
for i =1:dof         
   A{i}  = [Ax{i,1} Ax{i,2}  Ax{i,3} Ax{i,4} ];
end

%B matrisi yazýmý
z=[0 0 1].';
Bx=cell(dof,1); 
B=cell(1,dof);
Bx{1}= [0 0 0].';
Bx{2}= [0 0 0].';
Bx{3}= [0 0 0].';
Bx{4}= [0 0 0].';
%B={ 0 0 0; 0 0 0; 0 0 0}

Type = {'P','P','R','P'};   % Type of joints
for i =1 :dof
    if (Type{i} == 'R')
Bx{i}= tr{i}(1:3,1:3)*z;
    else
        Bx{i}=0*tr{i}(1:3,1:3)*z;
    end
    B{i}=[Bx{1}  Bx{2}  Bx{3} Bx{4}];
    
end

%Jakobiyen Matrisleri
D=eye(4)*0;
for i = 1:dof
    J{i}= [A{i};B{i}];

   
Dx{i}= m{i}*A{i}.'*A{i} + B{i}.'*I{i}*B{i};

D=D+Dx{i};

end
 
%finding C matrix
for i = 1:dof
    for k =1:dof
        for j = 1:dof
            temp(k,j) = diff(D(i,j),q{k})-1/2*(diff(D(k,j),q{i}));
        end
    end
    C{i} =temp
end


gravit=0;
grav=cell(dof,1);
dir=1 ; % gravity direction, if negative make 0 , if positive make 1

grav_type={'Y'}; % put the type of gravity direction X Y Z
if (grav_type{1} == 'X')
 for i=1:dof
        for j=1:dof
          gravit= gravit+m{j}*A{j}(1,i);
        end
        grav{i}=(g*(-1)^dir)*gravit;
        gravit=0;
        
 end
end
if (grav_type{1} == 'Y')
    for i=1:dof
        for j=1:dof
          gravit= gravit+m{j}*A{j}(2,i);
        end
        grav{i}=(g*(-1)^dir)*gravit;
        gravit=0;
        
    end
end
if (grav_type{1} == 'Z')
 for i=1:dof
        for j=1:dof
          gravit= gravit+m{j}*A{j}(3,i);
        end
        grav{i}=(g*(-1)^dir)*gravit;
        gravit=0;
        
 end
end
    cd_q=cell(dof,dof);
    for i=1: dof
        for j=1:dof
             cd_q{i,j}= d_q(i)*d_q(j);
        end 
    end
    
    
      cq=cell(dof,1);
      sayac=0;
      for i =1:dof
          for j=1:dof
              for k=1:dof
             temp(j,k) = C{i}(j,k)*cd_q(j,k) ;
             sayac=sayac + temp(j,k);
              end
          end
          cq{i}=sayac;
          sayac=0;
      end
 
 Tau=D*dd_q + cq +grav;
 simplify(Tau)
      
          