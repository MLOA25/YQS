%%--- Eigensolver algorithm 
%%% Using H=h omega(cos(fi)X+sin(fi)Z)/2 -Phase Damping
%%%Saves Data and plot Plotting of both, F(k) and W(k) varying Gamma---%%%
%%%---Mar 2025---%%%%
%%%---M.L.Olivera-Atencio---%%%%

clear 

%---Parameters/Constants---%
Tau=1;%2*pi;          % Real parameter.
r=0.9;          % Reward, heuristic number.
nu=2;           % nu=r*p. 
p=nu/r;         % Punishment, heuristic number.

omega=1;
fi=2*pi/3; 

NDT=5;           % Number of variation of Decohrence Time of the lambda.
N=600;          % Number of iterations. 
NR=1000;           % Number of repetitions.

%---Operators---%
Id=eye(2);          % Identity operator
X=[0 1;1 0];        % Sigma x operator
Y=1i*[0 -1;1 0];    % Sigma y operator
Z=[1 0;0 -1];       % Sigma z operator

%---Basis---%
c=[1;0];            % |0>
u=[0;1];            % |1>
cc=c*c';
cu=c*u';
uc=u*c';
uu=u*u';

%---Hamiltonian---%
H=omega*(sin(fi)*X+cos(fi)*Z)/2 ;  

%---Eigenvectors---%
[V,D]=eig(H);       % Eigenvalues
Au_mm=V(:,1);        % Eigenvector 1 
Au_pp=V(:,2);        % Eigenvector 2 

Au_e = cos(fi/2)*c+sin(fi/2)*u; %excited state
Au_g = -sin(fi/2)*c+cos(fi/2)*u; %excited state

%---Evol Operator---%
U=cos(omega*Tau/2)*Id-1i*sin(omega*Tau/2)*(2*H/omega);

%%Initial state
%IS = c ;
DT_base = 1; % Decoherence Time
DT = DT_base;

%---Auxiliar matrices and vectors---%
k_it=1:N;           

F_DT=zeros(N,NDT);   
w_DT=zeros(N,NDT);  
F_g=zeros(N,NDT); 
F_e=zeros(N,NDT);

%DecTime_NDT=zeros(1,NDT); % Values for legend.
DecTime = zeros(1, NDT);


lim_F = zeros(NDT);
lim_W = zeros(NDT);

%---Initialization for loop DT---%
%---Files for data---%
FileName1=sprintf('Fmean_Tau%1.0f_r%1.0f_nu%1.0f_NR%2.0f.txt',Tau,r*10,nu,NR);
fID1=fopen(FileName1,'wt');

FileName2=sprintf('wmean_Tau%1.0f_r%1.0f_nu%1.0f_NR%2.0f.txt',Tau,r*10,nu,NR);
fID2=fopen(FileName2,'wt');

FileName3=sprintf('Fvsk_g_Tau%1.0f_r%1.0f_nu%1.0f_NR%2.0f.txt',Tau,r*10,nu,NR);
fID3=fopen(FileName3,'wt');

FileName4=sprintf('Fvsk_e_Tau%1.0f_r%1.0f_nu%1.0f_NR%2.0f.txt',Tau,r*10,nu,NR);
fID4=fopen(FileName4,'wt');

FileName5=sprintf('DecTime_Tau%1.0f_r%1.0f_nu%1.0f_NR%2.0f.txt',Tau,r*10,nu,NR);
fID5=fopen(FileName5,'wt');

tic
rng(0); 
%%%---LOOP DT-Start---%%%  
for idt=1:NDT    
    %"Krauss operators E for Phase Damping"
    E0 = Au_e*Au_e'+exp(-Tau/DT)*Au_g*Au_g';
    E1 = sqrt(1-exp(-2*Tau/DT))*Au_g*Au_g';
    
    %---Inicialization for LOOP Repetition---%    
    F_gmean=zeros(1,N); 
    F_emean=zeros(1,N); % extra check 
    F_mean=zeros(1,N);   
    w_Nmean=0;          
    w_mean=zeros(1,N);       
    
    %%%---LOOP Repetitions-Start---%%%
    for ir=1:NR
        rng(ir);
        
        %---Initialization for LOOP Iteration---%   
        IS = c; 
        MS = IS ;
        ro_i = IS*IS'; % Initial state in density operator form, for example |0><0|
        nr=0;           % Number of iterations with reward.
        np=0;           % Number of iterqations with punishment.
        w=1;            % Inicial R/P rate 
        
        F_kg=zeros(1,N); % Fidelity the 1st eigenvector of H, ground state.
        F_ke=zeros(1,N); % Fidelity the 2d eigenvector of H, excited state.
        F_k=zeros(1,N);  % Fidelity for each repetition i.
        
        w_k=zeros(1,N); % R/T rate for each repetition i.
        
        %%%---LOOP Iterations-Start---%%%        
        for k=1:N
            ro_f = U*E0*ro_i*E0'*U'+U*E1*ro_i*E1'*U'; % State after the environment interaction (is above)
            PIS = abs(MS'*ro_f*MS); % Probability of obtaining initial state
            chi=rand;          % Uniformily distribuited random number between 0 and 1. 
            
            if chi<=PIS
                nr=nr+1; %Reward/Punishment control in numbers of iterations
                m=0;
            else
                np=np+1; %Reward/Punishment control in numbers of iterations
                m=1;
            end
                    
            %---Rotations---%            
            phi_y=2*w*pi*rand-w*pi; % Random angle in[-wpi,wpi]. 
            phi_z=2*w*pi*rand-w*pi ;
            phi_x=2*w*pi*rand-w*pi ;
            
            S_xj=1/2*X;             % (eq.3).
            S_yj=1/2*Y;             % (eq.3).
            S_zj=1/2*Z;             % (eq.3).
            
            Uj=(cos(phi_y/2)*Id-1i*sin(phi_y/2)*Y)*(cos(phi_z/2)*Id-1i*sin(phi_z/2)*Z)*(cos(phi_x/2)*Id-1i*sin(phi_x/2)*X); % Unitary transformation, pseudo random operator (eq.2).
            
            Ge=(1-m)*Id+m*Uj;    % Algorithm for aplication of unitary trasformation if it were needed, page 11.
            MS = Ge*MS;      % %Modified state, preparation of the new state for the netx iteration
            ro_i = MS*MS';      % Inicial State for next iteration step.
            
            %---Fidelities---%            
            F_kg(1,k) = abs(Au_g'*MS); % Fidelity with 1st eigenvalue of OE (on iteration k), negative.
            F_ke(1,k) = abs(Au_e'*MS); % Fidelity with 2st eigenvalue of OE (on iteration k), positive.
            
            %---Reward/Punishment---%                
            w_k(1,k)=w;         % R/P rate for each repetition (on iteration k).
            w=((1-m)*r+m*p)*w ; % R/P rate for next iteration. 
            
            if w>1
                w=1;
            end

        end %%%---LOOP Iterations-End---%%%
        F_gmean=F_gmean+F_kg/NR; % //Extra check: mean fidelity of <lE|De|j> for first eigenvalue of OE.
        F_emean=F_emean+F_ke/NR; % //Extra check: mean fidelity of <lE|De|j> for second eigenvalue of OE.
        F_mean=F_mean+max([F_kg;F_ke])/NR; % Mean fidelity(on iteration k) over repetitions: Take the maximun value of fidelities calculated for differents states and then calculate.
       
        w_mean=w_mean+w_k/NR;    % Mean R/P rate (on iteration k) over repetitions.
        w_N=r^(nr-np)*nu^np;     % R/T control (eq.32).
        w_Nmean=w_Nmean+w_N/NR;  % Mean WN R/P rate over repetition.
        
     end %%%--- LOOP Repetition-End---%%%        
        F_g(:,idt)=F_gmean; % Mean fidelity for first eigenvalue of OE for different Gammas. 
        F_e(:,idt)=F_emean; % Mean fidelity for second eigenvalue of OE for different Gammas.
        F_DT(:,idt)=F_mean;   % Mean fidelity global for different Gammas.
        w_DT(:,idt)=w_mean;   % R/T rate for different Gammas.

        lim_F(idt) = mean(F_mean(N-10:N)); 
        lim_W(idt) = mean(w_mean(N-10:N));
        
        fprintf(fID1,' %14.12d',F_DT(:,idt)); % Extra check: Save the fidelity. Each file is a different DT
        fprintf( fID1,'%6s\n','');
        
        fprintf(fID2,' %14.12d',w_DT(:,idt)); % Extra check: Save the fidelity. Each file is a different lDT
        fprintf( fID2,'%6s\n','');
      
        fprintf(fID3,' %14.12d',F_g(:,idt)); % Extra check: Save the fidelity. Each file is a different DT
        fprintf(fID3,'%6s\n','');  
    
        fprintf(fID4,' %14.12d',F_e(:,idt)); % Extra check: Save the fidelity. Each file is a different DT
        fprintf(fID4,'%6s\n','');
        
        DecTime_NDT(1,idt)=DT;

        fprintf(fID5,' %14.12d',DecTime_NDT(:,idt)); % Each file is a different lambda
        fprintf(fID5,'%6s\n','');
        
        DecTime(1,idt)=DT;
        DT=DT*10;
        
end %%%---LOOP DT-End---%%%
toc       
fclose('all');

%---Graphics---%
%---Pannels---%
figure
subplot(2,1,1) % Mean fidelity for different Gammas.
%plot(k_it, F_eig1); hold on
plot(k_it, F_DT,'o-'); hold on
title(['$\tau=$',num2str(Tau),',$r=$',num2str(r), ',$\nu=$',num2str(nu),'$NR=$',num2str(NR)], 'interpreter','latex','Fontsize',14)
lgd=legend( sprintfc('%g', DecTime_NDT) );
title(lgd,'\lambda')
xlabel('$k$', 'interpreter', 'latex'); 
ylabel('$F$', 'interpreter', 'latex');
%ylim([0.7 1]);

subplot(2,1,2)
plot(k_it,w_DT); 
%plot(k_it, F_eig1); hold on
%plot(k_it, F_eig2); hold on
lgd=legend( sprintfc('%g', DecTime_NDT) );
title(lgd,'\lambda')
xlabel('$k$', 'interpreter', 'latex'); 
ylabel('$W$', 'interpreter', 'latex');