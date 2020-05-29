

clear

%%% clearance and Frictional Effect   
%%%%%%%%%%%%%%%%%% 
%%%%% no cl no FE:  odeparams.tar_k=0.4861; odeparams.tau_k=1;  odeparams.alpha_hI=300;
%%%%% modest:       odeparams.tar_k=0.37499; odeparams.tau_k=1; odeparams.alpha_hI=50; 
%%%%% balance:      odeparams.tar_k=0.37499; odeparams.tau_k=1; odeparams.alpha_hI=15; 
%%%%% extreme:      odeparams.tar_k=0.37499; odeparams.tau_k=1; odeparams.alpha_hI=5; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%% begin of silulation of no clearnacs and no FE 


outfile='no_cl_no_FE.xlsx';
odeparams.tau_FE=2880; % no frctional effect odeparams.tau_FE=280000000000 and odeparams.tar_k=0.4861;
odeparams.hyperIbar=0.5; 
odeparams.hyperIk=4;
odeparams.alpha_hI=300; %%%%% modest:alpha_hI=50, balance:alpha_hI=15, extreme:alpha_hI=5 
odeparams.hyperI_sh=5; 
odeparams.hyperI_b=0.5;

odeparams.tar_k=0.4861;% modest, balance, and exteme: odeparams.tar_k=0.37499;
%odeparams.tar_k=0.37499;
odeparams.tau_k=1;

lw3=2;
lw4=1;
%total_t=7200;
total_t=2628000; %paper2
OGTT_period=7200;
nPeriods=total_t/OGTT_period;

%%%% NOTE THAT odeparams.HGP_no_si=1 and  HGP_no_si=1 are used only for HGP test with or
%%%% w/o si dependency, otherwise ALWAYS 0
odeparams.HGP_no_si=0;
HGP_no_si=0; % for recomputing HGP in full_run.m


odeparams.BW=75;
odeparams.mealbar=11.055; 
odeparams.meal=1;
odeparams.OGTT=0;
odeparams.IVGTT=0;
odeparams.r20=0.006;

odeparams.tar_si=0.1; %clearance and frictional effect
odeparams.tau_si=360000; 
odeparams.tar_hepasi=0.85;%clearance and frictional effect
odeparams.tau_hepasi=360000;


odeparams.Gs=100;%clearance and frictional effect
odeparams.ISRI_bar=0.5259; % beta-function defect for all pathways, default:1.4


init=[78.59, 5.63, 1533.91, -0.07663, 1, 0.8, 1, 1, 0.4861, 60.24, 443.39]; % IC   


options= odeset('RelTol',1e-5);
dt=10;
tspan=[0:dt:OGTT_period];
tspanOGTT=[0:1:120];
    
t0=0;
T=[]; %pre-allocate this array instead
Y=[];
OGTT(1+nPeriods).t=[];
OGTT(1+nPeriods).y=[];

%%%% DO an OGTT at baseline

%OGTT at the end of longitudinal simulations
    odeparams.meal=0;
    odeparams.OGTT=1;
    [t1,y1]=ode15s(@pathway_FIG9,tspanOGTT,init,options,odeparams);
    OGTT(1).t=t1;
    OGTT(1).y=y1;

%%%%%%%%%%%%%
%%%% Begin of longitudinal simulations and OGTTs

for i=1:nPeriods

   %inter-OGTTs
    odeparams.meal=1;
    odeparams.OGTT=0;
    odeparams.curt=(i-1)*OGTT_period; % show current time for forced improvement of si, tar_si2 and tau_si2
    
    [t,y]=ode15s(@pathway_FIG9,tspan,init,options,odeparams);


    T=[T;t+(i-1)*OGTT_period];
   
    Y=[Y;y];
    init=y(end,:);
    
    i %remove this to speed up
    
    %OGTT at the end of longitudinal simulations
    odeparams.meal=0;
    odeparams.OGTT=1;
    [t,y]=ode15s(@pathway_FIG9,tspanOGTT,init,options,odeparams);
    OGTT(i+1).t=t;
    OGTT(i+1).y=y;
       
end

%%%%% End of Simulations%%%

T0(1)=0;
for i=2:nPeriods + 1  
T0(i)=T(i-1+(i-1)*OGTT_period/dt); %%%% Don't be confused with # of array and real time

end 

for i=1:nPeriods + 1
    
    
    G0(i)=OGTT(i).y(tspanOGTT==0,1);
    I0(i)=OGTT(i).y(tspanOGTT==0,2);
    
    G30(i)=OGTT(i).y(tspanOGTT==30,1);
    I30(i)=OGTT(i).y(tspanOGTT==30,2);
    
    G60(i)=OGTT(i).y(tspanOGTT==60,1);
    I60(i)=OGTT(i).y(tspanOGTT==60,2);
    
    G90(i)=OGTT(i).y(tspanOGTT==90,1);
    I90(i)=OGTT(i).y(tspanOGTT==90,2);
    
    G120(i)=OGTT(i).y(tspanOGTT==120,1);
    I120(i)=OGTT(i).y(tspanOGTT==120,2);
    
    b(i)=OGTT(i).y(tspanOGTT ==0,3);
    gamma(i)=OGTT(i).y(tspanOGTT ==0,4);
    sigma(i)=OGTT(i).y(tspanOGTT ==0,5);
    si(i)=OGTT(i).y(tspanOGTT ==0,6);
    hepasi(i)=OGTT(i).y(tspanOGTT ==0,7);   
    FE(i)=OGTT(i).y(tspanOGTT ==0,8);  
    k(i)=OGTT(i).y(tspanOGTT ==0,9);


    m_G(i)=mean([G0(i),G30(i),G60(i),G90(i),G120(i)]);
    m_I(i)=mean([I0(i),I30(i),I60(i),I90(i),I120(i)]);
    
    IGI(i)=(I30(i)-I0(i))/(G30(i) - G0(i));
    matsuda(i)=10000/(sqrt(G0(i)*I0(i)*m_G(i)*m_I(i)));
    
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%% recompute HGP

     hepa_bar=15.443; hepa_k=0.27; hepa_b=-3.54277; con_si=0.8;
     hepa_max= hepa_bar./(hepa_k +OGTT(i).y(tspanOGTT==0,6).*(1-HGP_no_si) + con_si*HGP_no_si) + hepa_b;

     alpha_max=6; alpha_k=0.4; alpha_b=-0.5;
     alpha_HGP= alpha_max./(alpha_k + OGTT(i).y(tspanOGTT==0,6).*(1-HGP_no_si) + con_si*HGP_no_si) + alpha_b;

     HGP_b=0.104166;  
     HGP(i) = hepa_max./(alpha_HGP + OGTT(i).y(tspanOGTT==0,2)*OGTT(i).y(tspanOGTT==0,7)) + HGP_b;


%%%%%end of recompute HGP

end

%%% write data

long_OGTT=[T0(:),G0(:),G60(:),G120(:),I0(:),I60(:),I120(:),HGP(:),b(:),gamma(:),sigma(:),si(:),hepasi(:),FE(:),k(:)];
tmp_long_OGTT=array2table(long_OGTT,'VariableNAMES',{'t','G0','G60','G120','I0','I60','I120','HGP','b','gamma','sigma','si','hepasi','FE','k'});
writetable(tmp_long_OGTT,outfile) 


%%%% end of simulation of no clearance and no FE
%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%% begin of silulation of modest clearance and FE 

clear

outfile='modest_cl_FE.xlsx';
odeparams.tau_FE=2880; 
odeparams.hyperIbar=0.5; 
odeparams.hyperIk=4;
odeparams.alpha_hI=50; %%%%% modest:alpha_hI=50,balance:alpha_hI=15, extreme:alpha_hI=5 
odeparams.hyperI_sh=5; 
odeparams.hyperI_b=0.5;


 
%odeparams.tar_k=0.4861;% modest, balance, and exteme: odeparams.tar_k=0.37499;
odeparams.tar_k=0.37499;
odeparams.tau_k=1;

th_G60=155;
th_G120=140;
lw3=2;
lw4=1;
%total_t=7200;
total_t=2628000; %paper2
OGTT_period=7200;
nPeriods=total_t/OGTT_period;

%%%% NOTE THAT odeparams.HGP_no_si=1 and  HGP_no_si=1 are used only for HGP test with or
%%%% w/o si dependency, otherwise ALWAYS 0
odeparams.HGP_no_si=0;
HGP_no_si=0; % for recomputing HGP in full_run.m


odeparams.BW=75;
odeparams.mealbar=11.055; 
odeparams.meal=1;
odeparams.OGTT=0;
odeparams.IVGTT=0;
odeparams.r20=0.006;

odeparams.tar_si=0.1; %clearance and frictional effect
odeparams.tau_si=360000; 
odeparams.tar_hepasi=0.85;%clearance and frictional effect
odeparams.tau_hepasi=360000;


odeparams.Gs=100;%clearance and frictional effect
odeparams.ISRI_bar=0.5259; % beta-function defect for all pathways, default:1.4


init=[78.59, 5.63, 1533.91, -0.07663, 1, 0.8, 1, 1, 0.4861, 60.24, 443.39]; % IC   


options= odeset('RelTol',1e-5);
dt=10;
tspan=[0:dt:OGTT_period];
tspanOGTT=[0:1:120];
    
t0=0;
T=[]; %pre-allocate this array instead
Y=[];
OGTT(1+nPeriods).t=[];
OGTT(1+nPeriods).y=[];

%%%% DO an OGTT at baseline

%OGTT at the end of longitudinal simulations
    odeparams.meal=0;
    odeparams.OGTT=1;
    [t1,y1]=ode15s(@pathway_FIG9,tspanOGTT,init,options,odeparams);
    OGTT(1).t=t1;
    OGTT(1).y=y1;

%%%%%%%%%%%%%
%%%% Begin of longitudinal simulations and OGTTs

for i=1:nPeriods

   %inter-OGTTs
    odeparams.meal=1;
    odeparams.OGTT=0;
    odeparams.curt=(i-1)*OGTT_period; % show current time for forced improvement of si, tar_si2 and tau_si2
    
    [t,y]=ode15s(@pathway_FIG9,tspan,init,options,odeparams);


    T=[T;t+(i-1)*OGTT_period];
   
    Y=[Y;y];
    init=y(end,:);
    
    i %remove this to speed up
    
    %OGTT at the end of longitudinal simulations
    odeparams.meal=0;
    odeparams.OGTT=1;
    [t,y]=ode15s(@pathway_FIG9,tspanOGTT,init,options,odeparams);
    OGTT(i+1).t=t;
    OGTT(i+1).y=y;
       
end

%%%%% End of Simulations%%%

T0(1)=0;
for i=2:nPeriods + 1  

T0(i)=T(i-1+(i-1)*OGTT_period/dt); %%%% Don't be confused with # of array and real time

end 

for i=1:nPeriods + 1
    
    
    G0(i)=OGTT(i).y(tspanOGTT==0,1);
    I0(i)=OGTT(i).y(tspanOGTT==0,2);
    
    G30(i)=OGTT(i).y(tspanOGTT==30,1);
    I30(i)=OGTT(i).y(tspanOGTT==30,2);
    
    G60(i)=OGTT(i).y(tspanOGTT==60,1);
    I60(i)=OGTT(i).y(tspanOGTT==60,2);
    
    G90(i)=OGTT(i).y(tspanOGTT==90,1);
    I90(i)=OGTT(i).y(tspanOGTT==90,2);
    
    G120(i)=OGTT(i).y(tspanOGTT==120,1);
    I120(i)=OGTT(i).y(tspanOGTT==120,2);
    
    b(i)=OGTT(i).y(tspanOGTT ==0,3);
    gamma(i)=OGTT(i).y(tspanOGTT ==0,4);
    sigma(i)=OGTT(i).y(tspanOGTT ==0,5);
    si(i)=OGTT(i).y(tspanOGTT ==0,6);
    hepasi(i)=OGTT(i).y(tspanOGTT ==0,7);   
    FE(i)=OGTT(i).y(tspanOGTT ==0,8);  
    k(i)=OGTT(i).y(tspanOGTT ==0,9);


    m_G(i)=mean([G0(i),G30(i),G60(i),G90(i),G120(i)]);
    m_I(i)=mean([I0(i),I30(i),I60(i),I90(i),I120(i)]);
    
    IGI(i)=(I30(i)-I0(i))/(G30(i) - G0(i));
    matsuda(i)=10000/(sqrt(G0(i)*I0(i)*m_G(i)*m_I(i)));
    
   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%% recompute HGP

     hepa_bar=15.443; hepa_k=0.27; hepa_b=-3.54277; con_si=0.8;
     hepa_max= hepa_bar./(hepa_k +OGTT(i).y(tspanOGTT==0,6).*(1-HGP_no_si) + con_si*HGP_no_si) + hepa_b;

     alpha_max=6; alpha_k=0.4; alpha_b=-0.5;
     alpha_HGP= alpha_max./(alpha_k + OGTT(i).y(tspanOGTT==0,6).*(1-HGP_no_si) + con_si*HGP_no_si) + alpha_b;

     HGP_b=0.104166;  
     HGP(i) = hepa_max./(alpha_HGP + OGTT(i).y(tspanOGTT==0,2)*OGTT(i).y(tspanOGTT==0,7)) + HGP_b;


%%%%%%%%% end of recompute HGP

end

%%%%%% write data

long_OGTT=[T0(:),G0(:),G60(:),G120(:),I0(:),I60(:),I120(:),HGP(:),b(:),gamma(:),sigma(:),si(:),hepasi(:),FE(:),k(:)];
tmp_long_OGTT=array2table(long_OGTT,'VariableNAMES',{'t','G0','G60','G120','I0','I60','I120','HGP','b','gamma','sigma','si','hepasi','FE','k'});
writetable(tmp_long_OGTT,outfile) 



%%%%%%%%%% End of simulationof modest clearnace and FE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%% Begin of simulations of balance clearance and FE
clear

outfile='balance_cl_FE.xlsx';
odeparams.tau_FE=2880; 
odeparams.hyperIbar=0.5; 
odeparams.hyperIk=4;
odeparams.alpha_hI=15; %%%%% modest:alpha_hI=50,balance:alpha_hI=15, extreme:alpha_hI=5 
odeparams.hyperI_sh=5; 
odeparams.hyperI_b=0.5;
 
%odeparams.tar_k=0.4861;% modest, balance, and exteme: odeparams.tar_k=0.37499;
odeparams.tar_k=0.37499;
odeparams.tau_k=1;

th_G60=155;
th_G120=140;
lw3=2;
lw4=1;
%total_t=7200;
total_t=2628000; %paper2
OGTT_period=7200;
nPeriods=total_t/OGTT_period;

%%%% NOTE THAT odeparams.HGP_no_si=1 and  HGP_no_si=1 are used only for HGP test with or
%%%% w/o si dependency, otherwise ALWAYS 0
odeparams.HGP_no_si=0;
HGP_no_si=0; % for recomputing HGP in full_run.m


odeparams.BW=75;
odeparams.mealbar=11.055; 
odeparams.meal=1;
odeparams.OGTT=0;
odeparams.IVGTT=0;
odeparams.r20=0.006;

odeparams.tar_si=0.1; %clearance and frictional effect
odeparams.tau_si=360000; 
odeparams.tar_hepasi=0.85;%clearance and frictional effect
odeparams.tau_hepasi=360000;


odeparams.Gs=100;%clearance and frictional effect
odeparams.ISRI_bar=0.5259; % beta-function defect for all pathways, default:1.4


init=[78.59, 5.63, 1533.91, -0.07663, 1, 0.8, 1, 1, 0.4861, 60.24, 443.39]; % IC   


options= odeset('RelTol',1e-5);
dt=10;
tspan=[0:dt:OGTT_period];
tspanOGTT=[0:1:120];
    
t0=0;
T=[]; %pre-allocate this array instead
Y=[];
OGTT(1+nPeriods).t=[];
OGTT(1+nPeriods).y=[];

%%%% DO an OGTT at baseline

%OGTT at the end of longitudinal simulations
    odeparams.meal=0;
    odeparams.OGTT=1;
    [t1,y1]=ode15s(@pathway_FIG9,tspanOGTT,init,options,odeparams);
    OGTT(1).t=t1;
    OGTT(1).y=y1;

%%%%%%%%%%%%%
%%%% Begin of longitudinal simulations and OGTTs

for i=1:nPeriods

   %inter-OGTTs
    odeparams.meal=1;
    odeparams.OGTT=0;
    odeparams.curt=(i-1)*OGTT_period; % show current time for forced improvement of si, tar_si2 and tau_si2
    
    [t,y]=ode15s(@pathway_FIG9,tspan,init,options,odeparams);


    T=[T;t+(i-1)*OGTT_period];
   
    Y=[Y;y];
    init=y(end,:);
    
    i %remove this to speed up
    
    %OGTT at the end of longitudinal simulations
    odeparams.meal=0;
    odeparams.OGTT=1;
    [t,y]=ode15s(@pathway_FIG9,tspanOGTT,init,options,odeparams);
    OGTT(i+1).t=t;
    OGTT(i+1).y=y;
       
end

%%%%% End of Simulations%%%

T0(1)=0;
for i=2:nPeriods + 1  

T0(i)=T(i-1+(i-1)*OGTT_period/dt); %%%% Don't be confused with # of array and real time

end 

for i=1:nPeriods + 1
    
    
    G0(i)=OGTT(i).y(tspanOGTT==0,1);
    I0(i)=OGTT(i).y(tspanOGTT==0,2);
    
    G30(i)=OGTT(i).y(tspanOGTT==30,1);
    I30(i)=OGTT(i).y(tspanOGTT==30,2);
    
    G60(i)=OGTT(i).y(tspanOGTT==60,1);
    I60(i)=OGTT(i).y(tspanOGTT==60,2);
    
    G90(i)=OGTT(i).y(tspanOGTT==90,1);
    I90(i)=OGTT(i).y(tspanOGTT==90,2);
    
    G120(i)=OGTT(i).y(tspanOGTT==120,1);
    I120(i)=OGTT(i).y(tspanOGTT==120,2);
    
    b(i)=OGTT(i).y(tspanOGTT ==0,3);
    gamma(i)=OGTT(i).y(tspanOGTT ==0,4);
    sigma(i)=OGTT(i).y(tspanOGTT ==0,5);
    si(i)=OGTT(i).y(tspanOGTT ==0,6);
    hepasi(i)=OGTT(i).y(tspanOGTT ==0,7);   
    FE(i)=OGTT(i).y(tspanOGTT ==0,8);  
    k(i)=OGTT(i).y(tspanOGTT ==0,9);


    m_G(i)=mean([G0(i),G30(i),G60(i),G90(i),G120(i)]);
    m_I(i)=mean([I0(i),I30(i),I60(i),I90(i),I120(i)]);
    
    IGI(i)=(I30(i)-I0(i))/(G30(i) - G0(i));
    matsuda(i)=10000/(sqrt(G0(i)*I0(i)*m_G(i)*m_I(i)));
    
   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%% recompute HGP

     hepa_bar=15.443; hepa_k=0.27; hepa_b=-3.54277; con_si=0.8;
     hepa_max= hepa_bar./(hepa_k +OGTT(i).y(tspanOGTT==0,6).*(1-HGP_no_si) + con_si*HGP_no_si) + hepa_b;

     alpha_max=6; alpha_k=0.4; alpha_b=-0.5;
     alpha_HGP= alpha_max./(alpha_k + OGTT(i).y(tspanOGTT==0,6).*(1-HGP_no_si) + con_si*HGP_no_si) + alpha_b;

     HGP_b=0.104166;  
     HGP(i) = hepa_max./(alpha_HGP + OGTT(i).y(tspanOGTT==0,2)*OGTT(i).y(tspanOGTT==0,7)) + HGP_b;


%%%%%%% end of recompute HGP

end

%%% write data

long_OGTT=[T0(:),G0(:),G60(:),G120(:),I0(:),I60(:),I120(:),HGP(:),b(:),gamma(:),sigma(:),si(:),hepasi(:),FE(:),k(:)];
tmp_long_OGTT=array2table(long_OGTT,'VariableNAMES',{'t','G0','G60','G120','I0','I60','I120','HGP','b','gamma','sigma','si','hepasi','FE','k'});
writetable(tmp_long_OGTT,outfile) 


%%%%%%% End of simulations of balance clerance and FE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%% Begin of simulations of extreme and FE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

outfile='extreme_cl_FE.xlsx';
odeparams.tau_FE=2880; 
odeparams.hyperIbar=0.5; 
odeparams.hyperIk=4;
odeparams.alpha_hI=5; %%%%% modest:alpha_hI=50,balance:alpha_hI=15, extreme:alpha_hI=5 
odeparams.hyperI_sh=5; 
odeparams.hyperI_b=0.5;
 
%odeparams.tar_k=0.4861;% modest, balance, and exteme: odeparams.tar_k=0.37499;
odeparams.tar_k=0.37499;
odeparams.tau_k=1;

th_G60=155;
th_G120=140;
lw3=2;
lw4=1;
%total_t=7200;
total_t=2628000; %paper2
OGTT_period=7200;
nPeriods=total_t/OGTT_period;

%%%% NOTE THAT odeparams.HGP_no_si=1 and  HGP_no_si=1 are used only for HGP test with or
%%%% w/o si dependency, otherwise ALWAYS 0
odeparams.HGP_no_si=0;
HGP_no_si=0; % for recomputing HGP in full_run.m


odeparams.BW=75;
odeparams.mealbar=11.055; 
odeparams.meal=1;
odeparams.OGTT=0;
odeparams.IVGTT=0;
odeparams.r20=0.006;

odeparams.tar_si=0.1; %clearance and frictional effect
odeparams.tau_si=360000; 
odeparams.tar_hepasi=0.85;%clearance and frictional effect
odeparams.tau_hepasi=360000;


odeparams.Gs=100;%clearance and frictional effect
odeparams.ISRI_bar=0.5259; % beta-function defect for all pathways, default:1.4


init=[78.59, 5.63, 1533.91, -0.07663, 1, 0.8, 1, 1, 0.4861, 60.24, 443.39]; % IC   


options= odeset('RelTol',1e-5);
dt=10;
tspan=[0:dt:OGTT_period];
tspanOGTT=[0:1:120];
    
t0=0;
T=[]; %pre-allocate this array instead
Y=[];
OGTT(1+nPeriods).t=[];
OGTT(1+nPeriods).y=[];

%%%% DO an OGTT at baseline

%OGTT at the end of longitudinal simulations
    odeparams.meal=0;
    odeparams.OGTT=1;
    [t1,y1]=ode15s(@pathway_FIG9,tspanOGTT,init,options,odeparams);
    OGTT(1).t=t1;
    OGTT(1).y=y1;

%%%%%%%%%%%%%
%%%% Begin of longitudinal simulations and OGTTs

for i=1:nPeriods

   %inter-OGTTs
    odeparams.meal=1;
    odeparams.OGTT=0;
    odeparams.curt=(i-1)*OGTT_period; % show current time for forced improvement of si, tar_si2 and tau_si2
    
    [t,y]=ode15s(@pathway_FIG9,tspan,init,options,odeparams);


    T=[T;t+(i-1)*OGTT_period];
   
    Y=[Y;y];
    init=y(end,:);
    
    i %remove this to speed up
    
    %OGTT at the end of longitudinal simulations
    odeparams.meal=0;
    odeparams.OGTT=1;
    [t,y]=ode15s(@pathway_FIG9,tspanOGTT,init,options,odeparams);
    OGTT(i+1).t=t;
    OGTT(i+1).y=y;
       
end

%%%%% End of Simulations%%%

T0(1)=0;
for i=2:nPeriods + 1  

T0(i)=T(i-1+(i-1)*OGTT_period/dt); %%%% Don't be confused with # of array and real time

end 

for i=1:nPeriods + 1
    
    
    G0(i)=OGTT(i).y(tspanOGTT==0,1);
    I0(i)=OGTT(i).y(tspanOGTT==0,2);
    
    G30(i)=OGTT(i).y(tspanOGTT==30,1);
    I30(i)=OGTT(i).y(tspanOGTT==30,2);
    
    G60(i)=OGTT(i).y(tspanOGTT==60,1);
    I60(i)=OGTT(i).y(tspanOGTT==60,2);
    
    G90(i)=OGTT(i).y(tspanOGTT==90,1);
    I90(i)=OGTT(i).y(tspanOGTT==90,2);
    
    G120(i)=OGTT(i).y(tspanOGTT==120,1);
    I120(i)=OGTT(i).y(tspanOGTT==120,2);
    
    b(i)=OGTT(i).y(tspanOGTT ==0,3);
    gamma(i)=OGTT(i).y(tspanOGTT ==0,4);
    sigma(i)=OGTT(i).y(tspanOGTT ==0,5);
    si(i)=OGTT(i).y(tspanOGTT ==0,6);
    hepasi(i)=OGTT(i).y(tspanOGTT ==0,7);   
    FE(i)=OGTT(i).y(tspanOGTT ==0,8);  
    k(i)=OGTT(i).y(tspanOGTT ==0,9);


    m_G(i)=mean([G0(i),G30(i),G60(i),G90(i),G120(i)]);
    m_I(i)=mean([I0(i),I30(i),I60(i),I90(i),I120(i)]);
    
    IGI(i)=(I30(i)-I0(i))/(G30(i) - G0(i));
    matsuda(i)=10000/(sqrt(G0(i)*I0(i)*m_G(i)*m_I(i)));
    
   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%% recompute HGP

     hepa_bar=15.443; hepa_k=0.27; hepa_b=-3.54277; con_si=0.8;
     hepa_max= hepa_bar./(hepa_k +OGTT(i).y(tspanOGTT==0,6).*(1-HGP_no_si) + con_si*HGP_no_si) + hepa_b;

     alpha_max=6; alpha_k=0.4; alpha_b=-0.5;
     alpha_HGP= alpha_max./(alpha_k + OGTT(i).y(tspanOGTT==0,6).*(1-HGP_no_si) + con_si*HGP_no_si) + alpha_b;

     HGP_b=0.104166;  
     HGP(i) = hepa_max./(alpha_HGP + OGTT(i).y(tspanOGTT==0,2)*OGTT(i).y(tspanOGTT==0,7)) + HGP_b;


%%%%%%%% end of recompute HGP

end

%%% write data

long_OGTT=[T0(:),G0(:),G60(:),G120(:),I0(:),I60(:),I120(:),HGP(:),b(:),gamma(:),sigma(:),si(:),hepasi(:),FE(:),k(:)];
tmp_long_OGTT=array2table(long_OGTT,'VariableNAMES',{'t','G0','G60','G120','I0','I60','I120','HGP','b','gamma','sigma','si','hepasi','FE','k'});
writetable(tmp_long_OGTT,outfile) 


%%%%%%% End of simulations of extreme and FE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% plot data

clear 
DATA=readtable('no_cl_no_FE.xlsx');
%DATA2=readtable('cl_no_FE.xlsx');
DATA3=readtable('modest_cl_FE.xlsx');
DATA4=readtable('balance_cl_FE.xlsx');
DATA5=readtable('extreme_cl_FE.xlsx');


fs=12;
fn='arial';

%%% DATA

DATA_t=DATA{:,1};
DATA_FG=DATA{:,2};
DATA_2hG=DATA{:,4};
DATA_FI=DATA{:,5};
DATA_2hI=DATA{:,7};

DATA_si=DATA{:,12};
DATA_HGP=DATA{:,13};
DATA_sigma=DATA{:,11};
DATA_b=DATA{:,9};
DATA_FE=DATA{:,14};
DATA_k=DATA{:,15};

%%% DATA2

%DATA2_t=DATA2{:,1};
%DATA2_FG=DATA2{:,2};
%DATA2_2hG=DATA2{:,4};
%DATA2_FI=DATA2{:,5};
%DATA2_2hI=DATA2{:,7};

%DATA2_si=DATA2{:,12};
%DATA2_HGP=DATA2{:,13};
%DATA2_sigma=DATA2{:,11};
%DATA2_b=DATA2{:,9};
%DATA2_FE=DATA2{:,14};
%DATA2_k=DATA2{:,15};


%%% DATA3

DATA3_t=DATA3{:,1};
DATA3_FG=DATA3{:,2};
DATA3_2hG=DATA3{:,4};
DATA3_FI=DATA3{:,5};
DATA3_2hI=DATA3{:,7};

DATA3_si=DATA3{:,12};
DATA3_HGP=DATA3{:,13};
DATA3_sigma=DATA3{:,11};
DATA3_b=DATA3{:,9};
DATA3_FE=DATA3{:,14};
DATA3_k=DATA3{:,15};

%%% DATA4

DATA4_t=DATA4{:,1};
DATA4_FG=DATA4{:,2};
DATA4_2hG=DATA4{:,4};
DATA4_FI=DATA4{:,5};
DATA4_2hI=DATA4{:,7};

DATA4_si=DATA4{:,12};
DATA4_HGP=DATA4{:,13};
DATA4_sigma=DATA4{:,11};
DATA4_b=DATA4{:,9};
DATA4_FE=DATA4{:,14};
DATA4_k=DATA4{:,15};

%%% DATA5

DATA5_t=DATA5{:,1};
DATA5_FG=DATA5{:,2};
DATA5_2hG=DATA5{:,4};
DATA5_FI=DATA5{:,5};
DATA5_2hI=DATA5{:,7};

DATA5_si=DATA5{:,12};
DATA5_HGP=DATA5{:,13};
DATA5_sigma=DATA5{:,11};
DATA5_b=DATA5{:,9};
DATA5_FE=DATA5{:,14};
DATA5_k=DATA5{:,15};



% Frictional Effect
XI=linspace(0,300,3001);


hyperIbar=0.5; hyperIk=4; hyperI_sh=5; hyperI_b=0.5;
alpha_hI3=50;
alpha_hI4=15;
alpha_hI5=5;

hyperI_Si=1 + 0.*XI;
hyperI_Si2=1 + 0.*XI;

hyperI_Si3=hyperIbar.*(1- (XI-hyperI_sh).^hyperIk./(alpha_hI3.^hyperIk + (XI-hyperI_sh).^hyperIk) )+ hyperI_b;
hyperI_Si4=hyperIbar.*(1- (XI-hyperI_sh).^hyperIk./(alpha_hI4.^hyperIk + (XI-hyperI_sh).^hyperIk) )+ hyperI_b;
hyperI_Si5=hyperIbar.*(1- (XI-hyperI_sh).^hyperIk./(alpha_hI5.^hyperIk + (XI-hyperI_sh).^hyperIk) )+ hyperI_b;

%%% Note that t should be rescaled by 365*1440 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Begin of plot FIG
fpan=10;
fs=12;
fs2=8;
fn='arial';
lw=1.5;
lw2=1.5;
lw3=0.5;
%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%

tsi=subplot(3,2,1);

plot(DATA_t./(365*1440),DATA_si.*6.9444.*DATA_FE,'k', 'linewidth',lw);
hold('on')
plot(DATA3_t./(365*1440),DATA3_si.*6.9444.*DATA3_FE,'g', 'linewidth',lw);
plot(DATA4_t./(365*1440),DATA4_si.*6.9444.*DATA4_FE,'b', 'linewidth',lw);
plot(DATA5_t./(365*1440),DATA5_si.*6.9444.*DATA5_FE,'r', 'linewidth',lw);
%xlabel('time (year)','fontsize', fs, 'fontname',fn);
ylabel('{S_I} (10^{-4}ml/\muU/min)','fontsize', fs, 'fontname',fn);

%lh=legend('No cl and No FE','Modest cl and FE','Balanced cl and FE','Extreme cl and FE','location','northeast');
%set(lh,'FontSize',7.5); 
%legend('boxoff');


%legend('compensation',5,'decompensation',5);
%legend('boxoff');

%%%%% 10% of the length of y axis


text(0,6.6,'A','fontsize',fpan,'fontweight','bold');
axis ([0 5 0 6])
%%%%%%%%%%%%%%%%%%%%%%%%%%
tHGP=subplot(3,2,2);

plot(DATA_t./(365*1440),DATA_HGP,'k', 'linewidth',lw);
%hold('on')



%xlabel('time (year)','fontsize', fs, 'fontname',fn);

ylabel('{hepa_{S_I}}','fontsize', fs, 'fontname',fn);

%%%%%% 10% of the length of y axis

text(0,1.65,'B','fontsize',fpan,'fontweight','bold');
axis ([0 5 0 1.5])


%%%%% 10%  of the length of y axis



%%%%%%%%%%%%%%%%%%%%%%%%

tG=subplot(3,2,3);


plot(DATA_t./(365*1440),DATA_FG,'k', 'linewidth',lw);
hold('on')
plot(DATA_t./(365*1440),DATA_2hG,'k', 'linewidth',lw);

plot(DATA3_t./(365*1440),DATA3_FG,'g', 'linewidth',lw);
plot(DATA3_t./(365*1440),DATA3_2hG,'g', 'linewidth',lw);

plot(DATA4_t./(365*1440),DATA4_FG,'b', 'linewidth',lw);
plot(DATA4_t./(365*1440),DATA4_2hG,'b', 'linewidth',lw);

plot(DATA5_t./(365*1440),DATA5_FG,'r', 'linewidth',lw);
plot(DATA5_t./(365*1440),DATA5_2hG,'r', 'linewidth',lw);




tIGT=1;
tCGI=1.92;
tT2D=3.00;

th_IGT=140;
th_CGI=100;
th_T2D=200;

%plot([tIGT tIGT],[0 th_IGT],'k','linewidth',lw3);
%plot([tCGI tCGI],[0 th_CGI],'k','linewidth',lw3);
%plot([tT2D tT2D],[0 th_T2D],'k','linewidth',lw3);

%%%%%% longitudinal location at t=0 (NGT1), t=1*319680 (NGT2), t=2*319680 (IGT), t=4*319680 (CGI), t=7*319680

%plot(0,0,'.k','MarkerSize',15);
%plot(0.685,0,'.k','MarkerSize',15);
%plot(1.216, 0, '.k', 'MarkerSize',15);
%plot(2.432, 0, '.k', 'MarkerSize',15);
%plot(4.256, 0, '.k', 'MarkerSize',15);


%hold('on')
%line([1.1 1.1],[0 140]);
%hold('on');
%line([5 0],[5 200]);


%xlabel('time (year)','fontsize', fs, 'fontname',fn);
ylabel('{Glucose} (mg/dl)','fontsize', fs, 'fontname',fn);


%%%%%% 10% of the length of y axis

text(0,440,'C','fontsize',fpan,'fontweight','bold');
text(1,200,'2hPG', 'fontsize',fs2,'fontname',fn);
text(1,120, 'FPG','fontsize',fs2,'fontname',fn);

%text(6.02,280,{'2-h', 'glucose'}, 'fontsize',fs2,'fontname',fn);
%text(6.02,130, {'fasting', 'glucose'},'fontsize',fs2,'fontname',fn);

fs10=8;
%text(0.06,40,'NGT','fontsize',fs10,'fontname',fn);
%text(0.06,40,'NGT1,2','fontsize',fs10,'fontname',fn);
%text(1.3,40,'IGT','fontsize',fs10,'fontname',fn);
%text(2.2,40,'CGI','fontsize',fs10,'fontname',fn);
%text(4.1,40,'T2D','fontsize',fs10,'fontname',fn);

axis ([0 5 0 400])




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tI=subplot(3,2,4);
plot(DATA_t./(365*1440),DATA_FI,'k', 'linewidth',lw);

hold('on')
plot(DATA_t./(365*1440),DATA_2hI,'k', 'linewidth',lw);


plot(DATA3_t./(365*1440),DATA3_FI,'g', 'linewidth',lw);
plot(DATA3_t./(365*1440),DATA3_2hI,'g', 'linewidth',lw);

plot(DATA4_t./(365*1440),DATA4_FI,'b', 'linewidth',lw);
plot(DATA4_t./(365*1440),DATA4_2hI,'b', 'linewidth',lw);

plot(DATA5_t./(365*1440),DATA5_FI,'r', 'linewidth',lw);
plot(DATA5_t./(365*1440),DATA5_2hI,'r', 'linewidth',lw);


%xlabel('time (year)','fontsize', fs, 'fontname',fn);
ylabel('{Insulin} ({\mu}U/ml)','fontsize', fs, 'fontname',fn);


text(0,220,'D','fontsize',fpan,'fontweight','bold');
text(1,170,'2hPI', 'fontsize',fs2,'fontname',fn);
text(1,32, 'FPI','fontsize',fs2,'fontname',fn);


axis ([0 5 0 200])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



tsigma=subplot(3,2,5);

plot(DATA_t./(365*1440),DATA_sigma,'k', 'linewidth',lw);
hold('on')
plot(DATA3_t./(365*1440),DATA3_sigma,'g', 'linewidth',lw);
plot(DATA4_t./(365*1440),DATA4_sigma,'b', 'linewidth',lw);
plot(DATA5_t./(365*1440),DATA5_sigma,'r', 'linewidth',lw);


xlabel('time (year)','fontsize', fs, 'fontname',fn);
ylabel('{\sigma}','fontsize', fs, 'fontname',fn);


%%%%%%% 10% of the length of y axis

text(0,1.1,'E','fontsize',fpan,'fontweight','bold');
axis ([0 5 0 1])


%%%%%%%% 10%  of the length of y axis


tFE=subplot(3,2,6);

%%% I vs FE
plot(XI,hyperI_Si,'k', 'linewidth',lw);
hold('on')
plot(XI,hyperI_Si3,'g', 'linewidth',lw);
plot(XI,hyperI_Si4,'b', 'linewidth',lw);
plot(XI,hyperI_Si5,'r', 'linewidth',lw);

%%%% t vs FE
%plot(DATA_t./(365*1440),DATA_FE,'k', 'linewidth',lw);
%hold('on')
%plot(DATA3_t./(365*1440),DATA3_FE,'g', 'linewidth',lw);

%plot(DATA4_t./(365*1440),DATA4_FE,'b', 'linewidth',lw);
%plot(DATA5_t./(365*1440),DATA5_FE,'r', 'linewidth',lw);


%xlabel('time (year)','fontsize', fs, 'fontname',fn);
xlabel('{Insulin} ({\mu}U/ml)','fontsize', fs, 'fontname',fn);
ylabel('{Induced Effect}','fontsize', fs, 'fontname',fn);


%%%%%%% 10% of the length of y axis

text(0,1.32,'F','fontsize',fpan,'fontweight','bold');
axis ([5 300 0 1.2])
%axis ([0 5 0 1.2])




%% End of plot



