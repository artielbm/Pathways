clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% longitudinal OGTTs of IGT first pathway %%%%%%%%%%%


th_G60=155;
th_G120=140;
lw3=2;
lw4=1;
total_t=2877120; %paper2
%% For FIG 6, NGT at t=0, IGT at t=959040, CGI at t=1598400, T2D at t=2237760

OGTT_period=319680;
nPeriods=total_t/OGTT_period;
doplot=[0 959040 1598400 2237760 ];

%%%% NOTE THAT odeparams.HGP_no_si=1 and  HGP_no_si=1 are used only for HGP test with or
%%%% w/o si dependency, otherwise ALWAYS 0
odeparams.HGP_no_si=0;
HGP_no_si=0; % for recomputing HGP in full_run.m

odeparams.BW=75;
odeparams.mealbar=11.055; 
outfile='FIG6_IGT.xlsx';
odeparams.k=0.4861;

odeparams.meal=1;
odeparams.OGTT=0;
odeparams.IVGTT=0;
odeparams.r20=0.006;

odeparams.tar_si=0.1; %FIG1
odeparams.tau_si=360000; 
odeparams.tar_hepasi=0.85;%FIG1
odeparams.tau_hepasi=360000;
odeparams.Gs=100;
odeparams.ISRI_bar=0.5259; % beta-function defect for all pathways


init=[78.59, 5.63, 1533.91, -0.07663, 1, 0.8, 1, 60.24, 443.39]; % IC   


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
    [t1,y1]=ode15s(@pathway,tspanOGTT,init,options,odeparams);
    OGTT(1).t=t1;
    OGTT(1).y=y1;

%%%%%%%%%%%%%
%%%% Begin of longitudinal simulations and OGTTs

for i=1:nPeriods

   %inter-OGTTs
    odeparams.meal=1;
    odeparams.OGTT=0;
    odeparams.curt=(i-1)*OGTT_period; % show current time for forced improvement of si, tar_si2 and tau_si2
    
    [t,y]=ode15s(@pathway,tspan,init,options,odeparams);


    T=[T;t+(i-1)*OGTT_period];
   
    Y=[Y;y];
    init=y(end,:);
    
    i %remove this to speed up
    
    %OGTT at the end of longitudinal simulations
    odeparams.meal=0;
    odeparams.OGTT=1;
    [t,y]=ode15s(@pathway,tspanOGTT,init,options,odeparams);
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
    
    G10(i)=OGTT(i).y(tspanOGTT==10,1);
    I10(i)=OGTT(i).y(tspanOGTT==10,2);

    G20(i)=OGTT(i).y(tspanOGTT==20,1);
    I20(i)=OGTT(i).y(tspanOGTT==20,2);

    G30(i)=OGTT(i).y(tspanOGTT==30,1);
    I30(i)=OGTT(i).y(tspanOGTT==30,2);

    G40(i)=OGTT(i).y(tspanOGTT==40,1);
    I40(i)=OGTT(i).y(tspanOGTT==40,2);

    G50(i)=OGTT(i).y(tspanOGTT==50,1);
    I50(i)=OGTT(i).y(tspanOGTT==50,2);
    

    G60(i)=OGTT(i).y(tspanOGTT==60,1);
    I60(i)=OGTT(i).y(tspanOGTT==60,2);
    
    G70(i)=OGTT(i).y(tspanOGTT==70,1);
    I70(i)=OGTT(i).y(tspanOGTT==70,2);

    G80(i)=OGTT(i).y(tspanOGTT==80,1);
    I80(i)=OGTT(i).y(tspanOGTT==80,2);

    G90(i)=OGTT(i).y(tspanOGTT==90,1);
    I90(i)=OGTT(i).y(tspanOGTT==90,2);
    
    G100(i)=OGTT(i).y(tspanOGTT==100,1);
    I100(i)=OGTT(i).y(tspanOGTT==100,2);

    G110(i)=OGTT(i).y(tspanOGTT==110,1);
    I110(i)=OGTT(i).y(tspanOGTT==110,2);

    G120(i)=OGTT(i).y(tspanOGTT==120,1);
    I120(i)=OGTT(i).y(tspanOGTT==120,2);
    
    b(i)=OGTT(i).y(tspanOGTT ==0,3);
    gamma(i)=OGTT(i).y(tspanOGTT ==0,4);
    sigma(i)=OGTT(i).y(tspanOGTT ==0,5);
    si(i)=OGTT(i).y(tspanOGTT ==0,6);
    hepasi(i)=OGTT(i).y(tspanOGTT ==0,7);   


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


%% end of recompute HGP

end

%%% write data

long_OGTT=[T0(:),G0(:),G10(:),G20(:),G30(:),G40(:),G50(:),G60(:),G70(:),G80(:),G90(:),G100(:),G110(:),G120(:),I0(:),I10(:),I20(:),I30(:),I40(:),I50(:),I60(:),I70(:),I80(:),I90(:),I100(:),I110(:),I120(:),HGP(:),b(:),gamma(:),sigma(:),si(:),hepasi(:)];
tmp_long_OGTT=array2table(long_OGTT,'VariableNAMES',{'t','G0','G10','G20','G30','G40','G50','G60','G70','G80','G90','G100','G110','G120','I0','I10','I20','I30','I40','I50','I60','I70','I80','I90','I100','I110','I120','HGP','b','gamma','sigma','si','hepasi'});
writetable(tmp_long_OGTT,outfile) 



DATA=readtable('FIG6_IGT.xlsx');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% End of longitudinal OGTTs of IGT first pathway %%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% longitudinal OGTTs of IFG first pathway %%%%%%%%%%%


th_G60=155;
th_G120=140;
lw3=2;
lw4=1;
total_t=2877120; %paper2
%% For FIG 6, NGT at t=0, IGT at t=959040, CGI at t=1598400, T2D at t=2237760

OGTT_period=319680;
nPeriods=total_t/OGTT_period;
doplot=[0 959040 1598400 2237760 ];

%%%% NOTE THAT odeparams.HGP_no_si=1 and  HGP_no_si=1 are used only for HGP test with or
%%%% w/o si dependency, otherwise ALWAYS 0
odeparams.HGP_no_si=0;
HGP_no_si=0; % for recomputing HGP in full_run.m

odeparams.BW=75;
odeparams.mealbar=11.055; 
outfile='FIG6_IFG.xlsx';
odeparams.k=0.4861;

odeparams.meal=1;
odeparams.OGTT=0;
odeparams.IVGTT=0;
odeparams.r20=0.006;

odeparams.tar_si=0.5; %FIG1
odeparams.tau_si=360000; 
odeparams.tar_hepasi=0.1;%FIG1
odeparams.tau_hepasi=360000;
odeparams.Gs=100;
odeparams.ISRI_bar=0.5259; % beta-function defect for all pathways


init=[78.59, 5.63, 1533.91, -0.07663, 1, 0.8, 1, 60.24, 443.39]; % IC   


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
    [t1,y1]=ode15s(@pathway,tspanOGTT,init,options,odeparams);
    OGTT(1).t=t1;
    OGTT(1).y=y1;

%%%%%%%%%%%%%
%%%% Begin of longitudinal simulations and OGTTs

for i=1:nPeriods

   %inter-OGTTs
    odeparams.meal=1;
    odeparams.OGTT=0;
    odeparams.curt=(i-1)*OGTT_period; % show current time for forced improvement of si, tar_si2 and tau_si2
    
    [t,y]=ode15s(@pathway,tspan,init,options,odeparams);


    T=[T;t+(i-1)*OGTT_period];
   
    Y=[Y;y];
    init=y(end,:);
    
    i %remove this to speed up
    
    %OGTT at the end of longitudinal simulations
    odeparams.meal=0;
    odeparams.OGTT=1;
    [t,y]=ode15s(@pathway,tspanOGTT,init,options,odeparams);
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
    
    G10(i)=OGTT(i).y(tspanOGTT==10,1);
    I10(i)=OGTT(i).y(tspanOGTT==10,2);

    G20(i)=OGTT(i).y(tspanOGTT==20,1);
    I20(i)=OGTT(i).y(tspanOGTT==20,2);

    G30(i)=OGTT(i).y(tspanOGTT==30,1);
    I30(i)=OGTT(i).y(tspanOGTT==30,2);

    G40(i)=OGTT(i).y(tspanOGTT==40,1);
    I40(i)=OGTT(i).y(tspanOGTT==40,2);

    G50(i)=OGTT(i).y(tspanOGTT==50,1);
    I50(i)=OGTT(i).y(tspanOGTT==50,2);
    

    G60(i)=OGTT(i).y(tspanOGTT==60,1);
    I60(i)=OGTT(i).y(tspanOGTT==60,2);
    
    G70(i)=OGTT(i).y(tspanOGTT==70,1);
    I70(i)=OGTT(i).y(tspanOGTT==70,2);

    G80(i)=OGTT(i).y(tspanOGTT==80,1);
    I80(i)=OGTT(i).y(tspanOGTT==80,2);

    G90(i)=OGTT(i).y(tspanOGTT==90,1);
    I90(i)=OGTT(i).y(tspanOGTT==90,2);
    
    G100(i)=OGTT(i).y(tspanOGTT==100,1);
    I100(i)=OGTT(i).y(tspanOGTT==100,2);

    G110(i)=OGTT(i).y(tspanOGTT==110,1);
    I110(i)=OGTT(i).y(tspanOGTT==110,2);

    G120(i)=OGTT(i).y(tspanOGTT==120,1);
    I120(i)=OGTT(i).y(tspanOGTT==120,2);
    
    b(i)=OGTT(i).y(tspanOGTT ==0,3);
    gamma(i)=OGTT(i).y(tspanOGTT ==0,4);
    sigma(i)=OGTT(i).y(tspanOGTT ==0,5);
    si(i)=OGTT(i).y(tspanOGTT ==0,6);
    hepasi(i)=OGTT(i).y(tspanOGTT ==0,7);   


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


%% end of recompute HGP

end

%%% write data

long_OGTT=[T0(:),G0(:),G10(:),G20(:),G30(:),G40(:),G50(:),G60(:),G70(:),G80(:),G90(:),G100(:),G110(:),G120(:),I0(:),I10(:),I20(:),I30(:),I40(:),I50(:),I60(:),I70(:),I80(:),I90(:),I100(:),I110(:),I120(:),HGP(:),b(:),gamma(:),sigma(:),si(:),hepasi(:)];
tmp_long_OGTT=array2table(long_OGTT,'VariableNAMES',{'t','G0','G10','G20','G30','G40','G50','G60','G70','G80','G90','G100','G110','G120','I0','I10','I20','I30','I40','I50','I60','I70','I80','I90','I100','I110','I120','HGP','b','gamma','sigma','si','hepasi'});
writetable(tmp_long_OGTT,outfile) 

DATA2=readtable('FIG6_IFG.xlsx');

%%%%% Endo of longitudinal OGTTs of IFG first pathway %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fs=12;
fs2=8;
fn='arial';
lw=1.5;
lw2=1.5;
lw3=0.5;
lw4=0.75;
%%%%%%%%

OGTT_dt=20;

fpan=10;
%%OGTT_t=[0:10:120];
OGTT_t=[0:OGTT_dt:120];

OGTT_gdt=0.1*OGTT_dt;
%%%%%%%


OGTT1_G=subplot(2,2,1);

plot(OGTT_t,DATA{1,2:OGTT_gdt:14},'k', 'linewidth',lw);
hold('on')
plot(OGTT_t,DATA{3,2:OGTT_gdt:14},'k:', 'linewidth',lw);
plot(OGTT_t,DATA{5,2:OGTT_gdt:14},'k--', 'linewidth',lw);
plot(OGTT_t,DATA{8,2:OGTT_gdt:14},'k-.', 'linewidth',lw);


xlabel('time (min)','fontsize', fs, 'fontname',fn);
ylabel('G (mg/dl)','fontsize', fs, 'fontname',fn);

lh=legend('NGT 1','IGT','CGI','T2D','location','south');
set(lh,'FontSize',7.5); 
legend('boxoff');

%% 10% of the length of y axis

text(0,363,'A','fontsize',fpan,'fontweight','bold');
axis ([0 120 0 330])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OGTT2_G=subplot(2,2,2);


plot(OGTT_t,DATA2{1,2:OGTT_gdt:14},'k', 'linewidth',lw);
hold('on')
%plot(OGTT_IGT_first_NGT2_t,OGTT_IGT_first_NGT2_G,'k','linewidth',lw4);
plot(OGTT_t,DATA2{3,2:OGTT_gdt:14},'k:', 'linewidth',lw);
plot(OGTT_t,DATA2{5,2:OGTT_gdt:14},'k--', 'linewidth',lw);
plot(OGTT_t,DATA2{8,2:OGTT_gdt:14},'k-.', 'linewidth',lw);

xlabel('time (min)','fontsize', fs, 'fontname',fn);
ylabel('G (mg/dl)','fontsize', fs, 'fontname',fn);


lh=legend('NGT 1','IFG','CGI','T2D','location','south');
set(lh,'FontSize',7.5); 
legend('boxoff');


%% 10% of the length of y axis


text(0,363,'B','fontsize',fpan,'fontweight','bold');
axis ([0 120 0 330])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OGTT1_I=subplot(2,2,3);


plot(OGTT_t,DATA{1,15:OGTT_gdt:27},'k', 'linewidth',lw);
hold('on')
%plot(OGTT_IGT_first_NGT2_t,OGTT_IGT_first_NGT2_G,'k','linewidth',lw4);
plot(OGTT_t,DATA{3,15:OGTT_gdt:27},'k:', 'linewidth',lw);
plot(OGTT_t,DATA{5,15:OGTT_gdt:27},'k--', 'linewidth',lw);
plot(OGTT_t,DATA{8,15:OGTT_gdt:27},'k-.', 'linewidth',lw);


xlabel('time (min)','fontsize', fs, 'fontname',fn);
ylabel('I (\muU/ml)','fontsize', fs, 'fontname',fn);



%% 10% of the length of y axis

text(0,220,'C','fontsize',fpan,'fontweight','bold');
axis ([0 120 0 200])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OGTT2_I=subplot(2,2,4);

plot(OGTT_t,DATA2{1,15:OGTT_gdt:27},'k', 'linewidth',lw);
hold('on')
%plot(OGTT_IGT_first_NGT2_t,OGTT_IGT_first_NGT2_G,'k','linewidth',lw4);
plot(OGTT_t,DATA2{3,15:OGTT_gdt:27},'k:', 'linewidth',lw);
plot(OGTT_t,DATA2{5,15:OGTT_gdt:27},'k--', 'linewidth',lw);
plot(OGTT_t,DATA2{8,15:OGTT_gdt:27},'k-.', 'linewidth',lw);

xlabel('time (min)','fontsize', fs, 'fontname',fn);
ylabel('I (\muU/ml)','fontsize', fs, 'fontname',fn);


%% 10% of the length of y axis


text(0,220,'D','fontsize',fpan,'fontweight','bold');
axis ([0 120 0 200])


