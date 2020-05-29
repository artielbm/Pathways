
%%%%% simulations with the original model

clear

th_G60=155;
th_G120=140;
lw3=2;
lw4=1;
total_t=2628000; %paper2
OGTT_period=7200;
nPeriods=total_t/OGTT_period;

%%%% NOTE THAT odeparams.HGP_no_si=1 and  HGP_no_si=1 are used only for HGP test with or
%%%% w/o si dependency, otherwise ALWAYS 0
odeparams.HGP_no_si=0;
HGP_no_si=0; % for recomputing HGP in full_run.m

%%%% unit conversion from day to min in G equation, which is used to calculate glucose disposal
unit_con=0.0006944;

odeparams.Eg0=0.0118;%%%% :reduced Eg0 Default:0.0118
odeparams.BW=75;
odeparams.mealbar=11.055; 
outfile='Eg0_FIG2_control.xlsx';
odeparams.k=0.4861;

odeparams.meal=1;
odeparams.OGTT=0;
odeparams.IVGTT=0;
odeparams.r20=0.006;

odeparams.tar_si=0.1; %FIG1
odeparams.tau_si=360000; 
odeparams.tar_hepasi=0.85;%FIG1
odeparams.tau_hepasi=360000;
odeparams.Gs=90;


%%%%% Glucose disposal caculation
Eg0_control=0.0118;
Eg0_half=0.0059;

%%%% unit conversion from day to min in G equation, which is used to calculate glucose disposal


odeparams.ISRI_bar=0.5259; % beta-function defect for all pathways, default:1.4


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
    
    G30(i)=OGTT(i).y(tspanOGTT==30,1);
    I30(i)=OGTT(i).y(tspanOGTT==30,2);
    
    G60(i)=OGTT(i).y(tspanOGTT==60,1);
    I60(i)=OGTT(i).y(tspanOGTT==60,2);
    
    G90(i)=OGTT(i).y(tspanOGTT==90,1);
    I90(i)=OGTT(i).y(tspanOGTT==90,2);
    
    G120(i)=OGTT(i).y(tspanOGTT==120,1);
    I120(i)=OGTT(i).y(tspanOGTT==120,2);
    
    %%%%% si_d=si*I*G, Eg0_d=Eg0*G
    si_d0(i)=unit_con.*OGTT(i).y(tspanOGTT ==0,6).*OGTT(i).y(tspanOGTT==0,2).*OGTT(i).y(tspanOGTT==0,1);
    Eg0_d0(i)=Eg0_control.*OGTT(i).y(tspanOGTT==0,1);
    TGD0(i)=si_d0(i)+Eg0_d0(i);

    si_d60(i)=unit_con.*OGTT(i).y(tspanOGTT ==60,6).*OGTT(i).y(tspanOGTT==60,2).*OGTT(i).y(tspanOGTT==60,1);
    Eg0_d60(i)=Eg0_control.*OGTT(i).y(tspanOGTT==60,1);
    TGD60(i)=si_d60(i)+Eg0_d60(i);

    si_d120(i)=unit_con.*OGTT(i).y(tspanOGTT ==120,6).*OGTT(i).y(tspanOGTT==120,2).*OGTT(i).y(tspanOGTT==120,1);
    Eg0_d120(i)=Eg0_control.*OGTT(i).y(tspanOGTT==120,1);
    TGD120(i)=si_d120(i)+Eg0_d120(i);
    
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


%%%%%%% end of recompute HGP

end

%%% write data

long_OGTT=[T0(:),G0(:),G60(:),G120(:),I0(:),I60(:),I120(:),HGP(:),b(:),gamma(:),sigma(:),si(:),hepasi(:),si_d0(:),Eg0_d0(:),TGD0(:),si_d60(:),Eg0_d60(:),TGD60(:),si_d120(:),Eg0_d120(:),TGD120(:)];
tmp_long_OGTT=array2table(long_OGTT,'VariableNAMES',{'t','G0','G60','G120','I0','I60','I120','HGP','b','gamma','sigma','si','hepasi','si_d0','Eg0_d0','TGD0','si_d60','Eg0_d60','TGD60','si_d120','Eg0_d120','TGD120'});
writetable(tmp_long_OGTT,outfile) 




%%%%%%%%%%%%%%%% End of Eg0 FIG1 control 

clear

th_G60=155;
th_G120=140;
lw3=2;
lw4=1;
total_t=2628000; %paper2
OGTT_period=7200;
nPeriods=total_t/OGTT_period;

%%%% NOTE THAT odeparams.HGP_no_si=1 and  HGP_no_si=1 are used only for HGP test with or
%%%% w/o si dependency, otherwise ALWAYS 0
odeparams.HGP_no_si=0;
HGP_no_si=0; % for recomputing HGP in full_run.m

%%%% unit conversion from day to min in G equation, which is used to calculate glucose disposal
unit_con=0.0006944;

odeparams.Eg0=0.0059; %%%% :reduced Eg0 Default:0.0118
odeparams.BW=75;
odeparams.mealbar=11.055; 
outfile='Eg0_FIG2.xlsx';
odeparams.k=0.4861;

odeparams.meal=1;
odeparams.OGTT=0;
odeparams.IVGTT=0;
odeparams.r20=0.006;

odeparams.tar_si=0.1; %FIG1
odeparams.tau_si=360000; 
odeparams.tar_hepasi=0.85;%FIG1
odeparams.tau_hepasi=360000;
odeparams.Gs=90;


%%%%% Glucose disposal caculation

Eg0_control=0.0118;
Eg0_half=0.0059;


odeparams.ISRI_bar=0.5259; % beta-function defect for all pathways, default:1.4


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
    
    G30(i)=OGTT(i).y(tspanOGTT==30,1);
    I30(i)=OGTT(i).y(tspanOGTT==30,2);
    
    G60(i)=OGTT(i).y(tspanOGTT==60,1);
    I60(i)=OGTT(i).y(tspanOGTT==60,2);
    
    G90(i)=OGTT(i).y(tspanOGTT==90,1);
    I90(i)=OGTT(i).y(tspanOGTT==90,2);
    
    G120(i)=OGTT(i).y(tspanOGTT==120,1);
    I120(i)=OGTT(i).y(tspanOGTT==120,2);
    
    %%%%% si_d=si*I*G, Eg0_d=Eg0*G
   
    si_d0(i)=unit_con.*OGTT(i).y(tspanOGTT ==0,6).*OGTT(i).y(tspanOGTT==0,2).*OGTT(i).y(tspanOGTT==0,1);
    Eg0_d0(i)=Eg0_half.*OGTT(i).y(tspanOGTT==0,1);
    TGD0(i)=si_d0(i)+Eg0_d0(i);

    si_d60(i)=unit_con.*OGTT(i).y(tspanOGTT ==60,6).*OGTT(i).y(tspanOGTT==60,2).*OGTT(i).y(tspanOGTT==60,1);
    Eg0_d60(i)=Eg0_half.*OGTT(i).y(tspanOGTT==60,1);
    TGD60(i)=si_d60(i)+Eg0_d60(i);

    si_d120(i)=unit_con.*OGTT(i).y(tspanOGTT ==120,6).*OGTT(i).y(tspanOGTT==120,2).*OGTT(i).y(tspanOGTT==120,1);
    Eg0_d120(i)=Eg0_half.*OGTT(i).y(tspanOGTT==120,1);
    TGD120(i)=si_d120(i)+Eg0_d120(i);


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


%%%%%%% end of recompute HGP

end

%%% write data

long_OGTT=[T0(:),G0(:),G60(:),G120(:),I0(:),I60(:),I120(:),HGP(:),b(:),gamma(:),sigma(:),si(:),hepasi(:),si_d0(:),Eg0_d0(:),TGD0(:),si_d60(:),Eg0_d60(:),TGD60(:),si_d120(:),Eg0_d120(:),TGD120(:)];
tmp_long_OGTT=array2table(long_OGTT,'VariableNAMES',{'t','G0','G60','G120','I0','I60','I120','HGP','b','gamma','sigma','si','hepasi','si_d0','Eg0_d0','TGD0','si_d60','Eg0_d60','TGD60','si_d120','Eg0_d120','TGD120'});
writetable(tmp_long_OGTT,outfile) 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% plot data
%DATA=readtable('Eg0_FIG2_control.xlsx');
%DATA2=readtable('Eg0_FIG2.xlsx');
fs=12;
fn='arial';


figure(1); clf

%%%%% Glucose disposal caculation


Eg0_control=0.0118;
Eg0_half=0.0059;


%%% load longitudinal simulation data %%%%
%load Eg0.dat; 
Eg0=readtable('Eg0_FIG2_control.xlsx');

Eg0_t=Eg0{:,1};
Eg0_FG=Eg0{:,2};
Eg0_2hG=Eg0{:,4};
Eg0_FI=Eg0{:,5};
Eg0_2hI=Eg0{:,7};

Eg0_si_d0=Eg0{:,14};
Eg0_Eg0_d0=Eg0{:,15};
Eg0_TGD0=Eg0{:,16};

Eg0_si_d60=Eg0{:,17};
Eg0_Eg0_d60=Eg0{:,18};
Eg0_TGD60=Eg0{:,19};

Eg0_si_d120=Eg0{:,20};
Eg0_Eg0_d120=Eg0{:,21};
Eg0_TGD120=Eg0{:,22};


Eg0_si=Eg0{:,12};
Eg0_HGP=Eg0{:,13};
Eg0_sigma=Eg0{:,11};
Eg0_b=Eg0{:,9};
%%% Note that t should be rescaled by 365*1440 


%load Eg0_FIG2 data and name Eg0_2; 
Eg0_2=readtable('Eg0_FIG2.xlsx');
Eg0_2_t=Eg0_2{:,1};
Eg0_2_FG=Eg0_2{:,2};
Eg0_2_2hG=Eg0_2{:,4};
Eg0_2_FI=Eg0_2{:,5};
Eg0_2_2hI=Eg0_2{:,7};



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Eg0_2_si_d0=Eg0_2{:,14};
Eg0_2_Eg0_d0=Eg0_2{:,15};
Eg0_2_TGD0=Eg0_2{:,16};

Eg0_2_si_d60=Eg0_2{:,17};
Eg0_2_Eg0_d60=Eg0_2{:,18};
Eg0_2_TGD60=Eg0_2{:,19};

Eg0_2_si_d120=Eg0_2{:,20};
Eg0_2_Eg0_d120=Eg0_2{:,21};
Eg0_2_TGD120=Eg0_2{:,22};




Eg0_2_si=Eg0_2{:,12};
Eg0_2_HGP=Eg0_2{:,13};
Eg0_2_sigma=Eg0_2{:,11};
Eg0_2_b=Eg0_2{:,9};



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fs=12;
fs2=8;
fn='arial';
lw=1.5;
lw2=1.5;
lw3=0.5;
%%%%%%%%


fpan=10;
%title('FIG. 1', 'linewidth',15);

%%%%%%%%%%%%%%%%%%%%%

t_Sg=subplot(3,2,1);


plot([0 5],[Eg0_control Eg0_control],'k','linewidth',lw);
hold('on')
plot([0 5],[Eg0_half Eg0_half],'k--','linewidth',lw);


lh=legend('Control','Decreased Eg0','location','northwest');
set(lh,'FontSize',7.5); 
legend('boxoff');


%xlabel('time (year)','fontsize', fs, 'fontname',fn);
ylabel('Eg0 (1/min)','fontsize', fs, 'fontname',fn);

xticks([0 1 2 3 4 5])
xticklabels({'0','1','2','3','4','5'})
%%%%%%%%%%%%%%%%%%%%%% 10% of the length of y axis

text(0,0.044,'A','fontsize',fpan,'fontweight','bold');
axis ([0 5 0 0.04])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 10%  of the length of y axis

tsigma=subplot(3,2,2);

plot(Eg0_t./(365*1440),Eg0_sigma,'k', 'linewidth',lw);
hold('on')
plot(Eg0_2_t./(365*1440),Eg0_2_sigma,'k--', 'linewidth',lw);


%xlabel('time (year)','fontsize', fs, 'fontname',fn);
ylabel('{\sigma}','fontsize', fs, 'fontname',fn);

xticks([0 1 2 3 4 5])
xticklabels({'0','1','2','3','4','5'})
%%%%%%%%%%%%%%%%% 10% of the length of y axis

text(0,1.1,'B','fontsize',fpan,'fontweight','bold');
axis ([0 5 0 1])


%%%%%%%%%%%%%%%%%%%%%%%%

tG=subplot(3,2,3);


plot(Eg0_t./(365*1440),Eg0_FG,'k', 'linewidth',lw);

hold('on')
plot(Eg0_t./(365*1440),Eg0_2hG,'k', 'linewidth',lw);


plot(Eg0_2_t./(365*1440),Eg0_2_FG,'k--', 'linewidth',lw);
plot(Eg0_2_t./(365*1440),Eg0_2_2hG,'k--', 'linewidth',lw);


%xlabel('time (year)','fontsize', fs, 'fontname',fn);
ylabel('{Glucose} (mg/dl)','fontsize', fs, 'fontname',fn);

xticks([0 1 2 3 4 5])
xticklabels({'0','1','2','3','4','5'})
%%%%%%%%% 10% of %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%the length of y axis

text(0,550,'C','fontsize',fpan,'fontweight','bold');
text(2.2,180,'2hPG', 'fontsize',fs2,'fontname',fn);
text(1.5,75, 'FPG','fontsize',fs2,'fontname',fn);

fs10=8;

axis ([0.02 5 0 500])




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tI=subplot(3,2,4);
plot(Eg0_t./(365*1440),Eg0_FI,'k', 'linewidth',lw);

hold('on')
plot(Eg0_t./(365*1440),Eg0_2hI,'k', 'linewidth',lw);

plot(Eg0_2_t./(365*1440),Eg0_2_FI,'k--', 'linewidth',lw);
plot(Eg0_2_t./(365*1440),Eg0_2_2hI,'k--', 'linewidth',lw);




%xlabel('time (year)','fontsize', fs, 'fontname',fn);
ylabel('{Insulin} ({\mu}U/ml)','fontsize', fs, 'fontname',fn);

xticks([0 1 2 3 4 5])
xticklabels({'0','1','2','3','4','5'})

text(0,192.5,'D','fontsize',fpan,'fontweight','bold');
text(3,125,'2hPI', 'fontsize',fs2,'fontname',fn);
text(3.5,15, 'FPI','fontsize',fs2,'fontname',fn);


axis ([0 5 0 175])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


t_BGD=subplot(3,2,5);

plot(Eg0_t./(365*1440),Eg0_Eg0_d0,'k', 'linewidth',lw);
hold('on')
%plot(Eg0_t./(365*1440),(Eg0_si_d0 + Eg0_Eg0_d0),'k', 'linewidth',lw);
plot(Eg0_t./(365*1440),Eg0_TGD0,'r', 'linewidth',lw);

plot(Eg0_2_t./(365*1440),Eg0_2_Eg0_d0,'k--', 'linewidth',lw);
%plot(Eg0_2_t./(365*1440),(Eg0_2_si_d0 + Eg0_2_Eg0_d0),'k--', 'linewidth',lw);
plot(Eg0_2_t./(365*1440),Eg0_2_TGD0,'r--', 'linewidth',lw);

lh=legend('Control GE','Control TGD','Decreased Eg0 GE','Decreased Eg0 TGD','location','northwest');
set(lh,'FontSize',7.5); 
legend('boxoff');

xlabel('time (year)','fontsize', fs, 'fontname',fn);
ylabel('Basal GD(mg/dl/min)','fontsize', fs, 'fontname',fn);

xticks([0 1 2 3 4 5])
xticklabels({'0','1','2','3','4','5'})

%text(2.5,4, 'TGD','fontsize',fs2,'fontname',fn);
%text(2.7,2, 'GE','fontsize',fs2,'fontname',fn);
%%%%%%%%%%%%%%%%% 10% of the length of y axis

text(0,4.4,'E','fontsize',fpan,'fontweight','bold');
axis ([0.02 5 0 4])


%%%%%%%%%%%%%%%%%%%%%%%% 10%  of the length of y axis

t_PGD=subplot(3,2,6);

 %%%% G60
%plot(Eg0_t./(365*1440),Eg0_Eg0_d60,'k', 'linewidth',lw);
%hold('on')
%plot(Eg0_t./(365*1440),(Eg0_si_d60 + Eg0_Eg0_d60),'k', 'linewidth',lw);
%plot(Eg0_t./(365*1440),Eg0_TGD60,'k', 'linewidth',lw);

%plot(Eg0_2_t./(365*1440),Eg0_2_Eg0_d60,'k--', 'linewidth',lw);
%plot(Eg0_2_t./(365*1440),(Eg0_2_si_d60 + Eg0_2_Eg0_d60),'k--', 'linewidth',lw);
%plot(Eg0_t./(365*1440),Eg0_2_TGD60,'k', 'linewidth',lw);

%%%%% G120
plot(Eg0_t./(365*1440),Eg0_Eg0_d120,'k', 'linewidth',lw);
hold('on')
%plot(Eg0_t./(365*1440),(Eg0_si_d120 + Eg0_Eg0_d120),'k', 'linewidth',lw);
plot(Eg0_t./(365*1440),Eg0_TGD120,'r', 'linewidth',lw);

plot(Eg0_2_t./(365*1440),Eg0_2_Eg0_d120,'k--', 'linewidth',lw);
%plot(Eg0_2_t./(365*1440),(Eg0_2_si_d120 + Eg0_2_Eg0_d120),'k--', 'linewidth',lw);
plot(Eg0_t./(365*1440),Eg0_2_TGD120,'r--', 'linewidth',lw);

lh=legend('Control GE','Control TGD','Decreased Eg0 GE','Decreased Eg0 TGD','location','northwest');
set(lh,'FontSize',7.5); 
legend('boxoff');


xlabel('time (year)','fontsize', fs, 'fontname',fn);
ylabel('2h GD(mg/dl/min)','fontsize', fs, 'fontname',fn);


xticks([0 1 2 3 4 5])
xticklabels({'0','1','2','3','4','5'})

%text(1.3,6,'Total G Disposal','fontsize',fs2,'fontname',fn);
%text(1.3,3,'G Effectiveness','fontsize',fs2,'fontname',fn);


%%%%%%%%%%%%%%%%% 10% of the length of y axis

text(0,8.8,'F','fontsize',fpan,'fontweight','bold');
axis ([0.01 5 0 8])



%% End of plot

