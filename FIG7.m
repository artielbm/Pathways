clear
%%%% Longtudinal IGT IVGTTs
th_G60=155;
th_G120=140;
lw3=2;
lw4=1;
total_t=2877120; %paper2
%% For FIG 6, NGT at t=0, IGT at t=959040, CGI at t=1598400, T2D at t=2237760

IVGTT_period=319680;
nPeriods=total_t/IVGTT_period;

%%%% NOTE THAT odeparams.HGP_no_si=1 and  HGP_no_si=1 are used only for HGP test with or
%%%% w/o si dependency, otherwise ALWAYS 0
odeparams.HGP_no_si=0;
HGP_no_si=0; % for recomputing HGP in full_run.m

odeparams.BW=75;
odeparams.mealbar=11.055; 
outfile='FIG7_IGT.xlsx';
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
t0=0;
T=[]; %pre-allocate this array instead
Y=[];

tspan=[0:dt:IVGTT_period];
tspanIVGTT=[0:1:120];

IVGTT(1+nPeriods).t=[];
IVGTT(1+nPeriods).y=[];


%%%% DO an OGTT at baseline

%IVGTT at the end of longitudinal simulations
    odeparams.meal=0;
    odeparams.IVGTT=1;
    [t1,y1]=ode15s(@pathway,tspanIVGTT,init,options,odeparams);
    IVGTT(1).t=t1;
    IVGTT(1).y=y1;

%%%%%%%%%%%%%
%%%% Begin of longitudinal simulations and OGTTs

for i=1:nPeriods

   %inter-IVGTTs
    odeparams.meal=1;
    odeparams.IVGTT=0;
    odeparams.curt=(i-1)*IVGTT_period; % show current time for forced improvement of si, tar_si2 and tau_si2
    
    [t,y]=ode15s(@pathway,tspan,init,options,odeparams);


    T=[T;t+(i-1)*IVGTT_period];
   
    Y=[Y;y];
    init=y(end,:);
    
    i %remove this to speed up
    
    %OGTT at the end of longitudinal simulations
    odeparams.meal=0;
    odeparams.IVGTT=1;
    [t,y]=ode15s(@pathway,tspanIVGTT,init,options,odeparams);
    IVGTT(i+1).t=t;
    IVGTT(i+1).y=y;
       
end

%%%%% End of Simulations%%%

T0(1)=0;
for i=2:nPeriods + 1  

T0(i)=T(i-1+(i-1)*IVGTT_period/dt); %%%% Don't be confused with # of array and real time

end 
%%
for i=1:nPeriods + 1
    
    
    G0(i)=IVGTT(i).y(tspanIVGTT==0,1);
    I0(i)=IVGTT(i).y(tspanIVGTT==0,2);
    
    G1(i)=IVGTT(i).y(tspanIVGTT==1,1);
    I1(i)=IVGTT(i).y(tspanIVGTT==1,2);

    G2(i)=IVGTT(i).y(tspanIVGTT==2,1);
    I2(i)=IVGTT(i).y(tspanIVGTT==2,2);

    G3(i)=IVGTT(i).y(tspanIVGTT==3,1);
    I3(i)=IVGTT(i).y(tspanIVGTT==3,2);

    G4(i)=IVGTT(i).y(tspanIVGTT==4,1);
    I4(i)=IVGTT(i).y(tspanIVGTT==4,2);

    G5(i)=IVGTT(i).y(tspanIVGTT==5,1);
    I5(i)=IVGTT(i).y(tspanIVGTT==5,2);
    

    G6(i)=IVGTT(i).y(tspanIVGTT==6,1);
    I6(i)=IVGTT(i).y(tspanIVGTT==6,2);
    
    G7(i)=IVGTT(i).y(tspanIVGTT==7,1);
    I7(i)=IVGTT(i).y(tspanIVGTT==7,2);

    G8(i)=IVGTT(i).y(tspanIVGTT==8,1);
    I8(i)=IVGTT(i).y(tspanIVGTT==8,2);

    G9(i)=IVGTT(i).y(tspanIVGTT==9,1);
    I9(i)=IVGTT(i).y(tspanIVGTT==9,2);
    
    G10(i)=IVGTT(i).y(tspanIVGTT==10,1);
    I10(i)=IVGTT(i).y(tspanIVGTT==10,2);

    G11(i)=IVGTT(i).y(tspanIVGTT==11,1);
    I11(i)=IVGTT(i).y(tspanIVGTT==11,2);

    G12(i)=IVGTT(i).y(tspanIVGTT==12,1);
    I12(i)=IVGTT(i).y(tspanIVGTT==12,2);
    
    G13(i)=IVGTT(i).y(tspanIVGTT==13,1);
    I13(i)=IVGTT(i).y(tspanIVGTT==13,2);

    G14(i)=IVGTT(i).y(tspanIVGTT==14,1);
    I14(i)=IVGTT(i).y(tspanIVGTT==14,2);


    G15(i)=IVGTT(i).y(tspanIVGTT==15,1);
    I15(i)=IVGTT(i).y(tspanIVGTT==15,2);


    G20(i)=IVGTT(i).y(tspanIVGTT==20,1);
    I20(i)=IVGTT(i).y(tspanIVGTT==20,2); 

    G25(i)=IVGTT(i).y(tspanIVGTT==25,1);
    I25(i)=IVGTT(i).y(tspanIVGTT==25,2); 

    G30(i)=IVGTT(i).y(tspanIVGTT==30,1);
    I30(i)=IVGTT(i).y(tspanIVGTT==30,2); 

    G60(i)=IVGTT(i).y(tspanIVGTT==60,1);
    I60(i)=IVGTT(i).y(tspanIVGTT==60,2); 

    G90(i)=IVGTT(i).y(tspanIVGTT==90,1);
    I90(i)=IVGTT(i).y(tspanIVGTT==90,2); 

    G120(i)=IVGTT(i).y(tspanIVGTT==120,1);
    I120(i)=IVGTT(i).y(tspanIVGTT==120,2); 


    b(i)=IVGTT(i).y(tspanIVGTT ==0,3);
    gamma(i)=IVGTT(i).y(tspanIVGTT ==0,4);
    sigma(i)=IVGTT(i).y(tspanIVGTT ==0,5);
    si(i)=IVGTT(i).y(tspanIVGTT ==0,6);
    hepasi(i)=IVGTT(i).y(tspanIVGTT ==0,7);   
    RRP(i)=IVGTT(i).y(tspanIVGTT ==0,8);

    m_G(i)=mean([G0(i),G30(i),G60(i),G90(i),G120(i)]);
    m_I(i)=mean([I0(i),I30(i),I60(i),I90(i),I120(i)]);
    
    IGI(i)=(I30(i)-I0(i))/(G30(i) - G0(i));
    matsuda(i)=10000/(sqrt(G0(i)*I0(i)*m_G(i)*m_I(i)));
    
   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%% recompute HGP

    % hepa_bar=15.443; hepa_k=0.27; hepa_b=-3.54277; con_si=0.8;
    % hepa_max= hepa_bar./(hepa_k +IVGTT(i).y(tspanIVGTT==0,6).*(1-HGP_no_si) + con_si*HGP_no_si) + hepa_b;

     %alpha_max=6; alpha_k=0.4; alpha_b=-0.5;
     %alpha_HGP= alpha_max./(alpha_k + IVGTT(i).y(tspanIVGTT==0,6).*(1-HGP_no_si) + con_si*HGP_no_si) + alpha_b;

     %HGP_b=0.104166;  
     %HGP(i) = hepa_max./(alpha_HGP + IVGTT(i).y(tspanIVGTT==0,2)*IVGTT(i).y(tspanIVGTT==0,7)) + HGP_b;


%% end of recompute HGP

end

%%% write data

long_IVGTT=[T0(:),G0(:),G1(:),G2(:),G3(:),G4(:),G5(:),G6(:),G7(:),G8(:),G9(:),G10(:),G11(:),G12(:),G13(:),G14(:),G15(:),G20(:),G25(:),G30(:),G60(:),G90(:),G120(:),I0(:),I1(:),I2(:),I3(:),I4(:),I5(:),I6(:),I7(:),I8(:),I9(:),I10(:),I11(:),I12(:),I13(:),I14(:),I15(:),I20(:),I25(:),I30(:),I60(:),I90(:),I120(:),RRP(:)];
tmp_long_IVGTT=array2table(long_IVGTT,'VariableNAME',{'t','G0','G1','G2','G3','G4','G5','G6','G7','G8','G9','G10','G11','G12','G13','G14','G15','G20','G25','G30','G60','G90','G120','I0','I1','I2','I3','I4','I5','I6','I7','I8','I9','I10','I11','I12','I13','I14','I15','I20','I25','I30','I60','I90','I120','RRP'});
writetable(tmp_long_IVGTT,outfile) 

%% read saved data
DATA=readtable('FIG7_IGT.xlsx');

%%%% END of Longtudinal IGT IVGTTs







%%%% Longtudinal IFG IVGTTs
th_G60=155;
th_G120=140;
lw3=2;
lw4=1;
total_t=2877120; %paper2
%% For FIG 6, NGT at t=0, IGT at t=959040, CGI at t=1598400, T2D at t=2237760

IVGTT_period=319680;
nPeriods=total_t/IVGTT_period;

%%%% NOTE THAT odeparams.HGP_no_si=1 and  HGP_no_si=1 are used only for HGP test with or
%%%% w/o si dependency, otherwise ALWAYS 0
odeparams.HGP_no_si=0;
HGP_no_si=0; % for recomputing HGP in full_run.m

odeparams.BW=75;
odeparams.mealbar=11.055; 
outfile='FIG7_IFG.xlsx';
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
t0=0;
T=[]; %pre-allocate this array instead
Y=[];

tspan=[0:dt:IVGTT_period];
tspanIVGTT=[0:1:120];

IVGTT(1+nPeriods).t=[];
IVGTT(1+nPeriods).y=[];


%%%% DO an OGTT at baseline

%IVGTT at the end of longitudinal simulations
    odeparams.meal=0;
    odeparams.IVGTT=1;
    [t1,y1]=ode15s(@pathway,tspanIVGTT,init,options,odeparams);
    IVGTT(1).t=t1;
    IVGTT(1).y=y1;

%%%%%%%%%%%%%
%%%% Begin of longitudinal simulations and OGTTs

for i=1:nPeriods

   %inter-IVGTTs
    odeparams.meal=1;
    odeparams.IVGTT=0;
    odeparams.curt=(i-1)*IVGTT_period; % show current time for forced improvement of si, tar_si2 and tau_si2
    
    [t,y]=ode15s(@pathway,tspan,init,options,odeparams);


    T=[T;t+(i-1)*IVGTT_period];
   
    Y=[Y;y];
    init=y(end,:);
    
    i %remove this to speed up
    
    %OGTT at the end of longitudinal simulations
    odeparams.meal=0;
    odeparams.IVGTT=1;
    [t,y]=ode15s(@pathway,tspanIVGTT,init,options,odeparams);
    IVGTT(i+1).t=t;
    IVGTT(i+1).y=y;
       
end

%%%%% End of Simulations%%%

T0(1)=0;
for i=2:nPeriods + 1  

T0(i)=T(i-1+(i-1)*IVGTT_period/dt); %%%% Don't be confused with # of array and real time

end 
%%
for i=1:nPeriods + 1
    
    
    G0(i)=IVGTT(i).y(tspanIVGTT==0,1);
    I0(i)=IVGTT(i).y(tspanIVGTT==0,2);
    
    G1(i)=IVGTT(i).y(tspanIVGTT==1,1);
    I1(i)=IVGTT(i).y(tspanIVGTT==1,2);

    G2(i)=IVGTT(i).y(tspanIVGTT==2,1);
    I2(i)=IVGTT(i).y(tspanIVGTT==2,2);

    G3(i)=IVGTT(i).y(tspanIVGTT==3,1);
    I3(i)=IVGTT(i).y(tspanIVGTT==3,2);

    G4(i)=IVGTT(i).y(tspanIVGTT==4,1);
    I4(i)=IVGTT(i).y(tspanIVGTT==4,2);

    G5(i)=IVGTT(i).y(tspanIVGTT==5,1);
    I5(i)=IVGTT(i).y(tspanIVGTT==5,2);
    

    G6(i)=IVGTT(i).y(tspanIVGTT==6,1);
    I6(i)=IVGTT(i).y(tspanIVGTT==6,2);
    
    G7(i)=IVGTT(i).y(tspanIVGTT==7,1);
    I7(i)=IVGTT(i).y(tspanIVGTT==7,2);

    G8(i)=IVGTT(i).y(tspanIVGTT==8,1);
    I8(i)=IVGTT(i).y(tspanIVGTT==8,2);

    G9(i)=IVGTT(i).y(tspanIVGTT==9,1);
    I9(i)=IVGTT(i).y(tspanIVGTT==9,2);
    
    G10(i)=IVGTT(i).y(tspanIVGTT==10,1);
    I10(i)=IVGTT(i).y(tspanIVGTT==10,2);

    G11(i)=IVGTT(i).y(tspanIVGTT==11,1);
    I11(i)=IVGTT(i).y(tspanIVGTT==11,2);

    G12(i)=IVGTT(i).y(tspanIVGTT==12,1);
    I12(i)=IVGTT(i).y(tspanIVGTT==12,2);
    
    G13(i)=IVGTT(i).y(tspanIVGTT==13,1);
    I13(i)=IVGTT(i).y(tspanIVGTT==13,2);

    G14(i)=IVGTT(i).y(tspanIVGTT==14,1);
    I14(i)=IVGTT(i).y(tspanIVGTT==14,2);


    G15(i)=IVGTT(i).y(tspanIVGTT==15,1);
    I15(i)=IVGTT(i).y(tspanIVGTT==15,2);


    G20(i)=IVGTT(i).y(tspanIVGTT==20,1);
    I20(i)=IVGTT(i).y(tspanIVGTT==20,2); 

    G25(i)=IVGTT(i).y(tspanIVGTT==25,1);
    I25(i)=IVGTT(i).y(tspanIVGTT==25,2); 

    G30(i)=IVGTT(i).y(tspanIVGTT==30,1);
    I30(i)=IVGTT(i).y(tspanIVGTT==30,2); 

    G60(i)=IVGTT(i).y(tspanIVGTT==60,1);
    I60(i)=IVGTT(i).y(tspanIVGTT==60,2); 

    G90(i)=IVGTT(i).y(tspanIVGTT==90,1);
    I90(i)=IVGTT(i).y(tspanIVGTT==90,2); 

    G120(i)=IVGTT(i).y(tspanIVGTT==120,1);
    I120(i)=IVGTT(i).y(tspanIVGTT==120,2); 


    b(i)=IVGTT(i).y(tspanIVGTT ==0,3);
    gamma(i)=IVGTT(i).y(tspanIVGTT ==0,4);
    sigma(i)=IVGTT(i).y(tspanIVGTT ==0,5);
    si(i)=IVGTT(i).y(tspanIVGTT ==0,6);
    hepasi(i)=IVGTT(i).y(tspanIVGTT ==0,7);   
    RRP(i)=IVGTT(i).y(tspanIVGTT ==0,8);

    m_G(i)=mean([G0(i),G30(i),G60(i),G90(i),G120(i)]);
    m_I(i)=mean([I0(i),I30(i),I60(i),I90(i),I120(i)]);
    
    IGI(i)=(I30(i)-I0(i))/(G30(i) - G0(i));
    matsuda(i)=10000/(sqrt(G0(i)*I0(i)*m_G(i)*m_I(i)));
    


end

%%% write data

long_IVGTT=[T0(:),G0(:),G1(:),G2(:),G3(:),G4(:),G5(:),G6(:),G7(:),G8(:),G9(:),G10(:),G11(:),G12(:),G13(:),G14(:),G15(:),G20(:),G25(:),G30(:),G60(:),G90(:),G120(:),I0(:),I1(:),I2(:),I3(:),I4(:),I5(:),I6(:),I7(:),I8(:),I9(:),I10(:),I11(:),I12(:),I13(:),I14(:),I15(:),I20(:),I25(:),I30(:),I60(:),I90(:),I120(:),RRP(:)];
tmp_long_IVGTT=array2table(long_IVGTT,'VariableNAME',{'t','G0','G1','G2','G3','G4','G5','G6','G7','G8','G9','G10','G11','G12','G13','G14','G15','G20','G25','G30','G60','G90','G120','I0','I1','I2','I3','I4','I5','I6','I7','I8','I9','I10','I11','I12','I13','I14','I15','I20','I25','I30','I60','I90','I120','RRP'});
writetable(tmp_long_IVGTT,outfile) 

%% read saved data
DATA2=readtable('FIG7_IFG.xlsx');

%%%% END of Longtudinal IFG IVGTTs

%% %% plot Longitudinal IVGTTs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%
fs=12;
fs2=8;
fn='arial';
lw=1.5;
lw2=1.5;
lw3=0.75;
%%%%%%%%
fs10=10;
fpan=10;
t_f=60;

%%% RRP vector
 IGT_RRP=[DATA{1,46},DATA{2,46},DATA{3,46},DATA{5,46},DATA{8,46}];
 IFG_RRP=[DATA2{1,46},DATA2{2,46},DATA2{4,46},DATA2{5,46},DATA2{8,46}];

%%% plot for bargraph
IGT_first=categorical({'NGT','NGT2','IGT','CGI','T2D'});
IGT_first = reordercats(IGT_first,{'NGT','NGT2','IGT','CGI','T2D'});


IFG_first=categorical({'NGT','NGT2','IFG','CGI','T2D'});
IFG_first = reordercats(IFG_first,{'NGT','NGT2','IFG','CGI','T2D'});

%%%%%%%
IVGTT_t=[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 20 25 30 60 90 120];

OGTT1_I=subplot(2,2,1);

plot(IVGTT_t,DATA{1,24:45},'k', 'linewidth',lw);
hold('on')
plot(IVGTT_t,DATA{2,24:45},'k', 'linewidth',lw3);
plot(IVGTT_t,DATA{3,24:45},'k:', 'linewidth',lw);
plot(IVGTT_t,DATA{5,24:45},'k--', 'linewidth',lw);
plot(IVGTT_t,DATA{8,24:45},'k-.', 'linewidth',lw);


xlabel('time (min)','fontsize', fs, 'fontname',fn);
ylabel('I (\muU/ml)','fontsize', fs, 'fontname',fn);

lh=legend('NGT1','NGT2','IGT','CGI','T2D','location','east');
set(lh,'FontSize',8); 
legend('boxoff');

%% 10% of the length of y axis

text(0,110,'A','fontsize',fpan,'fontweight','bold');
axis ([0 t_f 0 100])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OGTT2_I=subplot(2,2,2);
plot(IVGTT_t,DATA2{1,24:45},'k', 'linewidth',lw);
hold('on')
plot(IVGTT_t,DATA2{2,24:45},'k', 'linewidth',lw3);
plot(IVGTT_t,DATA2{3,24:45},'k:', 'linewidth',lw);
plot(IVGTT_t,DATA2{5,24:45},'k--', 'linewidth',lw);
plot(IVGTT_t,DATA2{8,24:45},'k-.', 'linewidth',lw);



xlabel('time (min)','fontsize', fs, 'fontname',fn);
ylabel('I (\muU/ml)','fontsize', fs, 'fontname',fn);


lh=legend('NGT1','NGT2','IFG','CGI','T2D','location','east');
set(lh,'FontSize',8); 
legend('boxoff');


%% 10% of the length of y axis


text(0,110,'B','fontsize',fpan,'fontweight','bold');
axis ([0 t_f 0 100])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OGTT1_RRP=subplot(2,2,3);

bar(IGT_first, IGT_RRP,'k')

%% 10% of the length of y axis

text(0.5,88,'C','fontsize',fpan,'fontweight','bold');
ylim([0 80])

ylabel('RRP','fontsize', fs, 'fontname',fn);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OGTT2_RRP=subplot(2,2,4);

%% 10% of the length of y axis
bar(IFG_first, IFG_RRP,'k')

%text(0,66,'D','fontsize',fpan,'fontweight','bold');
%axis ([0 t_f 0 100])

text(0.5,88,'D','fontsize',fpan,'fontweight','bold');
ylim([0 80])

ylabel('RRP','fontsize', fs, 'fontname',fn);






