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

odeparams.BW=75;
odeparams.mealbar=11.055; 
outfile='FIG4.xlsx';
odeparams.k=0.4861;

odeparams.meal=1;
odeparams.OGTT=0;
odeparams.IVGTT=0;
odeparams.r20=0.006;

odeparams.tar_si=0.03; %FIG4
odeparams.tau_si=216000; %FIG4
odeparams.tar_hepasi=1.3;%FIG4
odeparams.tau_hepasi=1440; %FIG4
odeparams.Gs=100;
odeparams.ISRI_bar=0.5259; % beta-function defect, default:1.4


init=[78.59, 5.63, 1533.91, -0.07663, 1, 0.8, 1.3, 60.24, 443.39]; % IC   


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

long_OGTT=[T0(:),G0(:),G60(:),G120(:),I0(:),I60(:),I120(:),HGP(:),b(:),gamma(:),sigma(:),si(:),hepasi(:)];
tmp_long_OGTT=array2table(long_OGTT,'VariableNAMES',{'t','G0','G60','G120','I0','I60','I120','HGP','b','gamma','sigma','si','hepasi'});
writetable(tmp_long_OGTT,outfile) 



%% plot data
DATA=readtable('FIG4.xlsx');
fs=12;
fn='arial';
%% plot data

figure(1); clf


%%% load longitudinal simulation data %%%%
%load NGT_IGT_T2D.dat; 
NGT_IGT_T2D=readtable('FIG4.xlsx');



NGT_IGT_T2D_t=NGT_IGT_T2D{:,1};
NGT_IGT_T2D_FG=NGT_IGT_T2D{:,2};
NGT_IGT_T2D_2hG=NGT_IGT_T2D{:,4};
NGT_IGT_T2D_FI=NGT_IGT_T2D{:,5};
NGT_IGT_T2D_2hI=NGT_IGT_T2D{:,7};

NGT_IGT_T2D_si=NGT_IGT_T2D{:,12};
NGT_IGT_T2D_HGP=NGT_IGT_T2D{:,13};
NGT_IGT_T2D_sigma=NGT_IGT_T2D{:,11};
NGT_IGT_T2D_b=NGT_IGT_T2D{:,9};
%%% Note that t should be rescaled by 365*1440 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fs=12;
fs2=8;
fs10=8;
fn='arial';
lw=1.5;
lw2=1.5;
lw3=0.5;
fpan=10;


%%%%%%%%%%%%%%%%%%%%%

tsi=subplot(3,2,1);

plot(NGT_IGT_T2D_t./(365*1440),NGT_IGT_T2D_si.*6.9444,'k', 'linewidth',lw);
%hold('on')



%xlabel('time (year)','fontsize', fs, 'fontname',fn);
ylabel('{S_I} (10^{-4}ml/\muU/min)','fontsize', fs, 'fontname',fn);

%legend('compensation',5,'decompensation',5);
%legend('boxoff');

%% 10% of the length of y axis


text(0,6.6,'A','fontsize',fpan,'fontweight','bold');
axis ([0 5 0 6])
%%%%%%%%%%%%%%%%%%%%%%%%%%
tHGP=subplot(3,2,2);

plot(NGT_IGT_T2D_t./(365*1440),NGT_IGT_T2D_HGP,'k', 'linewidth',lw);
%hold('on')



%xlabel('time (year)','fontsize', fs, 'fontname',fn);

ylabel('{hepa_{S_I}}','fontsize', fs, 'fontname',fn);

%% 10% of the length of y axis

text(0,1.65,'B','fontsize',fpan,'fontweight','bold');
axis ([0 5 0 1.5])


%% 10%  of the length of y axis



%%%%%%%%%%%%%%%%%%%%%%%%
tG=subplot(3,2,3);


plot(NGT_IGT_T2D_t./(365*1440),NGT_IGT_T2D_FG,'k', 'linewidth',lw);

hold('on')
plot(NGT_IGT_T2D_t./(365*1440),NGT_IGT_T2D_2hG,'k', 'linewidth',lw);

tIGT=0.5;
tCGI=0.8;
tT2D=1.5;

th_IGT=140;
th_CGI=100;
th_T2D=200;

%plot([0 5],[140 140],'k','linewidth',lw);
%plot([0 5],[100 100],'b','linewidth',lw2);
%plot([0 5],[200 200],'r','linewidth',lw2);

%plot([tIGT tIGT],[0 th_IGT],'k','linewidth',lw);
%plot([tCGI tCGI],[0 th_CGI],'b','linewidth',lw2);
%plot([tT2D tT2D],[0 th_T2D],'r','linewidth',lw2);






%plot([0 5],[140 140],'k:','linewidth',lw);
%plot([0 5],[100 100],'k--','linewidth',lw2);
%plot([0 5],[200 200],'k-.','linewidth',lw2);

plot([tIGT tIGT],[0 th_IGT],'k','linewidth',lw3);
%plot([tCGI tCGI],[0 th_CGI],'k','linewidth',lw3);
plot([tT2D tT2D],[0 th_T2D],'k','linewidth',lw3);

%plot(0,0,'.k','MarkerSize',15);
%plot(1.1, 0, '.k', 'MarkerSize',15);
%plot(2.5, 0, '.k', 'MarkerSize',15);
%plot(3, 0, '.k', 'MarkerSize',15);


%hold('on')
%line([1.1 1.1],[0 140]);
%hold('on');
%line([5 0],[5 200]);


%xlabel('time (year)','fontsize', fs, 'fontname',fn);
ylabel('{Glucose} (mg/dl)','fontsize', fs, 'fontname',fn);


%% 10% of the length of y axis

text(0,440,'C','fontsize',fpan,'fontweight','bold');
text(3,325,'2hPG', 'fontsize',fs2,'fontname',fn);
text(3,165, 'FPG','fontsize',fs2,'fontname',fn);

%text(6.02,280,{'2-h', 'glucose'}, 'fontsize',fs2,'fontname',fn);
%text(6.02,130, {'fasting', 'glucose'},'fontsize',fs2,'fontname',fn);



text(0,40,'NGT','fontsize',fs10,'fontname',fn);
text(0.9,40,'IGT','fontsize',fs10,'fontname',fn);
%text(2.5,40,'CGI','fontsize',fs10,'fontname',fn);
text(3,40,'T2D','fontsize',fs10,'fontname',fn);

axis ([0 5 0 400])





%%% End of FIG4-C



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tI=subplot(3,2,4);
plot(NGT_IGT_T2D_t./(365*1440),NGT_IGT_T2D_FI,'k', 'linewidth',lw);

hold('on')
plot(NGT_IGT_T2D_t./(365*1440),NGT_IGT_T2D_2hI,'k', 'linewidth',lw);


%xlabel('time (year)','fontsize', fs, 'fontname',fn);
ylabel('{Insulin} ({\mu}U/ml)','fontsize', fs, 'fontname',fn);


text(0,440,'D','fontsize',fpan,'fontweight','bold');
text(2,200,'2hPI', 'fontsize',fs2,'fontname',fn);
text(2,48, 'FPI','fontsize',fs2,'fontname',fn);


axis ([0 5 0 400])



%%% End of FIG4-D

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tsigma=subplot(3,2,5);

plot(NGT_IGT_T2D_t./(365*1440),NGT_IGT_T2D_sigma,'k', 'linewidth',lw);
%hold('on')



xlabel('time (year)','fontsize', fs, 'fontname',fn);
ylabel('{\sigma}','fontsize', fs, 'fontname',fn);


%% 10% of the length of y axis

text(0,1.1,'E','fontsize',fpan,'fontweight','bold');
axis ([0 5 0 1])


%% 10%  of the length of y axis



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tb=subplot(3,2,6);

plot(NGT_IGT_T2D_t./(365*1440),NGT_IGT_T2D_b,'k', 'linewidth',lw);
%hold('on')



xlabel('time (year)','fontsize', fs, 'fontname',fn);
ylabel('{\beta} (mg)','fontsize', fs, 'fontname',fn);


%% 10% of the length of y axis

text(0,3300,'F','fontsize',fpan,'fontweight','bold');
axis ([0 5 0 3000])


%% 10%  of the length of y axis



