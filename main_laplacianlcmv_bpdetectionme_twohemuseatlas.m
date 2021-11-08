%%
clc
clear
%laplacian filter EEG
% load Laplacian_weight_64channel%Lap
% EEG_lsf=Lap*EEG;
load main_elec_new_aligned_64_add3signelectrode_mni
% x=elec_new1.chanpos(2:end-2,1);
% y=elec_new1.chanpos(2:end-2,2);
% z=elec_new1.chanpos(2:end-2,3);
% locs=[x y z];
% m=4;lambda=1e-5;
% [Et,Ep,L,S] = scalpef( locs , m , lambda );
load Et_left 
load Ep_left
load Et_right
load Ep_right
 
load main_elec_left_from64channel
for i=2:length(elec_left.label);numberleft(1,i-1)=find(strcmp(elec_new1.label,elec_left.label{1,i}))-1;end
load main_elec_right_from64channel
for i=2:length(elec_right.label);numberright(1,i-1)=find(strcmp(elec_new1.label,elec_right.label{1,i}))-1;end
load Laplacian_weight_64channel_left
load Laplacian_weight_64channel_right
 
cal=0
if cal==1
load simbio_vol_left_64_ftpreparevolsens
load simbio_volright_64_ftpreparevolsens
load simbio_volright_64_stiff
vol_right.transfer=vol2_right.transfer;
clear vol2_right;
 
if vol_right.transfer(1,:)==zeros(1,length(vol_right.transfer(1,:)));
vol_right.transfer=vol_right.transfer(2:end,:);
end
 
if vol_left.transfer(1,:)==zeros(1,length(vol_left.transfer(1,:)));
vol_left.transfer=vol_left.transfer(2:end,:);
end
 
for i=1:size(dip.pos,2);
dip.lfmright{1,i}=leadfield_simbio(dip.pos{1,i},vol_right);
end
 
for i=1:size(dip.pos,2);
dip.lfmleft{1,i}=leadfield_simbio(dip.pos{1,i},vol_left);
end
 
% save  dip_formeshmaskgreymatter_BrodMann4_usinggrayvol  dip
end%if cal
%dip.number is corresponded number of that point on mesh(1,1).pnt
load dip_formeshmaskgreymatter_BrodMann4_usinggrayvol;
 
% dip2=dip;
% load dip_leftright_meshmaskgrey_grayvol;
% dip.lfm=dip2.lfm;
% dip.pos=dip2.pos;
% clear dip2;
 
vol=[];
 
%%
% % % % % %[0.8 0.8 1] lightblue
%[0.5 0.5 0.5]is gray
% % cortex is[1 0.8 0.4]
% % % % % %brain is [0.79 0.39 0.39]
% % % % % %skin is [1 0.83 0.46]
%main_twosource_simulation_correlated_laplacian
%simulate one source and 2 correlated sources and reconstruct by lapalcian
%spatialfiltered EEG instead of EEG and lcmv beamformer
%%
%loading data
%electrode
 
% load main_elec_new_aligned_64_add3signelectrode_mni%electrode positions
sens=elec_new1;
load mesh_mask_near_brodmann4_mni_grey
    
%data
% load dip_formeshmaskgreymatter_BrodMann4_usinggrayvol;
load('mesh_3graylayers_mni.mat');
% load('main_elec_new_aligned_64_add3signelectrode_mni')
%%
%showing ROI on mesh
if 0
% NAI=0.3*zeros(size(mesh(1,1).pnt,1),3);%or zeros
% % % % % % % % load mesh_mask_near_brodmann1_4_6_8_mni_grey 
NAI2=[];
for i=1:size(mesh(1,1).pnt,1);NAI2(i,:)=[0.5 0.5 0.5];end%[1 0.8 0.4]is cortex
a=1;b=1000;%1000t0 10000 are c3 Fc3 C4 Fc4
%20000to 40000is for hand movement area of primary motor cortex
NAI2(mesh_mask_near.number(1,a:b),:)=[ones(b-a+1,1) zeros(b-a+1,1) zeros(b-a+1,1) ];
NAI2(mesh_mask_near.number(1,a+b:b+b),:)=[ones(b-a+1,1) zeros(b-a+1,1) zeros(b-a+1,1) ];
% plotting
figure;ft_plot_mesh(mesh(1,1),'vertexcolor',NAI2,'facealpha',0.4,'FaceColor','cortex','edgecolor','none','edgealpha',0.4)
% camlight;lighting gouraud
hold on;
select=0;
if select==1
ft_plot_mesh(mesh(1,1),'edgecolor','none','edgealpha',0.3,'FaceColor','cortex','facealpha',0.3);camlight;lighting gouraud
% always 'edgecolor','none' unless gives error and matlab stops
else
ft_plot_mesh(mesh(1,1),'surfaceonly','yes','vertexcolor','none','edgecolor','none','FaceColor',[0.5 0.5 0.5],'facealpha',0.7);
camlight;lighting gouraud
end
elec.pnt=elec_new1.chanpos;
elec.tri=convhulln(elec.pnt);
elec.label=elec_new1.label;
hold on;ft_plot_sens(elec,'style','k*');%'sr'
% hold on;ft_plot_mesh(elec.pnt,'facecolor','no')
% end%if 0
legend('primarymotor area-centralsulcus')
end
%%
load EEG_rest_S1_cleaned
% % % %EEG at *rest*
% % % % load('main_elec_new_aligned_64_add3signelectrode_mni')
% % % exist2=1;
% % % if exist2==0
% % % load('Rest_07_09_2016_12_50_16_0001.mat');
% % % fs=256;
% % % 
% % % % load data
% % % load 64_chan_mont;
% % % exist=1;
% % % if exist==0
% % % %load data
% % % InputSignal=Rest;
% % % data1=reshape(InputSignal,size(InputSignal,1),size(InputSignal,3));
% % % data=[data1 ];
% % % if 0
% % % save S1_L1_rest data
% % % end
% % % else
% % % load S1_L1_rest
% % % end
% % % % if 0
% % % % InputSignal1=Soo1R01;
% % % % data1=reshape(InputSignal1.Data,size(InputSignal1.Data,2),size(InputSignal1.Data,3));
% % % % end
% % % EEG_rest=data(2:end-1,:);
% % % figure;plot(EEG_rest')
% % % for i=1:size(Channel_Labels,2)
% % %     number2(1,i)=find(strcmp(elec_new1.label,Channel_Labels{1,i} ))-1;%because first one is left electrode
% % % EEG_rest(number2(1,i),:)=EEG1(i,:);
% % % end
% % % EEG_marker1=data(end,100*fs:end);
% % % 
% % % EEG_rest_S1_L1_eeglab=[EEG_rest;EEG_marker1];
% % % save EEG_rest_S1_L1_eeglab %with marker
% % % %after preprocessing with eeglab
% % % % EEG_rest=real(EEG.icawinv*EEG.icaact);
% % % save EEG_rest_S1_L1_cleaned  EEG_rest
% % % EEG_marker1=EEG_rest(end,:);
% % % 
% % % EEG=EEG_rest(1:end-1,:);
% % % save EEG_rest_S1_L1_cleaned  EEG_rest
% % % 
% % % 
% % % 
% % % % number=find(EEG_marker1~=0);
% % % % length(number)
% % % % counter=1;
% % % % markernumber=[];
% % % % for i=2:length(number);
% % % %     k=i-1;
% % % %     if number(1,i)-number(1,k)>100;%because the length of triger is 80 samples
% % % %         counter=counter+1;
% % % %         markernumber(1,counter)=number(1,i);
% % % % %         EEG_marker2(1,number(1,i))=1;
% % % %     else
% % % %                 markernumber(1,counter)=number(1,k)-100;%length of trigger is 100atmost
% % % % %         EEG_marker2(1,number(1,k))=1;
% % % % 
% % % %     end
% % % % end
% % % % if markernumber(1,1)==0;markernumber(1,1)=1;end
% % % % length(markernumber)
% % % % if 0
% % % % figure;plot(EEG_marker1)
% % % % end
% % % % EEG_marker=[];EEG_marker(1,markernumber)=1;
% % % % if 0
% % % % hold on;plot(EEG_marker,'r')
% % % % hold on;plot(EEG(10,:))
% % % % end
% % % % if 0
% % % % save EEG_S1_L1              EEG %used for preprocessing with eeglab
% % % % save EEG_marker_S1_L1       EEG_marker
% % % % save markernumber_S1_L1     markernumber
% % % % end
% % % save EEG_rest_S1_L1   EEG_rest
% % % else
% % %     eeglabprocessed=1;
% % %     if      eeglabprocessed==1;
% % % % EEG_rest=real(EEG.icawinv*EEG.icaact);
% % % % EEG_rest=real(ALLEEG(1, 2).icawinv*ALLEEG(1, 2).icaact );%or last line
% % % % save EEG_rest_S1_L1_cleaned  EEG_rest
% % % if 0
% % % %we can use *automatic channel rejection of eeglab* after baseline vandering
% % % %so we can use average of near electrodes to replace by bad channel signal
% % % %Fpz Fp2 F2 are bad channels
% % % 
% % % EEG(22,:)=(EEG(21,:)+EEG(26,:))/2;
% % % EEG(23,:)=(EEG(27,:)+EEG(28,:))/2;
% % % EEG(34,:)=(EEG(26,:)+EEG(27,:)+EEG(33,:)+EEG(35,:)+EEG(4,:))/5;
% % % 
% % % EEG_rest=EEG;
% % % save EEG_rest_S1_L1_cleaned   EEG_rest
% % % end
% % % load EEG_rest_S1_L1_cleaned
% % % 
% % %     else   
% % %  load   EEG_rest_S1_L1
% % %     end
% % % % load   EEG_marker_S1_L1
% % % % load   markernumber_S1_L1
% % % end
% % % if 0
% % % figure;plot(EEG_rest')
% % % counter=1;
% % % for i=+1:size(EEG,1)/2%size(EEG,1)
% % %     
% % % figure(1);hold on;subplot(4,8,counter);plot(EEG(i,:));title(i)
% % % xlabel(elec_new1.label{1,i+1})
% % % counter=counter+1;
% % % end
% % % %by watching carefully we see it's needed to preprocess EEG by eeglab
% % % %Af8 is not good electrode--> 28 elec_new1.label{1,28+1}
% % % %Fp2 Af4 F6 F8
% % % n1=find(strcmp(elec_new1.label,'Fp2'))-1;
% % % n2=find(strcmp(elec_new1.label,'Af4'))-1;
% % % n3=find(strcmp(elec_new1.label,'F6'))-1;
% % % n4=find(strcmp(elec_new1.label,'F6'))-1;
% % % EEG(28,:)=(EEG(n1,:)+EEG(n2,:)+EEG(n3,:)+EEG(n4,:))/4;
% % % save EEG_S1_L1_cleaned EEG;
% % % end
%%
% load EEG_S1_L1_cleaned
load EEG_S1_L1
load EEG_marker_S1_L1
load markernumber_S1_L1
 
if 0
%ffttest
%%%%%%%%%%%%%%%
data=detrend(EEG(10,:));
timeexist=0;
if timeexist==1
t=InputSignal.Time;
Fs = round(length(t)/t(end))% Sampling frequency
else
Fs=256;
end
T = 1/Fs;                     % Sample time
L = length(data);                     % Length of signal
 
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(data,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);
 
% Plot single-sided amplitude spectrum.
figure;plot(f,2*abs(Y(1:NFFT/2+1))) 
title('Single-Sided Amplitude Spectrum of recorded signal')
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')
 
% print(gcf,'-djpeg','-r300','Single-Sided Amplitude Spectrum of recorded signal');
end
%%
%bandpass filtering all EEG and see template
fs=256;
w=[0.1/fs 10/fs]%[0.05/fs 10/fs];%[16/fs 22/fs]is beta [11 14]is mu
[a,b]=butter(3,w(1,1),'high');
[a2,b2]=butter(3,w(1,2),'low');
signal1=[]; signal2=[];
m=size(EEG,1)
signal2_2=[];
for i=1:m
signal1(i,:)=filtfilt(a,b,EEG(i,:));%highpass
signal2(i,:)=filtfilt(a2,b2,signal1(i,:));%low pass
% signal2(i,:)=downsample(signal2(i,:),2,2);
signal2_2(i,:)=filtfilt(a2,b2,EEG(i,:));%low pass
end
 
EEG3=signal2;
% % % for i=1:m
% % % EEG3(i,:)=decimate(EEG(i,:)',2,10);
% % % end
if 0
figure;plot(signal2_2(:,:)');hold on;plot(EEG_marker(1,:),'g')
 
title('FC4  bandpass filtered')
end
%align electrodes on sphere for calculating spherical spline laplacian
%%%%%%%%%%find template of BP for right hand movement%%%%%%%%%%%%%%%%%%%%%%%%
 
markernumber;BP_templateright=[];
l_marker=size(markernumber,2);
BP_temp=[];
BP=[];
for num_chan=1:size(signal2,1)
BP=zeros(1,3/2*fs+1+0*fs);
for  i=1:l_marker
    if markernumber(1,i)==1;i=1+i;end
BP_temp=EEG3(num_chan,markernumber(1,i)-3/2*fs-1:markernumber(1,i)+0*fs-1)+BP;
BP=BP_temp;
end
BP_templateright(num_chan,:)=BP_temp/size(markernumber,2);
%averaging over length(marker_number)= 29 epochs 
end
% 8 10 and 12 are C3 Cz C4 %32 is F3 %36 is F4
% 2 is FC3 %Fz34
%for laplacian left 4 1 14  are c3  Fc3  F3 
%32 F3 36 F4
 
figure;plot(-3/2-1/fs:(3/2+1/fs)/((size(BP_templateright(1,:),2))-1):0,BP_templateright');ylabel('micro volt')%8 c3 10 Cz 12 C4
 
xlabel('time(s)');ylabel(' voltage(uV) ')
title('C4 template for left hand achieved by synchronized averaging of bandpass filtered-EEG fs=256Hz ')
title('RP template before movement onset')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
%decreasing the data samples
if 0
p=randperm(size(BP_templateright,2),100);
p2=sort(p)
BP_templateright2=BP_templateright(:,p2);
 
figure;plot(-2:4/(length(BP_templateright2(1,:))-1):2,BP_templateright2(1,:));ylabel('micro volt')%8 c3 10 Cz 12 C4
title('reduced samples by randperm')
end
%%
if 0
%process data by lcmv and narrow band filter
%inorder to achieve ERDERS curve and
%NAI over the ROI and also signals of positions of maximum NAI
%for premotor and SMA(BM6) and primary motor and centralsulcus(BM4)
%data
% if 0
% EEG;
% EEG3;
% EEG_marker;
% markernumber;
% end
% load dip_formeshmaskgreymatter_BrodMann4_usinggrayvol;
% alldip=dip;
%dip.pos & dip.lfm for primarymotorcortex(2000points of 80000points) 
% dip.number is corresponding number of points on mesh
% load mesh_3graylayers_mni
vol=[];%not needed
if 0
load simbio_vol_stiff_3graylayers_mni
load simbio_vol_ft_prepare_vol_sens_64add3refchan_3greylayers_mni%transfer for grey matter 
load simbio_vol_stiff_3graylayers_mni;
vol.transfer=vol2.transfer;
clear vol2;
end
 
 
%
start=1;%start of process
stop=size(EEG,2);
l=round((stop-start)/(fs/2));
 
train=0;
%if train data this program finds location and signal of active sources
%by searching the whole ROI and finds best locations 1(s) before movement
length=fs/2;%length of EEG in each trial it should be stationary for beamforming
countertrial=1;
 
if train==1
for t=start:fs/2:l*fs/2%time 
EEG4=EEG3(:,t:t+length);
C=EEG4*EEG4';
[eigv,eigval]=eig(C,'nobalance');
[a,b]=find(eigval==max(max(eigval)));
landa=max(max(eigval));
eigvector=eigv(:,a);
[U,S,V]=svd(C);
C=C+0.03*landa*U;
 
for locnumber=1:size(dip.pos,2)-80
lfm1=dip.lfm{1,locnumber};
if size(lfm1,1)>64;lfm1=lfm1(2:end-2,:);end
for k=1:3;lfm1(:,k)=lfm1(:,k)/norm(lfm1(:,k));end
    pos=dip.pos{1,locnumber};
for n=1:80;node(n,:)=dip.pos{1,locnumber+n-1};end
X=node(:,1);Y=node(:,2);Z=node(:,3);
    dt = delaunayTriangulation(X,Y,Z);
    [Tfb,Xfb] = freeBoundary(dt);
    TR = triangulation(Tfb,Xfb);
%     vn = faceNormal(TR);
 vn = vertexNormal(TR);
 dip.momanatomyright(locnumber,:)=vn(1,:);
 momanatomy=vn(1,:)';
optmom  = optimummomentum_lcmv_sekihara2015( C,lfm1,EEG4,pos,vol,momanatomy );
%   optmom.mom and optmom.NAI
dip.optmom{countertrial,locnumber}=optmom.mom;
dip.NAI(countertrial,locnumber)=optmom.NAI ; 
 
end
aa=find(dip.NAI(countertrial,:)==max(dip.NAI(countertrial,:)));
dip.numbermaxNAI(countertrial,1)=aa(1,1);
dip.numbermaxNAI2(countertrial,1)=dip.number{1,dip.numbermaxNAI(countertrial,1)};
dip.mommaxNAI(countertrial,:)= dip.optmom{countertrial,dip.numbermaxNAI(countertrial,1)}';
 
countertrial=countertrial+1;
end
if 0
save dip_S1_L1  dip
end
dip2=dip;
load dip_formeshmaskgreymatter_BrodMann4_usinggrayvol;
% dip.pos=dip2.pos;
% dip.number=dip2.number;
% dip.lfm=dip2.lfm;
dip.numbermaxNAI=dip2.numbermaxNAI;
dip.numbermaxNAI2=dip2.numbermaxNAI2;
dip.mommaxNAI=dip2.mommaxNAI;
dip.NAI=dip2.NAI;
 
save dip_S1_L2 dip
load  dip_S1_L2   %number of probable positions
a=1;
counter=1;
for j=1:size(markernumber,2);%start+(countertrial-1)*fs/2;%start
countertrial_marker(1,j)= round((markernumber(1,j)-start)*2/fs)+1;   
 
maxNAInumber(1,counter)=dip.numbermaxNAI(countertrial_marker(1,j)-a,1);
maxNAInumber2(1,counter)=dip.numbermaxNAI2(countertrial_marker(1,j)-a,1);
mom_maxNAInumber2(counter,:)=dip.mommaxNAI(countertrial_marker(1,j)-a,:);
NAItrial(counter,:)=dip.NAI(countertrial_marker(1,j)-a,:);
 
maxNAInumber(1,counter+1)=dip.numbermaxNAI(countertrial_marker(1,j)-a,1);
maxNAInumber2(1,counter+1)=dip.numbermaxNAI2(countertrial_marker(1,j)-a,1);
mom_maxNAInumber2(counter+1,:)=dip.mommaxNAI(countertrial_marker(1,j)-a,:);
NAItrial(counter+1,:)=dip.NAI(countertrial_marker(1,j)-a,:);
 
maxNAInumber(1,counter+2)=dip.numbermaxNAI(countertrial_marker(1,j)-a,1);
maxNAInumber2(1,counter+2)=dip.numbermaxNAI2(countertrial_marker(1,j)-a,1);
mom_maxNAInumber2(counter+2,:)=dip.mommaxNAI(countertrial_marker(1,j)-a,:);
NAItrial(counter+2,:)=dip.NAI(countertrial_marker(1,j)-a,:);
 
counter=counter+3;
end
save maxNAInumber_S1_L2   maxNAInumber
save maxNAInumber2_S1_L2  maxNAInumber2%for positioning between whole dipole positions
save mom_maxNAInumber2_S1_L2  mom_maxNAInumber2
save NAItrial_S1_L2  NAItrial
 
else
     load  dip_S1_L2   %dip.pos number lfm lfm_right lfm_left
 
    load maxNAInumber_S1_L2
    load maxNAInumber2_S1_L2
    load mom_maxNAInumber2_S1_L2
    load NAItrial_S1_L2
end
%finding anatomical location and momentum for comparison
exist=1;
if exist==0
momanatomy=[];
counter=1;
for i=1:size(maxNAInumber,2)%location of maximum NAI
 pos(i,:)=dip.pos{1,maxNAInumber(1,i)};
 for n=1:80;node(n,:)=dip.pos{1,maxNAInumber(1,i)+n-1};end
X=node(:,1);Y=node(:,2);Z=node(:,3);
    dt = delaunayTriangulation(X,Y,Z);
    [Tfb,Xfb] = freeBoundary(dt);
    TR = triangulation(Tfb,Xfb);
%     vn = faceNormal(TR);
 vn = vertexNormal(TR);
 momanatomy(i,:)=vn(1,:);
% %  momanatomy=vn(1,:)';
end
save momanatomy_S1_L2  momanatomy
save pos_S1_L2   pos
else
load momanatomy_S1_L2
load pos_S1_L2
end
 
%showing momentum on mesh
show=0;
if show==1
load mesh2;%with less nodes for better visualization
figure(1);
ft_plot_mesh(mesh2(1,1),'edgecolor','none','edgealpha',0.5,'FaceColor','cortex','facealpha',0.8);camlight;lighting gouraud
counter=1;
for i=1:size(pos,1)-3
figure(1);
hold on;quiver3(pos(i+3,1),pos(i+3,2),pos(i+3,3),mom_maxNAInumber2(i+3,1), mom_maxNAInumber2(i+3,2),mom_maxNAInumber2(i+3,3) ,8,'r');
hold on;quiver3(pos(i+3,1),pos(i+3,2),pos(i+3,3),momanatomy(i+3,1), momanatomy(i+3,2),momanatomy(i+3,3) ,8,'g');
end
title('optimum direction of sources green is direction byanatomy')
end
end%if 0
%%
load L_right
load L_left
load Laplacian_weight_64channel
%%
%choosing predefined locations of primarymotor cortex
poses=[];dip.lfm2=[];
a3=454;b3=460;%some points of near right hand and with NS' component
counter=1;
for i=a3:b3;
    poses(counter,:)=dip.pos{1,i};
    dip.lfmright2{1,counter}=dip.lfmright{1,i};
maxNAInumber(1,counter)=i;
    counter=counter+1;
end
% poses(counter,:)=[2.5,-3.5,67.5];%for SMA
% dip.lfmright2{1,counter}=dip.lfmright{1,i};
% maxNAInumber(1,counter)=i;
% counter=counter+1;
for i=a3+1010:b3+1010;poses(counter,:)=dip.pos{1,i};
    dip.lfmleft2{1,counter}=dip.lfmleft{1,i};
    maxNAInumber(1,counter)=i;
    counter=counter+1;
end
pos=poses;
tic
%see sources of predefined locations of maximum activity for test data
% and also ERDERS test
% load EEG at rest for ERDERS
mainsourcetrial=[];
mainERDERStrial=[];
posmainsource=[];
 
%make preprocessed and bandpass filtered (EEG3)
fs=256;
length=round(3*fs*1/2);
length1=fs/2;
start=1;%start of process
stop=size(EEG,2);
l2=round((stop-start)/(fs/2))
countertrial=1;
for t=start:fs/2:l2*fs/2%time 
    if 0
EEG4=EEG3(:,t:t+length);
C=EEG4*EEG4';
[eigv,eigval]=eig(C,'nobalance');
[a,b]=find(eigval==max(max(eigval)));
landa=max(max(eigval));
eigvector=eigv(:,a);
[U,S,V]=svd(C);
C=C+0.03*landa*U;
    end
alleegtrial{1,countertrial}=Lap*EEG3(:,t:t+length);

    %%%%%%%%%%%%%for ERDERS
%first bandpass filter each epoch
 
w=[16/fs 22/fs]%[0.05/fs 10/fs];%[16/fs 30/fs]is beta [8 14]is mu
[a,b]=butter(3,w(1,1),'high');
[a2,b2]=butter(3,w(1,2),'low');
signal1=[]; signal2=[];
m=size(EEG,1)
for i2=1:m
signal1(i2,t:t+length1)=filtfilt(a,b,EEG(i2,t:t+length1));%highpass %startsub:t is 2seconds
signal2(i2,t:t+length1)=filtfilt(a2,b2,signal1(i2,t:t+length1));%low pass
signal2_2(i2,t:t+length1)=filtfilt(a2,b2,signal1(i2,t:t+length1));%low pass
end
% Ct=signal2_2(:,t:t+length)*signal2_2(:,t:t+length)';
 
w=[16/fs 22/fs];%[0.05/fs 10/fs];%[16/fs 22/fs]is beta [11 14]is mu
[a,b]=butter(3,w(1,1),'high');
[a2,b2]=butter(3,w(1,2),'low');
signal1=[]; signal2=[];
m=size(EEG,1)
for i2=1:m
signal1(i2,1:2*fs)=filtfilt(a,b,EEG_rest(i2,1:2*fs));%,markernumber(1,1)-6*fs:markernumber(1,1)-4*fs-1));%highpass
%calculate at rest (1:2*fs)
signal2(i2,1:2*fs)=filtfilt(a2,b2,signal1(i2,1:2*fs));%low pass
signal2_3(i2,1:2*fs)=filtfilt(a2,b2,signal1(i2,1:2*fs));%low pass
end
if 0
Ct=signal2_2(:,t:t+length)*signal2_2(:,t:t+length)';
Cbg=signal2_3(:,1:2*fs)*signal2_3(:,1:2*fs)';
Ctot=Ct+Cbg;
[U1,S1,V1]=svd(Ctot);
[eigv,eigval]=eig(Ctot);
[a,b]=find(eigval==max(max(eigval)));
landa=max(max(eigval));
Ctot=Ctot+0.03*landa*U1;
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
counter1=1;
for i=1:size(poses,1)%-3   size(pos,1)
    pos1(counter1,:)=poses(i,:);
    %lfm
C=[];
Laplac_right=Lap_right;%Et_right+Ep_right;%Lap_right
    if poses(i,1)>0;lfm1=Laplac_right*dip.lfmright2{1,i};%dip.lfmright{1,maxNAInumber(1,i)}
% Et_right %Lap_right %Ep_right %L_right
        C=(Laplac_right*EEG3(numberright,t:t+length))*(Laplac_right*EEG3(numberright,t:t+length))';
EEG4=Laplac_right*EEG3(numberright,t:t+length);
    [eigv,eigval]=eig(C,'nobalance');
[a,b]=find(eigval==max(max(eigval)));
landa=max(max(eigval));
eigvector=eigv(:,a);
[U,S,V]=svd(C);
C=C+0.03*landa*U; 
Ct=(Laplac_right*signal2_2(numberright,t:t+length1))*(Laplac_right*signal2_2(numberright,t:t+length1))';
Cbg=(Laplac_right*signal2_3(numberright,1:2*fs))*(Laplac_right*signal2_3(numberright,1:2*fs))';
Ctot=Ct+Cbg;
[U1,S1,V1]=svd(Ctot);
[eigv,eigval]=eig(Ctot);
[a,b]=find(eigval==max(max(eigval)));
landa=max(max(eigval));
Ctot=Ctot+0.03*landa*U1;
    end
    %maxNAInumber2 equivalent to dip.number from 80000nodes 
Laplac_left=Lap_left;%Et_left+Ep_left;%Lap_left
    if poses(i,1)<0;lfm1=Laplac_left*dip.lfmleft2{1,i};%dip.lfmleft{1,maxNAInumber(1,i)}
% Et_left %Lap_left or L_left %Ep_left
 
        C=(Laplac_left*EEG3(numberleft,t:t+length))*(Laplac_left*EEG3(numberleft,t:t+length))';
EEG4=Laplac_left*EEG3(numberleft,t:t+length);
[eigv,eigval]=eig(C,'nobalance');
[a,b]=find(eigval==max(max(eigval)));
landa=max(max(eigval));
eigvector=eigv(:,a);
[U,S,V]=svd(C);
C=C+0.03*landa*U; 
%%%%%%%%%%%%%%%%%%%%%%    
Ct=(Laplac_left*signal2_2(numberleft,t:t+length1))*(Laplac_left*signal2_2(numberleft,t:t+length1))';
Cbg=(Laplac_left*signal2_3(numberleft,1:2*fs))*(Laplac_left*signal2_3(numberleft,1:2*fs))';
Ctot=Ct+Cbg;
[U1,S1,V1]=svd(Ctot);
[eigv,eigval]=eig(Ctot);
[a,b]=find(eigval==max(max(eigval)));
landa=max(max(eigval));
Ctot=Ctot+0.03*landa*U1;
   
    end
 
%     lfm1=dip.lfm{1,maxNAInumber(1,i)};
%maxNAInumber is number of activenode for 2000 locations
if size(lfm1,1)>64;lfm1=lfm1(2:end-2,:);end
for k=1:3;lfm1(:,k)=lfm1(:,k)/norm(lfm1(:,k));end
 
    for n=1:80;node(n,:)=dip.pos{1,maxNAInumber(1,i)+n-1};end
X=node(:,1);Y=node(:,2);Z=node(:,3);
    dt = delaunayTriangulation(X,Y,Z);
    [Tfb,Xfb] = freeBoundary(dt);
    TR = triangulation(Tfb,Xfb);
%     vn = faceNormal(TR);
 vn = vertexNormal(TR);
 momanatomy1=vn(1,:)';
    % optimum momentum
optmom  = optimummomentum_lcmv_sekihara2015( C,lfm1,EEG4,pos1,vol,momanatomy1 );
%optmom.mom;optmom.l;optmom.W; optmom.NAI; optmom.Wmomanatomy; optmom.NAIlcmvmomanatomy
sources(counter1,t:t+length)=(optmom.W'*EEG4);
NAI(1,counter1)=optmom.NAI;
l=lfm1*optmom.mom;
W=(inv(Ctot)*l)/(l'*inv(Ctot)*l);
moms(counter1,:)=optmom.mom';
F_factor=(W'*Ctot*W - W'*Cbg*W)/(W'*Cbg*W);
totalERDERS(counter1,t:t+length1)=F_factor;%t-fs/2:t
counter1=counter1+1;%for position
end%predefined locations
aa=find(NAI==max(NAI));
numbermainsource(1,countertrial)=aa(1,1);
mainsource(1,t:t+length)=sources(numbermainsource(1,countertrial),t:t+length);%without attention to position
mainsourcetrial{1,countertrial}=mainsource(1,t:t+length);%the same as last line but incell
allsourcetrial{1,countertrial}=sources(:,t:t+length);
allpowerbandtrial{1,countertrial}=powerband(sources(:,t:t+length),fs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mainERDERS(1,t:t+length1)=totalERDERS(numbermainsource(1,countertrial),t:t+length1);
mainERDERStrial{1,countertrial}=mainERDERS(1,t:t+length1);
allERDERS(:,t:t+length1)=totalERDERS(:,t:t+length1);%for all positions
posmainsource(countertrial,:)=pos1(numbermainsource(1,countertrial),:);
mommainsource(countertrial,:)=moms(numbermainsource(1,countertrial),:)
allNAItrial{1,countertrial}=NAI';
countertrial=countertrial+1;
end%time instant
toc
%%
save alleegtrial_S1_lap       alleegtrial 

save allsourcetrial_S1_L1_lap allsourcetrial
save posmainsource_s1_L1_lap  posmainsource
save mommainsource_s1_L1_lap   mommainsource
save allNAItrial_S1_L1_lap  allNAItrial
save allERDERS_S1_L1_lap allERDERS
save allpowerbandtrial_S1_L1_lap  allpowerbandtrial
save pos1_S1_L1_lap pos1
%%
if 0
%showing signals of active points of all trials before movement
size(mainsourcetrial)
size (markernumber)
for n=1:size (markernumber,2)
a=3
countertrial1=round((markernumber(1,n)-1)*2/fs)+1;
%most powerful sourceposition and signal
figure(11);hold on;subplot(4,7,n)
plot(mainsourcetrial{1,countertrial1-a})%a=-3 is 1.5s beforemovement onset
% /norm(mainsourcetrial{1,countertrial1-a})
hold on;
plot(mainsourcetrial{1,countertrial1-a-1},'g')%a=-3 is 1.5s beforemovement onset
% /norm(mainsourcetrial{1,countertrial1-a-1})
title(posmainsource(countertrial1-a,:))
 
end
 
countertrial1=round((markernumber(1,n)-1)*2/fs)+1
posmainsource(countertrial1,:)%most powerful sourceposition
%see active sources in consequent epochs of a trial of movement before
%movement
color={'r' 'g' 'b' 'k'}
counter=1
for i=countertrial1-4:countertrial1
%     start+i*fs/2:start+i*fs/2+length
%/norm(mainsourcetrial{1,i})
figure(1);hold on;plot( start+i*fs/2:start+i*fs/2+length,mainsourcetrial{1,i}/norm(mainsourcetrial{1,i}),color{1,counter})
counter=counter+1;
end
end
%%
%see signals of all predefinedlocations for one trial of movement
n=14;size (markernumber,2)
a=3
countertrial1=round((markernumber(1,n)-start)*2/fs)+1;
posmainsource(countertrial1-a,:)%most powerful sourceposition
counter2=1;
for counter1=1:size(pos1,1)
    figure(5);hold on;subplot(4,9,counter2);
    plot(allsourcetrial{1,countertrial1-a+1}(counter1,:),'k');
%     plot(allsourcetrial{1,countertrial1-a-1}(counter1,:),'m');
%  hold on;
%    plot(allsourcetrial{1,countertrial1-a-2}(counter1,:)/norm(allsourcetrial{1,countertrial1-a-2}(counter1,:)),'g');
counter2=counter2+1;
title(pos1(counter1,:))
end
% title('trial 1 of movement -3 consequent epochs g m r start from 2.5 2 and 1.5s before onset duration=1.5s')
%%
% all sources of all predefinedpositions together for trials of movement a=3 meams ->start1.5 s before movementonset  red is sum of signals
% n=2%trial number
% countertrial1=round((markernumber(1,n)-1)*2/fs)+1;
a=3
sum2=zeros(size(allsourcetrial{1,countertrial1-a},1),size(allsourcetrial{1,countertrial1-a},2));
for n=1:size (markernumber,2)
a=3
countertrial1=round((markernumber(1,n)-1)*2/fs)+1;
 
figure(8);hold on;subplot(4,7,n)
plot(allsourcetrial{1,countertrial1-a}');
sum1=sum(allsourcetrial{1,countertrial1-a});%sum over locations
sum2=sum2+allsourcetrial{1,countertrial1-a};%sum of signals of each location over trials 
%  hold on;
% plot(sum1,'r')
%     
title(posmainsource(countertrial1-a,:))%(posmainsource(countertrial1-a,:))
end
%plot sum of time series
for i=1:size(pos,1)
figure(9);hold on;subplot(4,7,i)
plot(sum2(i,:)/28);title(pos(i,:))
end
% template of source space for predefined locations start=1.5s befre movement onset using laplacian=Et+Ep
%template of source space
 
%a  NS %1.5s before movement
%a-1 and a-2 are BP and maybe NS in tail
%a-3 above rest
%choose good trials for BP and NS'-->1 and rest 
 
%%
%saving all good sources or sources of some positions for classifying
load  alleegtrial_S1_lap     

load allsourcetrial_S1_L1_lap
load markernumber_S1_L1
load allpowerbandtrial_S1_L1_lap
good_NS=[(1:size(markernumber,2))];%a
% good_BP1=[1 4 5 6 7 9 14 117 18 21 22 24 26 28];%a-2
% good_BP2=[4 5 7 8 10 17 20 22 24 25 27 28 ];%a-1
good_rest=[(1:size (markernumber,2))];%a-5
 
size (markernumber)
counter1=1;
for n1=1:size(good_NS,2)%1:size (markernumber,2)
a=3
n=good_NS(1,n1);
countertrial1=round((markernumber(1,n)-1)*2/fs)+1;
%most powerful sourceposition and signal
class_NS{1,counter1}=(allsourcetrial{1,countertrial1-a});%a=-3 is 1.5s beforemovement onset
class_NS{1,counter1+1}=(allsourcetrial{1,countertrial1-a-1});%a=-3 is 1.5s beforemovement onset
class_eegNS{1,counter1}=alleegtrial{1,countertrial1-a};
class_eegNS{1,counter1+1}=alleegtrial{1,countertrial1-a-1};

% allpowerbandtrial
class_NS_power{1,counter1}=(allpowerbandtrial{1,countertrial1-a});
class_NS_power{1,counter1+1}=(allpowerbandtrial{1,countertrial1-a-1});
 
class_rest{1,counter1}=(allsourcetrial{1,countertrial1-a-3})%a=-3 is 1.5s beforemovement onset
class_rest{1,counter1+1}=(allsourcetrial{1,countertrial1-a-4})%a=-3 is 1.5s beforemovement onset
class_eegrest{1,counter1}=alleegtrial{1,countertrial1-a-3};
class_eegrest{1,counter1+1}=alleegtrial{1,countertrial1-a-4};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class_rest_power{1,counter1}=(allpowerbandtrial{1,countertrial1-a-3});
class_rest_power{1,counter1+1}=(allpowerbandtrial{1,countertrial1-a-4});
 
counter1=counter1+2;
% title(posmainsource(countertrial1-a,:))
end
%%
save class_eegrest_S1_L1_lapme   class_eegrest
save class_eegNS_S1_L1_lapme    class_eegNS

save class_NS_S1_L1_lap  class_NS;
% % save class_BP2_S1_L1  class_BP2
save class_rest_S1_L1_lap  class_rest
save class_rest_power_S1_L1_lap  class_rest_power
save class_NS_power_S1_L1_lap  class_NS_power
 
% figure;plot(allsourcetrial{1,countertrial1-a}')
%%
if 0
if 0
    figure;plot(mainERDERStrial{1,countertrial1-3})
end
% allERDERS(:,t:t+length1)=totalERDERS(:,t:t+length1);%for all positions
for n=1:size (markernumber,2)
a=3
countertrial1=round((markernumber(1,n)-1)*2/fs)+1;
 
figure(11);hold on;subplot(4,7,n)
plot(allERDERS(:,markernumber(1,n)-4*fs:markernumber(1,n)+fs)')
sum1=sum(allERDERS(:,markernumber(1,n)-4*fs:markernumber(1,n)+fs));
hold on;
plot(sum1,'r')
%     
title(posmainsource(countertrial1-a,:))%(posmainsource(countertrial1-a,:))
end
 
 
 
% % figure;plot(1e-6*mainsource);
% % % hold on;plot(EEG_marker,'g');%markernumber
if 0
 figure;plot((mainERDERS));
 hold on;plot(20*EEG_marker,'g');%markernumber
end
% figure;plot(1e-5*sources(1,:));hold on;plot(EEG_marker,'g')
if 0
n=16
countertrial1=round((markernumber(1,n)-1)*2/fs)+1;
numbermainsource(1,countertrial1)
s1=markernumber(1,n)-10*fs;
s2=markernumber(1,n)+10*fs;
figure;plot(totalERDERS(numbermainsource(1,countertrial1),s1:s2),'b')
hold on;plot(10*EEG_marker(1,s1:s2),'g');
end
 
 
 
% counter=1
% for counter1=1:size(pos1,1)
%     figure(1);hold on;subplot(4,7,counter);
%     plot(totalERDERS(counter1,:))
%     counter=counter+1;
% end
end
%%
%first showing predefinedpositions of sources from traindata
load mesh2
figure(3);ft_plot_mesh(mesh2(1,1),'surfaceonly','yes','vertexcolor','none','edgecolor','none','FaceColor',[0.5 0.5 0.5],'facealpha',0.8,'edgealpha',0.7);camlight;lighting gouraud
for n=1:size(pos,1)
figure(3);hold on;
quiver3(pos(n,1),pos(n,2),pos(n,3),0,0,0,'rp')
end
 
%%
%showing predefined points of 2000 points of primarymotor cortex for analysis of handmovement
dip.pos
poses=[];dip.lfm2=[];
a3=440;b3=455;
counter=1;
for i=a3:b3;poses(counter,:)=dip.pos{1,i};dip.lfm2{1,counter}=dip.lfm{1,i};counter=counter+1;end
poses(counter,:)=[2.5,-3.5,67.5];%for SMA
dip.lfm2{1,counter}=dip.lfm{1,i};
counter=counter+1;
for i=a3+1010:b3+1010;poses(counter,:)=dip.pos{1,i};dip.lfm2{1,counter}=dip.lfm{1,i};counter=counter+1;end
load mesh2
figure(5);ft_plot_mesh(mesh2(1,1),'surfaceonly','yes','vertexcolor','none','edgecolor','none','FaceColor',[0.5 0.5 0.5],'facealpha',0.8,'edgealpha',0.7);camlight;lighting gouraud
for n=1:size(poses,1)
figure(5);hold on;
quiver3(poses(n,1),poses(n,2),poses(n,3)+7,0,0,0,'rp')
end
%%
%showing allpositions of maximum activity before movement on mesh for all
%trials of movement
load mesh2;
a=3
size (markernumber,2)
figure(4);ft_plot_mesh(mesh2(1,1),'surfaceonly','yes','vertexcolor','none','edgecolor','none','FaceColor',[0.5 0.5 0.5],'facealpha',0.8,'edgealpha',0.7);camlight;lighting gouraud
for n=1:size (markernumber,2)
figure(4);hold on;
countertrial1=round((markernumber(1,n)-1)*2/fs)+1
posmainsource(countertrial1-a,:)%3 is 1.5 s befoem movementonset %most powerful sourceposition
quiver3(posmainsource(countertrial1-a,1),posmainsource(countertrial1-a,2),posmainsource(countertrial1-a,3),mommainsource(countertrial1-a,1),mommainsource(countertrial1-a,2),mommainsource(countertrial1-a,3),10,'rp')%+10for showing better
end
title('source positions for trials of each 1.5s before movement each position=one of 28 trials so some positions are aligned')
%%
%showing 3 consequent epochs and their corresponding source position
% posmainsource(countertrial,:)
for i=1:size (markernumber,2)
figure(8);subplot(4,7,i)
ft_plot_mesh(mesh2(1,1),'surfaceonly','yes','vertexcolor','none','edgecolor','none','FaceColor',[0.5 0.5 0.5],'facealpha',0.8,'edgealpha',0.7);camlight;lighting gouraud
end 
for n=1:size (markernumber,2)
subplot(4,7,n);hold on;
countertrial1=round((markernumber(1,n)-1)*2/fs)+1
%posmainsource(countertrial1-a,:)%3 is 1.5 s befoem movementonset %most powerful sourceposition
positions(1,:)=posmainsource(countertrial1-a,:);
positions(2,:)=posmainsource(countertrial1-a-1,:);
positions(3,:)=posmainsource(countertrial1-a-2,:);
momentums(1,:)=mommainsource(countertrial1-a,:);
momentums(2,:)=mommainsource(countertrial1-a-1,:);
momentums(3,:)=mommainsource(countertrial1-a-2,:);
 
quiver3(positions(:,1),positions(:,2),positions(:,3),momentums(:,1),momentums(:,2),momentums(:,3),'rp')%+10for showing better
end
%%
%prepare data for classify
exist=0
if exist==0
load class_rest_power_S1_L1_lap
load class_NS_power_S1_L1_lap
load class_NS_S1_L1_lap
load class_rest_S1_L1_lap
load class_eegrest_S1_L1_lap
load class_eegNS_S1_L1_lap
data3=[]
for  i=1:size(class_NS,2)
   counter=1;
    for j=1:size(class_NS{1,1},2):size(class_NS{1,1},1)*size(class_NS{1,1},2)
   data3(1,i)=1;
   data3(j+1:j+1+size(class_NS{1,1},2)-1,i)=class_NS{1,i}(counter,:);
   dataeeg3(1,i)=1;
   dataeeg3(j+1:j+1+size(class_eegNS{1,1},2)-1,i)=class_eegNS{1,i}(counter,:);
   
   counter=counter+1;
   end
end
data4=[]
for  i=1:size(class_NS_power,2)
   counter=1;
    for j=1:size(class_NS_power{1,1},2):size(class_NS_power{1,1},1)*size(class_NS_power{1,1},2)
   data4(1,i)=1;
   data4(j+1:j+1+size(class_NS_power{1,1},2)-1,i)=class_NS_power{1,i}(counter,:);
   counter=counter+1;
   end
end
% data3=([class_NS{1,1};class_NS{1,2};class_NS{1,3};class_NS{1,4};class_NS{1,5};class_NS{1,6};class_NS{1,7};class_NS{1,8};class_NS{1,9};...
% class_NS{1,10};class_NS{1,11};class_NS{1,12};class_NS{1,13};...
% class_NS{1,14};class_NS{1,15};class_NS{1,16};class_NS{1,17};...
% class_NS{1,18};class_NS{1,19};class_NS{1,20};class_NS{1,21};...
% class_NS{1,22};class_NS{1,23};class_NS{1,24};class_NS{1,25};class_NS{1,26};...
% class_NS{1,27};class_NS{1,28};...
% class_NS{1,29};class_NS{1,30};class_NS{1,31};class_NS{1,32};...
% class_NS{1,33};class_NS{1,34};class_NS{1,35};class_NS{1,36};class_NS{1,37};...
% class_NS{1,38};class_NS{1,39};class_NS{1,40};class_NS{1,41};class_NS{1,42};...
% class_NS{1,43};class_NS{1,44};class_NS{1,45};class_NS{1,46};class_NS{1,47};...
% class_NS{1,48};class_NS{1,49};class_NS{1,50};class_NS{1,51};class_NS{1,52};...
% class_NS{1,53};class_NS{1,54};class_NS{1,55};class_NS{1,56}]);
% 
% data3=data3';
% 
% for i=1:size(data3,2)
% data6(1,i)=1;
% data6(2:size(data3,1)+1,i)=data3(:,i);    
% end
data2=[]
for  i=1:size(class_rest,2)
   counter=1;
    for j=1:size(class_rest{1,1},2):size(class_rest{1,1},1)*size(class_rest{1,1},2)
   data2(1,i)=0;
   data2(j+1:j+1+size(class_rest{1,1},2)-1,i)=class_rest{1,i}(counter,:);
   dataeeg2(1,i)=0;
   dataeeg2(j+1:j+1+size(class_eegrest{1,1},2)-1,i)=class_eegrest{1,i}(counter,:);
   
   counter=counter+1;
   end
end
data5=[];
for  i=1:size(class_rest_power,2)
   counter=1;
    for j=1:size(class_rest_power{1,1},2):size(class_rest_power{1,1},1)*size(class_rest_power{1,1},2)
   data5(1,i)=0;
   data5(j+1:j+1+size(class_rest_power{1,1},2)-1,i)=class_rest_power{1,i}(counter,:);
   counter=counter+1;
   end
end
% data2=([class_rest{1,1};class_rest{1,2};class_rest{1,3};class_rest{1,4};class_rest{1,5};class_rest{1,6};class_rest{1,7};class_rest{1,8};class_rest{1,9};...
%     class_rest{1,10};class_rest{1,11};class_rest{1,12};class_rest{1,13};...
%     class_rest{1,14};class_rest{1,15};class_rest{1,16};class_rest{1,17};...
%     class_rest{1,18};class_rest{1,19};class_rest{1,20};class_rest{1,21};...
%     class_rest{1,22};class_rest{1,23};class_rest{1,24};class_rest{1,25};class_rest{1,26};...
%     class_rest{1,27};class_rest{1,28};...
%     class_rest{1,29};class_rest{1,30};class_rest{1,31};class_rest{1,32};...
%     class_rest{1,33};class_rest{1,34};class_rest{1,35};class_rest{1,36};class_rest{1,37};...
%     class_rest{1,38};class_rest{1,39};class_rest{1,40};class_rest{1,41};class_rest{1,42};...
%     class_rest{1,43};class_rest{1,44};class_rest{1,45};class_rest{1,46};class_rest{1,47};...
%     class_rest{1,48};class_rest{1,49};class_rest{1,50};class_rest{1,51};class_rest{1,52};...
%     class_rest{1,53};class_rest{1,54};class_rest{1,55};class_rest{1,56}]);
% data2=data2';
% data7=[];
% for i=1:size(data2,2)
% data7(1,i)=0;
% data7(2:size(data2,1)+1,i)=data3(:,i);    
% end
% data=[data6 data7];
data=[];
data=[data3 data2];
dataeeg= [dataeeg3 dataeeg2];
data_power=[data4 data5];
save data_forclassify_S1_L1_lap  data
save dataeeg_forclassify_S1_L1_lap  dataeeg
save data_forclassify_S1_L1_power_lap data_power
else
    load data_forclassify_S1_L1_lap
%     load  data_forclassify_S1_L1_power_etep
end

%%
%classification
clear
clc
% load dataeeg_forclassify_S1_L1_lapme

load data_forclassify_S1_L1_lap%etep%data_forclassify_S1_L1_lap
data1=data;
% load data_forclassify_S1_L1_lap
% data2=data;
% % load  data_forclassify_S1_L1_power_et
%  data=[data1 data2];
fs=256;
if size(data,1)>10*(1.5*fs+1)+1;
x1=data(2:(10*(1.5*fs+1)+1),:);%first row is classes 
else
    x1=data(2:(10*(1.5*fs+1)+1),:);
end
wavelet1=0;
if wavelet1==1
wave=1*DWT(1.5*fs+1);
% figure;plot(wave(1,:))
len=round((size(x1,1)-1)/(1.5*fs+1));
 
for i=1:(1.5*fs+1):len*(1.5*fs+1)
    for k=1:size(x1,2)
x(i:i+1.5*fs,k)=wave*x1(i:i+1.5*fs,k);
    end
end
else 
    x=x1;
end
len=round((size(x1,1)-1)/(1.5*fs+1));

count=1;k=1;
for i=1:(1.5*fs+1):len*(1.5*fs+1)
   
% figure(4);hold on;
% % subplot(10,1,count) ;
% plot(x1(i:i+1.5*fs,k)/norm(x1(i:i+1.5*fs,k)));
source(count,:)=x1(i:i+1.5*fs,k);
    count=count+1;

end
figure(4);plot(-1.5:1.5/(size(source,2)-1):0,source');
title('source activities of primary motor cortex for one trial of 1.5s before movementonset')
xlabel('time (s)');ylabel('current sources(UA)')

for k1=1:size(x,2);x(:,k1)=x(:,k1)/norm(x(:,k1));end
y=data(1,:);
% datawav=[data(1,:);x];
% save datawav_forclassify_S1_L1_etep  datawav
if 0
figure;plot(x1(:,1)/norm(x1(:,1)));hold on;plot(x(:,1),'r');for i=1:size(x1,1)/(1.5*fs+1);markersource(1,i)=i*(1.5*fs+1);end
EEGmarkersource(1,markersource)=1;
hold on;plot(0.1*EEGmarkersource,'k--')
legend('activities of sources on primary motor cortex for one trial of movement','sources in waveletdomain','marker to separate activities of sources ')
end


fs=256;
power=0;
if power==1;n_channels=round((size(x,1))/5);
else
n_channels=round((size(x,1))/round(1.5*fs));
end
%prepare index sets for cross-validation
n_permutations = 2;
n_epochs = size(x,2);%942
testsetsize = round(n_epochs / 10);
[trainsets, testsets] = crossValidation(1:n_epochs, testsetsize,...
n_permutations);                                      
 
est_y=[]
correct = [];figure(1);
for i = 1:n_permutations
   
    % draw data from CV index sets
    train_x = x(:, trainsets(i,:));%x(1:2944,848random numbers between1:942 defines the columns ) it defines which columns to be train
% it means we train with 848 signals of length( 23 channels*128 epochs)
    train_y = y(:, trainsets(i,:));%y is logical output for x
    test_x  = x(:, testsets(i,:));
    test_y  = y(:, testsets(i,:));    
    
    % train classifier and apply to test data
    l =LogitBoost(50, 0.05, 1);%LogitBoost number of iterations=50,
    l = train(l, train_x, train_y, n_channels);%n_channels=23
    p = classify(l, test_x); 
%********************************
 
 
    % evaluate classification accuracy 
    i0 = find(p <= 0.5);
    i1 = find(p > 0.5);
    est_y = zeros(size(p));
    est_y(i0) = 0;
    est_y(i1) = 1;
    for j = 1:size(est_y,1)
        n_correct(j) = length(find(est_y(j,:) == test_y));% each electrode is concidered for all testing epochs to find which one has the best feature
        %n_correct=n_correct1(j)+n_correct;
    end
    p_correct = n_correct / size(est_y,2);% it's an array (not a number) containing number of corrects
    correct = [correct ; p_correct];    
    
    % plot number of steps vs. classification accuracy 
    if (i>1)
        plot(mean(correct));
        xlabel('number of boosting iterations');
        ylabel('classification accuracy');
        drawnow;
    end
    
end
hold on
 
% counttrue=0;%number of true answers with atleast 85% accuracy  
% for i2=1:round(1/10*size(x,2))%94 %number of test_input epochs or signals(each signal has23*128samples or features)
% % disp('enter i2='); testsets(i,i2)
% pause(0.5);
%   p2=round(classify(l,x(:,testsets(i,i2))));
%   %answer to a specific data with 50 iteration
%   y2=length(find(p2==y(testsets(i,i2)) ));
%   % number of true answers of classifier for each data for 50 iteration of classification
% %   pause(0.5);
%   if y2>=.85*50% 50 is number of iterationsor answers of logitboost for each input
%       counttrue=counttrue+1;
%   end
% end
% disp('number of true answers with atleast 85% accuracy=')
% counttrue
% counttrue/(0.1*size(x,2))*100
% 
% num=10;
% y(1,num);
% p3=round(classify(l,x(:,num)));
% y3=length(find(p3==y(1,num)))%of 50 iteration

% save logit_L1_lap l

%%
%TP
num1=find(y==1);
numberpositive=size(num1,2);
counter=1;
for i=1:size(num1,2);
 p3=round(classify(l,x(:,num1(1,i))));
y3=length(find(p3==y(1,num1(1,i))));%of 50 iteration
if y3>=0.85*50;counter=counter+1;end
end
sensitivity=counter/size(num1,2)*100
%%%%%%selectivity=numberof TP/number of detections as positive

counterdetection=0;
for i=1:size(y,2);
 p3=round(classify(l,x(:,i)));
 y4=length(find(p3==1));if y4>=0.85*50;counterdetection=counterdetection+1;end

end
selectivity=counter/counterdetection*100

%TN
num0=find(y==0);
size(num0,2);
counter=1;
for i=1:size(num0,2);
 p3=round(classify(l,x(:,num0(1,i))));
y3=length(find(p3==y(1,num0(1,i))));%of 50 iteration
if y3>=0.85*50;counter=counter+1;end
end
specificity=counter/size(num0,2)*100
% accuracy=numberofTP/(totalnumber of detectionofposiives)
 save logitnew l
%%
%latency
clear 
eegselect=1;
if eegselect==1
load  dataeeg_forclassify_S1_L1_lap
data1=dataeeg;
else
load allsourcetrial_S1_L1_etep% allsourcetrial lap or etep
data1=data;
end

load markernumber_S1_L1
load EEG_marker_S1_L1
load logit_etep
% load  logitnew
% for n=1:size (markernumber,2)
% a=3
% countertrial1=round((markernumber(1,n)-1)*2/fs)+1;
% %markernumber(1,n) is time sample of movement onset
% 
% figure(8);hold on;subplot(4,7,n)
% plot(allsourcetrial{1,countertrial1-a}');
% end
counter=1;
fs=256;
wave=1*DWT(1.5*fs+1);

 for i=2:length(allsourcetrial)
i
     count=1;
    for k=1:size(allsourcetrial{1,1},1)
        x(count:count+size(allsourcetrial{1,1},2)-1,1)=allsourcetrial{1,i}(k,:)/norm(allsourcetrial{1,i}(k,:));
% % % % % % % % % % % % % % % %        x(count:count+size(allsourcetrial{1,1},2)-1,1)=wave*x(count:count+size(allsourcetrial{1,1},2)-1,1);
       x(count:count+size(allsourcetrial{1,1},2)-1,1)=x(count:count+size(allsourcetrial{1,1},2)-1,1)/norm(x(count:count+size(allsourcetrial{1,1},2)-1,1));
      count=count+size(allsourcetrial{1,1},2);
    end
    
 x=x/norm(x);
 p4=round(classify(l,x));
%  if isnan(p4)==1;disp('alarm nan');end
for k2=1:size(p4,1);
if isnan(p4(k2,1))==1;p4(k2,1)=0;end
end 
 
 positive=length(find(p4==1))
 negative=length(find(p4==0));
 if positive>negative;prob(1,i)=1;end
 if negative>positive;prob(1,i)=0;end
if positive==negative;prob(1,i)=0;end
 if prob(1,i)==1;
%        if positive>=.8*50
     markernumberdetection(1,counter)=round((i-1)*fs/2+1) ;
    EEG_markerdetection(1,markernumberdetection(1,counter))=1;
    counter=counter+1;
%    end
 end
end
figure;stem((1:(length(EEG_markerdetection)))/60/fs,1/2*EEG_markerdetection,'m');hold on;plot((1:(length(EEG_marker)))/60/fs,EEG_marker);legend('detectiontime','movementonset')
xlabel('time minute')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% n1=150;
% n=300;
% end1=(n-1)*fs/2+1+1.5*fs;
% for i=n1:n
% figure(5);hold on;plot(((i-1)*fs/2+1:(i-1)*fs/2+1+1.5*fs)/fs,(allsourcetrial{1,i}(1,:)));%/norm(allsourcetrial{1,i}(1,:))
% end
% hold on;stem(((n1-1)*fs/2+1:end1)/fs,100*EEG_markerdetection(1,(n1-1)*fs/2+1:end1),'bO')
% hold on;stem(((n1-1)*fs/2+1:end1)/fs,100*EEG_marker(1,(n1-1)*fs/2+1:end1),'r*')
% xlabel('time(s) ');legend('red star is movementonset and blue circle is movementdetection')
% selectivity=TP/total number of detection(as positive)=20/counter*100
%  sensitivity=TP/positive=(20)/28*100
%%
counterdetection=0;
totallatency=0;
distance=[];
latency=[];
totallatency=0;
for i=1:length(markernumberdetection)
    for k1=1:length(markernumber) 
        dis2=markernumber(1,k1)- markernumberdetection(1,i);
    distance(1,k1)=abs(dis2)/fs   
    end
   num=find(distance==min(distance));
   if distance(1,num(1,1))<=3;
      counterdetection=counterdetection+1;
       latency(1,counterdetection)=(markernumberdetection(1,i)-markernumber(1,num(1,1)))/fs;
   totallatency=totallatency+ latency(1,counterdetection);
   end
end
totallatency=totallatency /counterdetection*1000 
disp(totallatency );disp('ms') 
% % % % % % % % % % % selectivity=counterdetection/length(markernumberdetection)*100
% % % % % % % % % % % sensitivity=counterdetection/length(markernumber)*100
% % % % % % % % % % % % TPR=sensitivity;
% % % % % % % % % % % % FP=(length(markernumber)-counterdetection)/length(markernumber)*100
figure;stem(latency/10);ylabel('sec');title('latency');xlabel('trials')
