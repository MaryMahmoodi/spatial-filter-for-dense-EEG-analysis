function [ optmom ] = optimummomentum_lcmv_sekihara2015( C,lfm,EEG,pos,vol,momanatomy)
%it gets leadfieldmatrix(lfm)for a point in space and also covariance
%matrix (C) of the data in a trial or epoch and gives the best momentum for the
%source of that place also the neural activity index
%optmom.mom and %optmom.NAI  optmom.W %wight vector for lcmv
% This is according to sekihara book 2015 for lcmv beamformer
%and by this optimum orientation used in lcmv of fieldtrip we try it
%besides  this one to get an other NAI and weight vector for EEG
% vol is simbio vol; volume conductor model from simbio finite element model, it is only used in beamformer_lcmv.m 
% momanatomy: facenormal direction on each source location on the brain surface
%%
parameter1=lfm'*(inv(C))*lfm;
parameter2=lfm'*lfm;
%% define optimum direction with eigen value%%%
[eigv,eigval]=eig(parameter1);
for i=1:size(eigval,1);eigvalue(1,i)=eigval(i,i);end
a=find(eigvalue==min(eigvalue));
landa_param1=eigvalue(1,a);
eigvector_param1=eigv(:,a);
eigvector_param1=eigvector_param1/norm(eigvector_param1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[eigv,eigval]=eig(parameter2);
for i=1:size(eigval,1);eigvalue(1,i)=eigval(i,i);end
a=find(eigvalue==min(eigvalue));
landa_param2=eigvalue(1,a);

eigvector_param2=eigv(:,a);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if landa_param1<landa_param2
    optmom.mom=real(eigvector_param1);
% % % % % %     optmom.mom=mom1;
    
    l=lfm*optmom.mom;
        l=l/norm(l);
        optmom.l=l;
    optmom.W=(inv(C)* (optmom.l))/( (optmom.l)'*inv(C)* (optmom.l));
optmom.NAI=trace((optmom.W)'*(C)*optmom.W)%landa_param1;
 optmom.NAI=abs(1/landa_param1);

else
    optmom.mom=real(eigvector_param2);
% % % %         optmom.mom=mom1;

    l=lfm*optmom.mom;
        l=l/norm(l);
     optmom.l=l;
    optmom.W=(inv(C)* (optmom.l))/( (optmom.l)'*inv(C)* (optmom.l));
 optmom.NAI=(optmom.W)'*(C)*optmom.W%landa_param2;
 optmom.NAI=abs(1/landa_param2);

end
%%now find optmom.Wlcmv and  optmom.NAIlcmv for pos1 and mom1
if 0
grad=[];
dip2.pos=pos;
dip2.ori=optmom.mom;
dip2.inside=1;dip2.outside=[];

dip2.leadfield{1,1}=lfm;
dipout1 = beamformer_lcmv(dip2, grad, vol, EEG, C,'lambda',1,'powmethod','trace','feedback','text','fixedori' ,'no','projectnoise' ,'yes','projectmom','yes','keepfilter' ,'yes','keepleadfield','yes','keepmom','yes','keepcov','yes', 'reducerank' ,'no'  );
optmom.Wlcmv=dipout1.filter{1,1}';
optmom.NAIlcmv=abs(optmom.Wlcmv'*C*optmom.Wlcmv);

% figure;plot(optmom.Wlcmv'*EEG)
%try with momanatomy
dip2.ori=momanatomy;
dipout1 = beamformer_lcmv(dip2, grad, vol, EEG, C,'lambda',1,'powmethod','trace','feedback','text','fixedori' ,'no','projectnoise' ,'yes','projectmom','yes','keepfilter' ,'yes','keepleadfield','yes','keepmom','yes','keepcov','yes', 'reducerank' ,'no'  );
optmom.Wlcmvanatomy=dipout1.filter{1,1}';
optmom.NAIlcmvanatomy=abs(optmom.Wlcmv'*C*optmom.Wlcmv);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end%if 0
l=lfm*momanatomy;
      l=l/norm(l);
W=(inv(C)*l)/(l'*inv(C)*l);
optmom.Wmomanatomy=W;
optmom.NAIlcmvmomanatomy=abs((l'*l)/(l'*((C)\l)));



end

