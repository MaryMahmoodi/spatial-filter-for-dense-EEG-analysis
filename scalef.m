%scalpef.m ? Spline estimation of the components of the scalp Electric Field
 %
 % Usage: [Et,Ep,L,S] = scalpef( locs , m , lambda );
 %
 % Required Inputs:
 % locs = sensor locations in Cartesian coordinates ([X Y Z])
 % m = interpolation order (small integer larger than 2)
 % lambda = smoothing parameter
 %
 % Output:
 % Et, Ep, L = spherical components of the tangential electric field and the surface Laplacian
 % S = smoothing matrix
 %
 % Author: Claudio G. Carvalhaes, May 01, 2013
 %

function [Et,Ep,L,S] = scalpef( locs , m , lambda);% , data , channels )

 % Handle arguments
 if nargin~=3
 help ssef.m;
 return;
 end

 [N,d] = size(locs);

 if m < 3
 error('scalpef:m','m should be greater than 3.');
 end

 if d~=3
 error('scalpef:locs','LOCS must be a Nx3 matrix.');
 end

 M = nchoosek(m+2,d);
 if M>=N
 error('scalpef:M','M (=m+2 choose 3) must be less than %d.',N);
 end

 if lambda<0 || ~isscalar(lambda)
 error('scalpef:lambda','LAMBDA must be a non?negative scalar.');
 end

 h = sqrt(mean(dot(locs,locs,2)));
 locs = h*normr(locs);

 [theta,phi] = cart_to_sph(locs(:,1),locs(:,2),locs(:,3));

 ct = cos(theta);
 st = sin(theta);
 cp = cos(phi);
 sp = sin(phi);

 cpp = cos(bsxfun(@minus,phi,phi'));
 spp = sin(bsxfun(@minus,phi,phi'));

 cos_gamma = ct*ct' + (st*st').*cpp;
 cos_gamma = max(cos_gamma,-1);
 cos_gamma = min(cos_gamma,1);

 K = (2^(m-3/2))*(h^(2*m-3))*((1-cos_gamma).^(m-3/2));
 K_Et = (2^(m-3/2))*(h^(2*m-4))*(m-3/2)*((1-cos_gamma).^(m-5/2)).*((ct*st').*cpp - st*ct');
 K_Ep = -(2^(m-3/2))*(h^(2*m-4))*(m-3/2)*((1-cos_gamma).^(m-5/2)).*(ones(N,1)*st').*spp;
 K_Lap = (2^(m-3/2))*(h^(2*m-5))*(m-3/2)*((1-cos_gamma).^(m-5/2)).*((m-1/2)*cos_gamma + m - 5/2);

 b = combnk(repmat(0:m-1,1,3),3);
 t = sort(b,2,'descend')==b;
 n = b(all(t==1,2),:)'; % b = [i;j;k]

 i = repmat(n(1,:),N,1);
 j = repmat(n(2,:),N,1);
 k = repmat(n(3,:),N,1);

 st = repmat(st,1,M);
 ct = repmat(ct,1,M);
 sp = repmat(sp,1,M);
 cp = repmat(cp,1,M);

 st2 = st.*st;
 ct2 = ct.*ct;
 ct4 = ct2.*ct2;
 sp2 = sp.*sp;
 sp4 = sp2.*sp2;
 cp2 = cp.*cp;
 cp4 = cp2.*cp2;

 st2cp2 = st2.*cp2;
 st2sp2 = st2.*sp2;
 sp2cp2 = sp2.*cp2;
 ct2sp2 = ct2.*sp2;

 % The polynomial term and its derivatives
 % For m=3: T = [1 x y z x^2 x*y x*z y^2 y*z z^2]

 subs = @(F,x) F(x);
 hi1 = h.^(i-1);
 T = h*hi1.*(st.^(i-k)).*(ct.^k).*(sp.^(j-k)).*(cp.^(i-j));
 T_Et = -hi1.*(st.^(i-k-1)).*(ct.^(k-1)).*(sp.^(j-k)).*(cp.^(i-j)).*(i.*ct2 - k);
 T_Et(k==i) = subs(i.*hi1.*st.*(ct.^(i-1)),k==i);
 T_Et(k==0) = -subs(i.*hi1.*(st.^(i-1)).*ct.*(sp.^j).*(cp.^(i-j)),k==0);
 T_Et(i==0) = 0;
 T_Ep = -hi1.*(st.^(i-k-1)).*(ct.^k).*(sp.^(j-k-1)).*(cp.^(i-j-1)).*(j - i.*sp2 - k.*cp2);
 T_Ep(j==i) = -subs((i-k).*hi1.*(st.^(i-k-1)).*(ct.^k).*(sp.^(i-k-1)).*cp,j==i);
 T_Ep(k==j) = subs((i-j).*hi1.*(st.^(i-j-1)).*(ct.^j).*sp.*(cp.^(i-j-1)),k==j);
 T_Ep(j==i & k==j) = 0;
 ii1 = i.*(i+1);
 hi2 = hi1/h;
 T_Lap = hi2.*(st.^(i-k-2)).*(ct.^(k-2)).*(sp.^(j-k-2)).*(cp.^(i-j-2))...
 .* (sp2cp2 .* (ii1.*ct4 + k.*(k-1)) + ct2.*( (i-j).*(i-j-1).*sp4 ...
 + (j-k).*(j-k-1).*cp4 - 2.*(i + i.*j - j.*j + j.*k - k).*sp2cp2 ));

 evtl = @(f,e) evalc('T_Lap(e) = f(e)');
 evtl(hi2.*(st.^(i-k-2)).*(ct.^(k-2)).*(sp.^(i-k-2)) ...
 .*(sp2.*(ii1.*ct4 + k.*k - k) + ct2.*((i-k).*(i-k-1).*cp2 - 2*(i + i.*k - k).*sp2)) , j==i);
 evtl(hi2.*(st.^(i-k-2)).*(ct.^(k-2)).*(sp.^(i-k-3)).*cp ...
 .*(sp2.*(ii1.*ct4 + k.*k - k) + ct2.*((i-k-1).*(i-k-2).*cp2 - 2*(2*i + i.*k - 2.*k - 1).*sp2 ...
)),...
 j==i-1);
 evtl(i.*hi2.*(ct.^(i-2)).*(i.*st2 - ct2 -1) , k==i);
 evtl(hi2.*(st.^(i-j-2)).*(ct.^(j-2)).*(cp.^(i-j-2))...
 .*(cp2.*(ii1.*ct4 + j.*j - j) + ct2.*( (i-j).*(i-j-1).*sp2 - 2*(i + i.*j - j).*cp2 )) , k==j);
 evtl(hi2.*(st.^(i-j-1)).*(ct.^(j-3)).*sp.*(cp.^(i-j-2)) ...
 .*(cp2.*(ii1.*ct4 + (j-1).*(j-2)) + ct2.*( (i-j).*(i-j-1).*sp2 - 2*(i + i.*j - 2*j + 1).*cp2 ...
)),...
 k==j-1);
 evtl(hi2.*(st.^(i-2)).*(sp.^(j-2)).*(cp.^(i-j-2))...
 .*(ii1.*ct2sp2.*cp2 + (i-j).*(i-j-1).*sp4 + j.*(j-1).*cp4 - 2*(i + i.*j - j.*j).*sp2cp2 ) , ...
k==0);
 evtl(hi2.*(st.^(i-3)).*ct.*(sp.^(j-3)).*(cp.^(i-j-2))...
 .*(ii1.*ct2sp2.*cp2 + (i-j).*(i-j-1).*sp4 + (j-1).*(j-2).*cp4 - 2*(i + i.*j - j.*j + j ...
-1).*sp2cp2),...
 k==1);
 evtl(hi2.*st.*(ct.^(i-3)).*sp.*(ii1.*st2 - 4*i+2) , j==i & k==i-1);
 evtl(hi2.*st.*(ct.^(i-3)).*cp.*(ii1.*st2 - 4*i+2), j==i-1 & k==j);
 evtl(hi2.*st2.*(ct.^(i-4)).*sp.*cp.*(ii1.*st2 - 6*i + 6) , j==i-1 & k==i-2);
 evtl(hi2.*(st.^(i-2)).*(sp.^(i-3)).*cp.*((i-1).*(i-2) -ii1.*st2sp2) , j==i-1 & k==0);
 evtl(i.*hi2.*(st.^(i-2)).*(sp.^(i-2)).*(i-1 - (i+1).*st2sp2) , j==i & k==0);
 evtl(hi2.*(st.^(i-3)).*ct.*(sp.^(i-3)).*(ii1.*(cp2 + ct2sp2) - 4*i + 2), j==i & k==1);
 evtl(hi2.*(st.^(i-3)).*ct.*(sp.^(i-4)).*cp.*(ii1.*(cp2 + ct2sp2) - 6*i + 6), j==i-1 & k==1);
 evtl(i.*hi2.*(st.^(i-2)).*(cp.^(i-2)).*(i-1 - (i+1).*st2cp2), j==0 & k==0);
 evtl(hi2.*(st.^(i-2)).*sp.*(cp.^(i-3)).*((i-1).*(i-2) - ii1.*st2cp2) , j==1 & k==0);
 evtl(hi2.*(st.^(i-3)).*ct.*(cp.^(i-3)).*((i-1).*(i-2) - ii1.*st2cp2) , j==1 & k==1);
 evtl(hi2.*(st.^(i-3)).*ct.*sp.*(cp.^(i-4)).*((i-2).*(i-3) - ii1.*st2cp2), j==2 & k==1);
 T_Lap(i==0) = 0;
 evtl(-2*st.*cp/h, i==1 & j==0 & k==0);
 evtl(-2*st.*sp/h, i==1 & j==1 & k==0);
 evtl(-2*ct/h, i==1 & j==1 & k==1);
 evtl(-6*st2.*sp.*cp, i==2 & j==1 & k==0);
 evtl(-6*st.*ct.*cp, i==2 & j==1 & k==1);
 evtl(-6*st.*ct.*sp, i==2 & j==2 & k==1);
 evtl(-12*h*st2.*ct.*sp.*cp,i==3 & j==2 & k==1);

 [Q1,Q2,R] = qr2(T);
 C = Q1/(Q1'*(K + N*lambda*eye(N))*Q1)*Q1';% primary C = Q2*pinv(Q2'*(K + N*lambda*eye(N))*Q2)*Q2';
%changed by mary
[U1,S1,V1]=svd(C);
[eigv,eigval]=eig(C);
[a,b]=find(eigval==max(max(eigval)));
landa=max(max(eigval));
C=C+0.03*landa*U1;
 
D = pinv(R)*(eye(N) - K*C - N*lambda*C);%primary D = pinv(R)*Q1'*(eye(N) - K*C - N*lambda*C);
%changed by mary
%to solve rank deficient we add unitary matrix to each matrix
[U1,S1,V1]=svd(D);
[eigv,eigval]=eig(V1);
[a,b]=find(eigval==max(max(eigval)));
landa=real(max(max(eigval)));
D=D+0.03*landa*S1;
 


S = K*C + T*D; % The smoothing matrix
 Et = K_Et*C + T_Et*D;
 Ep = K_Ep*C + T_Ep*D;
 L = K_Lap*C + T_Lap*D;

 end
