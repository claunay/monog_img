% Script to analyze one image using monogenic and structure tensors
% Estimation of the Hurst parameter and the coherence index of a lighthouse
% field


close all
rng(0)          % uncomment for reproducibility

M=1024;         % Size of the realization
H=0.5;          % Hurst parameter
alp=pi/3;       % Degree of anisotropy (If alp=pi/2, X is an isotropic FBF)
or = 0;%pi/8;   % Rotation offset of the field
L = 7;          % Depth (number of scales) of the analysis
scale = [3,4];  % Choice of the scales for the estimations


% Image generation
[ang,c,h] = PC_h_constant_or(H,alp,or);
[X] = PC_TurningBandsV4(M,500,ang,c,h,0);


%Multiscale analysis
bord=floor(M*15/100); %border correction
[orient,AM,phase,Dir,Vap,coherence,Sigma,Corr,dim] = calcul_mwt_classique(X,L,bord); 


%Coherence estimation
coherence_est = mean(coherence(scale));


%Hurst parameter estimation using the squared amplitude of the monogenic tensor
AM_sq_est = reshape(sum(AM.^2,[1,2])/((M-2*bord+1)^2),[1,L]);
Lreg = [ones(length(scale),1), (scale-1)'];
reg_amp = Lreg\ log(AM_sq_est(scale)');
H_est_amp = (1/2)*(reg_amp(2)/log(2)-2);


%Hurst parameter estimation using the structure tensor
Vapmax=reshape(Vap(2,2,:),1,L);
Vapmin=reshape(Vap(1,1,:),1,L);
Lreg = [ones(length(scale),1), zeros(length(scale),1), (scale-1)'; zeros(length(scale),1), ones(length(scale),1), (scale-1)'];
Vapreg = [log(Vapmin(scale))';log(Vapmax(scale))'];
reg_riesz = Lreg\ Vapreg;
H_est_riesz = (1/2)*(reg_riesz(3)/log(2)-2);



% Display results

figure; imagesc(X); colormap("gray")
title('One realization')

figure; 
for k=1:L
subplot(2,L,k)
imagesc(AM(:,:,k))
colorbar('horiz')
title(['Scale ',num2str(k)])
subplot(2,L,L+k)
hist(reshape(AM(:,:,k),dim(1)*dim(2),1),100)
end
sgtitle('Amplitude of the monogenic tensor')
 
 
figure; 
for k=1:L
subplot(2,L,k)
imagesc(orient(:,:,k))
colorbar('horiz')
title(['Scale ',num2str(k)])
subplot(2,L,L+k)
hist(reshape(orient(:,:,k),dim(1)*dim(2),1),100)
end
sgtitle('Orientation of the monogenic tensor')
 
  
figure; 
for k=1:L
subplot(2,L,k)
imagesc(phase(:,:,k))
colorbar('horiz')
title(['Scale ',num2str(k)])
subplot(2,L,L+k)
hist(reshape(phase(:,:,k),dim(1)*dim(2),1),100)
end
sgtitle('Phase of the monogenic tensor')


figure
hold on
plot(1:L,coherence_est*ones(L),'b:')
plot(1:L,coherence,'xb')
plot(scale,coherence(scale),'o','MarkerFaceColor','b','MarkerEdgeColor','b')
plot(1:L, (sin(2*alp)./(2*alp))'*ones(1,L), 'k:');
axis([1 L 0 1])
xlabel('scale')
ylabel('coherence')
xlabel('scale j')
title('Estimation of the coherence index')
hold off


figure
hold on
plot(1:L,log(AM_sq_est),'or')
plot(scale,log(AM_sq_est(scale)),'o','MarkerFaceColor','r','MarkerEdgeColor','r')
plot(1:L,(0:L-1)*reg_amp(2)+reg_amp(1),'r:');hold on
plot(1:L,(-2:L-3)*log(2)*(2*H+2)+log(AM_sq_est(3)),'k:')
xlabel('scale j')
ylabel('2log(V(u_j))')
title('Estimation of the Hurst index (Squared amplitude - Monogenic tensor)')
hold off


figure
hold on
plot(1:L,log(Vapmin),'or')
plot(scale,log(Vapmin(scale)),'o','MarkerFaceColor','r','MarkerEdgeColor','r')
plot(1:L,(0:L-1)*reg_riesz(3)+reg_riesz(1),'r:')
plot(1:L,log(Vapmax),'xb')
plot(scale,log(Vapmax(scale)),'o','MarkerFaceColor','b','MarkerEdgeColor','b')
plot(1:L,(0:L-1)*reg_riesz(3)+reg_riesz(2),'b:')
plot(1:L,(-2:L-3)*log(2)*(2*H+2)+log(Vapmin(3)),'k:')
plot(1:L,(-2:L-3)*log(2)*(2*H+2)+log(Vapmax(3)),'k:')
xlabel('scale j')
ylabel('log(\lambda^\pm(u_j))')
title('Estimation of the Hurst index (Structure tensor)')
hold off

