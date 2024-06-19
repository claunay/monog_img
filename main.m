%Script to analyze one image using monogenic and structure tensors

close all
rng(0)  % uncomment for reproducibility

M=1024;
H=0.5; 
alp=pi/8;
or = pi/3;
L = 4;

[ang,c,h] = PC_h_constant_or(H,alp,or);
[X] = PC_TurningBandsV4(M,500,ang,c,h,0);
bord=floor(M*15/100); %border correction

[orient,AM,phase,Dir,Vap,coherence,Sigma,Corr,dim] = calcul_mwt_classique(X,L,bord); 

%Coherence estimation
coherence_est = mean(coherence);

%Hurst parameter estimation
Vapmax=reshape(Vap(2,2,:),1,L);
Vapmin=reshape(Vap(1,1,:),1,L);
Lreg = [ones(L,1), zeros(L,1), (0:L-1)'; zeros(L,1), ones(L,1), (0:L-1)'];
Vapreg = [log(Vapmin)';log(Vapmax)'];
reg = Lreg\ Vapreg;
H_est = (1/2)*(reg(3)/log(2)-2);

% Display results
figure; imagesc(X); colormap("gray")

figure; 
for k=1:L
subplot(2,L,k)
imagesc(AM(:,:,k))
colorbar('horiz')
subplot(2,L,L+k)
hist(reshape(AM(:,:,k),dim(1)*dim(2),1),100)
end
 title('Amplitude')
 
 
figure; 
for k=1:L
subplot(2,L,k)
imagesc(orient(:,:,k))
colorbar('horiz')
subplot(2,L,L+k)
hist(reshape(orient(:,:,k),dim(1)*dim(2),1),100)
end
 title('Orientation')
 
  
figure; 
for k=1:L
subplot(2,L,k)
imagesc(phase(:,:,k))
colorbar('horiz')
subplot(2,L,L+k)
hist(reshape(phase(:,:,k),dim(1)*dim(2),1),100)
end
 title('Phase')


figure
plot(0:L-1,coherence_est*ones(L),':b')
hold on
plot(0:L-1,coherence,'xb')
axis([0 L-1 0 1])
xlabel('scale')
ylabel('coherence')
xlabel('scale j')

figure
plot(0:L-1,log(Vapmin),'or')
hold on
plot(0:L-1,(0:L-1)*reg(3)+reg(1),'r:')
plot(0:L-1,log(Vapmax),'xb')
plot(0:L-1,(0:L-1)*reg(3)+reg(2),'b:')
xlabel('scale j')
ylabel('log(\lambda^\pm(u_j))')
