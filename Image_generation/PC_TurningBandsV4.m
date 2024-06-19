%Author: Frederic Richard
%Institution: Aix-Marseille University, LATP (UMR CNRS 6632).
%Date:january, 8th, 2012
%Version:4
%Developed and tested on qoctave version 0.10.1

%Simulation of anisotropic fractional Brownian fields
%by a Turning Band method

%The field X to be simulated is Gaussian with stationary increments,
%and its variogram is characterized by a spectral density of the form
%f(w)=c(arg(w)) |w|^h(arg(w)), w in R^2,
%where c and h are pi-periodic functions which depend both on the direction arg(w) of w
%and are both piecewise constant on the unit sphere of R^2.

%Note: in this version, orientations of bands are sampled using dynamic programming to satisfy a tradeoff between approximation accuracy and computational cost.

%Reference: A turning-band method for the simulation of anisotropic
%fractional Brownian fields, by H. Bierme, L. Moisan and F. Richard, 2012


%INPUT
%N: NxN = image size
%K determines the accuracy of the simulation which is measured as pi/K.
%ang, c, h: h(i+1) and c(i+1) are respective values of h and c on interval [ang(i),ang(i+1)[ 
%The [ang(i),ang(i+1)[ must form a partition of [-pi/2,pi/2[. In particular ang(1)=-pi/2 and ang(end)=pi/2. 
%display: option for displaying results 
%coord: a set of coordinates where to simulate the field (the coordinates must be entered as positive integers)
%If this argument is missing, the field is simulated 
%on the grid {0,1/N,...,(N-1)/N}x{0,1/N,...,(N-1)/N} by default. 
%vario: option for computing the variogram of the synthetic field.

%OUTPUT
%X: simulated field
%Kangle,Pangle,Qangle: angles used for defining turning bands (Kangle=angle, tan(Kangle)=Qangle/Pangle).
%coordx,coordy: cartesian coordinates where the field was simulated
%v: variogram of the simulation field


% For the simulation of a fractional Brownian field
% of Hurst index H (=0.5),
% [ang,cexpe,hexpe] = PC_h_constant(H,pi/2);
% X=PC_TurningBandsV4(512,1000,ang,cexpe,hexpe,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% MAIN FUNCTION %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [X,coordx,coordy,v] = PC_TurningBandsV4(N,K,ang,c,h,display,vario,coord)


if nargin<1, N=512; end
if nargin<2, K=1000; end
if nargin<3, ang=[-pi/2 pi/2]; end
if nargin<4, c=1; end
if nargin<5, h=0.5; end
if nargin<6, display=1; end
if nargin<7, vario=0; end


%

%I. INITIALIZATION OF TURNING BANDS

%I.A) definition of coordinates where to simulate the field
%By default, if the input coord is not given, 
%the coordinates are on a uniform grid of [0,1]x[0,1].
if nargin<8,
    %Definition of the uniform grid
    coord=repmat(0:(N-1),N,1);
    coord=[reshape(coord,N^2,1),reshape(coord',N^2,1)];
else
    %Definition of coordinates as given by coord.
    N=floor(max(coord(:))/2)*2+2;
end

%I.b) definition of band orientations with rational slopes using dynamic programming.
[Kangle, Pangle, Qangle]=BuildOptVec(K); 

%II. SIMULATION BY TURNING BANDS
X=zeros(size(coord,1),1);
if (vario==1), v=zeros(size(coord,1),1); else v=[]; end
K=length(Kangle);
for k=2:K,
    kang=Kangle(k);
    p=Pangle(k); q=Qangle(k); 
    nbsamp=(N-1)*(p+abs(q))+1; %number of samples to be simulated on the turning band
    [Cangle,H]=Find_parameters_angle(ang,h,c,Kangle(k)); %values of parameters h and c at angle Kangle
    if (Cangle~=0)
        %Simulation of fractional Brownian motion of order H on the band
        [Y,vl]=FBMsimulation(H,nbsamp);
       %Computations of the projection coordinates of grid points on the band
        indr=p*coord(:,1)+q*coord(:,2);   
        %Weight the process on the turning band (based on rectangular integral approximation of the variogram)
        weig=Weight_TB(p,N,Cangle,H,Kangle(k-1),kang); 
        %Computation of the simulation field variogram (optional)
        if (vario==1), v=v+0.5*weig^2*abs(indr).^(2*H); end
 
        %Update the simulation field
        indr=indr-min(indr)+(1); 
        X=X+weig*Y(sub2ind(size(Y),indr));
    end
end


if nargin<8
   %reshape vectors in an image format
   coordx=reshape(coord(:,1),N,N)/N;
   coordy=reshape(coord(:,2),N,N)/N;
   X=reshape(X,N,N);
   if (vario==1), v=reshape(v,N,N); end
else 
   coordx=coord(:,1)/N;
   coordy=coord(:,2)/N;
end

if nargin<8, 
    if display,
        disp(['Simulation with ' num2str(K) ' turning bands'])
        figure(display), clf
        colormap(gray)
        imagesc(X);
        axis xy, axis off, axis equal, colorbar
    end
end

%Computation of the weight of a turning band using a rectangular integral approximation
 function weig=Weight_TB(p,N,c,h,ang0,ang1)   
       weig=sqrt((ang1-ang0)*Kval(h)*c);
        if p~=0, 
                weig=weig*(cos(ang1)/(N*p))^h;
        else, 
               weig=weig/N^h; 
        end
        
%Finding the values of parameters c and h of the fields at a given orientation.
%The parameters c and h of the fields are entered via ang, h and c. The given orientation
%is kang
function [cangle, hangle] = Find_parameters_angle(ang,h,c,kang)
   ang=[ang, pi/2+ang(2)-ang(1)]; c=[c, c(1)]; h=[h, h(1)];
   aux=sort(find(ang>kang));
   %value of the topothesis constant C of the Brownian motion simulated on the lth band
   cangle=c(aux(1)-1);
   %value of the Hurst index H of the Brownian motion simulated on the lth band
   hangle=h(aux(1)-1);

function [Kangle, Pangle, Qangle] = BuildOptVec(N,display)

if nargin<1, N=1000; end
if nargin<2, display=0; end
prec=pi/N;
Nmax=30000;
% Definition of a set S of possible angles and their associated costs.
if N<Nmax/10, [Kangle, Pangle, Qangle] = Sampling_uniformly_angles(Nmax);
else, [Kangle, Pangle, Qangle] = Sampling_uniformly_angles(10*N); end
% Possible angles s are chosen to be approximately
% sampled on [-pi/2,pi/2] while satisfying tan(s)=q/p for q in Z and p in N.
% To save time, the set S can be loaded as follows
%load('uniformly-sampled-rat-angles.mat');



Cost=abs(Qangle)+Pangle;
nvec=length(Kangle);

%Dynamic programming for the selection of an optimal subset
cost = zeros(nvec,1); %partial costs
pos  = zeros(nvec,1);
cost(nvec) = Cost(nvec);
pos(nvec) = nvec;
for i=(nvec-1):-1:1
    bound = Kangle(i)+prec;
    bestj = i+1;
    mini = cost(bestj);
    %Seek among upper angles at distance below prec of the current angle
    %for the one with minimal partial cost.
    for j=(i+2):nvec
        if Kangle(j)>bound, break; end
        if cost(j)<mini
            bestj = j;
            mini = cost(bestj);
        end
    end
    %Define the partial cost of the current angle and the best associated upper angle.
    cost(i) = Cost(i)+mini;
    pos(i) = bestj;
end

%Build the set S' by finding the best path.
i=1;
Select=zeros(nvec,1);
while i<nvec
    Select(i)=1;
    i=pos(i);
end
Select(nvec)=1;

Kangle=Kangle(Select==1);
Pangle=Pangle(Select==1);
Qangle=Qangle(Select==1);
nvec=length(Qangle);
acc=max(diff(Kangle)); %accuracy associated to the selected set of angles



%Function for building the set S of possible angles
%S is formed with K angles which are approximately sampled on [-pi/2,pi/2]
%and with tangents satisfying  tan(s)=q/p, where q is in Z and p in N.

function [Kangle, Pangle, Qangle, Cost] = Sampling_uniformly_angles(K)

if nargin<1, K=100000; end

Kangle = linspace(-pi/2,pi/2,K);
K=length(Kangle);
Pangle = zeros(1,K); Qangle = zeros(1,K); Cost=zeros(1,K);
delta=Kangle(2)+pi/2;
for l=2:(K-1),
    prec=delta*(1+tan(Kangle(l))^2)*0.5;
    %rational approximation of the slope tan(Kangle(l)) of orientation Kangle(l)
    %at precision prec
    [q,p] = rat(tan(Kangle(l)),prec);
    %update the orientation
    Kangle(l)=atan2(q,p);
    Pangle(l)=p; Qangle(l)=q;
end
%plot(1:length(Kangle),Kangle);
Pangle(1)=-1; Qangle(1)=0; Cost(1)=1;
Pangle(end)=1; Qangle(end)=0; Cost(end)=1;
%removal of possible doublons
[Kangle, I]=unique(Kangle);
nvec=length(Kangle);
Pangle=Pangle(I); Qangle=Qangle(I);
save('uniformly-sampled-rat-angles.mat','Kangle','Qangle','Pangle');

% non-negative greatest common divisor of two integers (positive or not)
function v = gcd(x,y)

if (x<y),  v=gcd(y,x); return; end
if (y<0),  v=gcd(x,-y); return; end
if (y==0), v=x; return; end
v=mod(x,y);




% Computation of the weight of a turning band using a rectangular integral 
% approximation
%INPUT:
%H: Hurst index
%OUTPUT :
%val : quantity part of the weight of each turning band
function val=Kval(H)
val=pi/(H*gamma(2*H)*sin(H*pi));




%Simulation of a fractional Brownian motion of Hurst index H
%Method described in E. Perrin, R. Harba, R. Jennane, and I. Irribaren,
%Fast and Exact Synthesis for 1-D fractional Brownian motion
%and fractional Gaussian noise, IEEE SIGNAL PROCESSING LETTERS,
%9(11):382-384, NOVEMBER 2002
%INPUT:
%H: Hurst index
%nbsamp: number of samples
%OUTPUT
%y: a realization of a FBM of Hurst index H
%on points of the uniform grid 0,1,...,nbsamp-1
%s: variogram of the Gaussian fractional noise on 0,1,...,nbsamp-1
function [y,s]=FBMsimulation(H,nbsamp,s)

if nargin<3,
    %Definition of an embedded variogram of the fractional Brownian motion of Hurst index H
    k=0:(nbsamp-1);
    H2=2*H;
    s=0.5*((k+ones(size(k))).^H2+abs(k-ones(size(k))).^H2-2*abs(k).^H2);
    s=[s,fliplr(s(2:(length(s)-1)))]';
    nbsamp2=length(s);
    %Computation of square-root of spectral density.
    s=fft(s);
    if (min(real(s(:)))<=0), disp('warning: negative eigenvalue'); end
    s=sqrt(s);
end


%Realization of a fractional Gaussian noise
y=real(fft(s.*(randn(nbsamp2,1)+sqrt(-1)*randn(nbsamp2,1)))/sqrt(nbsamp2));

%Realization of a fractional Brownian motion
y=y(1:nbsamp);
y=cumsum(real(y)); %y=y-(y(1));


