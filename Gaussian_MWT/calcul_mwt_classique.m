function [orient,AM,phase,Dir,Vap,coherence,Sigma,Corr,dim] = calcul_mwt_classique(X,L,bord)
    M = size(X,1);

    param.typ = {'GaussianHP'}; 
    param.norma = true;   % normalize subbands
    param.sampling = 'undec'; % 'undec' or 'pyramid'
    param.noRzHF = false; % remove Riesz transform at first scale; % monogenic wavelet transform
    
    img=X(1:M,1:M); 
    mwt = mwt_radial('a',img,L,param);

    for l=1:L
        s = mwt{l,1};
        sy= mwt{l,2};
        sx= mwt{l,3};

        % Border correction
        prim=s(bord:M-bord,bord:M-bord);riez1=sx(bord:M-bord,bord:M-bord);riez2=sy(bord:M-bord,bord:M-bord);
        dim = size(prim);

        riez=riez1+i*riez2;
       
        % Monogenic transform characteristics
        th = atan2(riez2,riez1);
        am = sqrt(prim.^2+riez1.^2+riez2.^2);
        ph = acos(prim./am);
        AM(:,:,l) = am;
        phase(:,:,l)=ph;
        orient(:,:,l)= th;
  
        % Computation of the structure tensor
        a=riez1.^2;
        b=riez1.*riez2;
        c=riez2.^2; 
        J=[mean(a(:)) mean(b(:))
            mean(b(:)) mean(c(:))];
        [V D]=eig(J);
 
        % Structure tensor caracteristics
        Dir(:,:,l)=V;
        Vap(:,:,l)=D;

        coherence(l)=(max(D(1,1),D(2,2))-min(D(1,1),D(2,2)))/(max(D(1,1),D(2,2))+min(D(1,1),D(2,2)));
        Sigma(:,:,l)=cov([reshape(prim(:),dim(1)*dim(2),1),reshape(riez1(:),dim(1)*dim(2),1),reshape(riez2(:),dim(1)*dim(2),1)]);
        Corr(:,:,l)=corrcoef([reshape(prim(:),dim(1)*dim(2),1),reshape(riez1(:),dim(1)*dim(2),1),reshape(riez2(:),dim(1)*dim(2),1)]);
    end
end