 function nima=RemoveStripesVertical(ima,decNum,wname,sigma)
 
 if nargin<2
     decNum = 4;
 end
 
 if nargin<3
     wname = 'db10'; %db10, db3
 end
 
 if nargin<4
     sigma = 6;
 end

 oma = ima;
 ima = imrotate(ima,10);
 ima = single(ima);
 
 % wavelet decomposition
 for ii=1:decNum
    [ima,Ch{ii},Cv{ii},Cd{ii}]=dwt2(ima,wname);
 end

 % FFT transform of horizontal frequency bands
 for ii=1:decNum
     % FFT
     fCv=fftshift(fft(Cv{ii}));
     [my,mx]=size(fCv);

     % damping of vertical stripe information
     damp=1-exp(-(-floor(my/2):-floor(my/2)+my-1).^2/(2*sigma^2));
     fCv=fCv.*repmat(damp',1,mx);

     % inverse FFT
     Cv{ii}=ifft(ifftshift(fCv));
 end

 % wavelet reconstruction
 nima=ima;
 for ii=decNum:-1:1
    nima=nima(1:size(Ch{ii},1),1:size(Ch{ii},2));
    nima=idwt2(nima,Ch{ii},Cv{ii},Cd{ii},wname);
 end
 
 nima = uint16(nima);
nima = imrotate(nima,-10);

 imshow(imadjust(nima))
 return