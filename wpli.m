function PhaseLV = wpli(A, B,tw,Fc)
% This function computes the phase locking values between two 
%signals both in time frequencies( in time-freq)
% The inputs are expected to be multi channels and multi trials series. 
% inputs: the two signals whose plv are to be computed, the signals should
% be in the form of (samples X channel X trials)
% Input: tw is the wavelet time (of form 
% tw=-1/4+1/fs:1/fs:2-1/fs)
% Fc is the wavelet center frequency of form Fc=linspace(6,14,80);

% Output: The time-frequency phase locking values of the two signals
if nargin<4
    error('specify four inputs')
end

if length(A)~=length(B)
    error('A and B must have equal length')
end
p=size(A,1);
numtrials=size(A,3);
m=size(A,2);
n=size(B,2);
% PhaseLV=zeros(m-1,n,4,length(Fc),numtrials);
for Achannel=1:m
    Asig=squeeze(A(:,Achannel,:));    
    
for Bchannel=Achannel+1:n
    Bsig=squeeze(B(:,Bchannel,:));
    
    wav_A=Wavetrans(Asig,tw,Fc);
    wav_B=Wavetrans(Bsig,tw,Fc);
    
    phi_A=angle(wav_A);
    phi_B=angle(wav_B);
    
            
    % Calculate the imaginary part of the cross-spectrum
for ii=1:size(wav_B,2)
    for jj=1:size(wav_A,3)
        imaginary_part(:,ii,jj) = imag(squeeze(wav_A(:,ii,jj)) .* conj(squeeze(wav_B(:,ii,jj))));
        wpli_result(ii,jj) = abs(mean(squeeze(imaginary_part(:,ii,jj)), 'omitnan')) / mean(abs(squeeze(imaginary_part(:,ii,jj))), 'omitnan');
    end
end
    PhaseLV(Achannel,Bchannel,:,:)=wpli_result;
    
end
end

end

