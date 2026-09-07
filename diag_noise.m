% diag_noise.m  v4：raw相位 ad 的 直接FFT / pwelch 检测，跨 SNR 是否稳在 25/75
thisdir = fileparts(mfilename('fullpath'));
s = load(fullfile(thisdir,'put_vmd_Mix_40s.mat')); Mix = s.Mix;
fs = 25; Ns = 64;
R_rms = sqrt(mean(abs(Mix(:)).^2));
fprintf('%-6s %-12s %-12s %-12s %-12s\n','SNR','BR_fft_bpm','BR_pw_bpm','HR_fft_bpm','HR_pw_bpm');
for SNR=[inf 40 30 25 20 15 10]
    if isinf(SNR), MN=Mix; lbl='clean'; else rng(SNR); sig=1/sqrt(2)*randn(size(Mix))+1i/sqrt(2)*randn(size(Mix)); MN=Mix+R_rms/10^(SNR/20)*sig; lbl=num2str(SNR); end
    mix1 = squeeze(MN(1,:,1,:));
    RF = fft(mix1 .* hanning(Ns), Ns, 1);
    [~, Mi] = max(mean(abs(RF),2));
    ad = diff(unwrap(angle(RF(Mi,:))));
    M=numel(ad); fr=(-floor(M/2):ceil(M/2)-1)*(fs/M); sp=abs(fftshift(fft(ad)));
    % direct fft
    mR = fr>=0.1 & fr<=0.5; mH = fr>=0.8 & fr<=2;
    [~,iR]=max(sp.*mR); BRf=fr(iR)*60;
    [~,iH]=max(sp.*mH); HRf=fr(iH)*60;
    % pwelch (ad, detrend)
    [Pw,fw]=pwelch(ad(:)-mean(ad),[],[],[],fs);
    mRw=fw>=0.1 & fw<=0.5; mHw=fw>=0.8 & fw<=2;
    [~,iRw]=max(Pw.*mRw); BRp=fw(iRw)*60;
    [~,iHw]=max(Pw.*mHw); HRp=fw(iHw)*60;
    fprintf('%-6s %-12.1f %-12.1f %-12.1f %-12.1f\n', lbl, BRf, BRp, HRf, HRp);
end
