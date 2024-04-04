% 
% Process dataset and extract information
%system('del *.mat');
clear all;
warning off MATLAB:iofun:UnsupportedEncoding;


curr_dir = pwd;
save_dir = [pwd, '\MATfiles'];

if exist(save_dir, 'dir') ~= 7
    % If the directory does not exist, create it
    mkdir(save_dir);
    disp(['Directory created: ' save_dir]);
else
    disp(['Directory already exists: ' save_dir]);
end

% The folder with .WAV and .TextGrid files
cd 'SAMPLES';
F = dir;

cutoff_flag = 1;
normalize_flag = 1;
plot_flag = 0;

% Process files
for folder = 3:length(F) % skip . and ..
    cd([pwd '\' F(folder).name]);
    A = dir('*.TextGrid');
    dline = " ";
    for i = 1:length(A)
        %fid = fopen(A(i).name, 'r','n', 'Unicode');
        fid = fopen(A(i).name, 'r');
        fprintf('Processing %s\n', A(i).name);
        %txtData = textscan(fid,'%[^\n\r]');
        %txtData = txtData{:};
        prev3line = fgetl(fid);
        prev2line = fgetl(fid);
        prev1line = fgetl(fid);
        dline = fgetl(fid);
        % Get start and end of consonant /s/ and succeeding vowel
        while ~contains(dline, 'text = "s_c') 
            prev3line = prev2line;
            prev2line = prev1line;
            prev1line = dline;
            dline = fgetl(fid);
            stress = strfind(dline, '_s');
            if isempty(stress)
                firstVowelstress = 1;
            else
                firstVowelstress = 0;
            end
        end
        % Get start and end of consonant
        cons_start = str2double(prev2line(20:end));
        cons_end = str2double(prev1line(20:end));
        % Mark interval number
        intervalNo = str2double(prev3line(20));

        % Move on to the suceeding vowel
        prev2line = fgetl(fid);
        prev1line = fgetl(fid);
        dline = fgetl(fid);
        while ~contains(dline, 'text = "e2') && ~contains(dline, 'text = "a2') ...
                && ~contains(dline, 'text = "o2') && ~contains(dline, 'text = "i2') ...
                && ~contains(dline, 'text = "e2') && ~contains(dline, 'text = "''a2') ...
                && ~contains(dline, 'text = "i2') && ~contains(dline, 'text = "''i2') ...
                && ~contains(dline, 'text = "u2') && ~contains(dline, 'text = "''u2')

            prev2line = prev1line;
            prev1line = dline;
            dline = fgetl(fid);
            if feof(fid)
                disp('No second vowel identified! Skipping');
                break;
            end
        end
        % Get start and end of vowel
        vowel_start = str2double(prev2line(20:end));
        vowel_end = str2double(prev1line(20:end));
        idxs = strfind(dline, """");
        vowel_ID = dline(idxs(1):idxs(2));

        % Get previous vowel (before stop consonant)
        frewind(fid);
        while ~contains(dline, ['intervals [' num2str(intervalNo-1) ']'])
            dline = fgetl(fid);
            continue;
        end
        dline = fgetl(fid);
        prevVow_start = str2double(dline(20:end));
        dline = fgetl(fid);
        prevVow_end = str2double(dline(20:end));
        dline = fgetl(fid);
        idxs = strfind(dline, """");
        prevVow_ID = dline(idxs(1):idxs(2));


        % Cut stuff
		% =====================================
        % Define sample values
        [s,fs] = audioread([A(i).name(1:end-9) '.wav']);
        prevVow_start_sample = round(prevVow_start*fs);
        prevVow_end_sample = round(prevVow_end*fs);
        cons_start_sample = round(cons_start*fs);
        cons_end_sample = round(cons_end*fs);
        vowel_start_sample = round(vowel_start*fs);
        vowel_end_sample = round(vowel_end*fs);
		
		% ===========================================
		%
		%   Start Computing Features
		% 
		% ==========================================
		
		
		% =====================================
        % (a) Normalized duration = Dur(s)/Dur(VsV)
        % ========================================
        Consnum = cons_end_sample - cons_start_sample + 1;
        VsVden = vowel_end_sample - prevVow_start_sample + 1;
        NormDur = Consnum/VsVden;
		% ========================================
		
		
        % =====================================
        % (b) Spectral peak location
        % ==========================
        winsize_ms = 40;
        winsize_smp = round(winsize_ms*10^(-3)*fs);
        % sanity check in case window is larger than fricative
        if winsize_smp > (cons_end_sample - cons_start_sample + 1)
            disp('Error in spectral peak: window size is larger than fricative size!');
            disp('Reducing window size to 25 ms');
            winsize_ms = 25;
            winsize_smp = round(winsize_ms*10^(-3)*fs);
        end
        % Create a normalized 40 ms Hamming window
        % and place it in the middle of the fricative
        if mod(winsize_smp,2) == 0
            winsize_smp = winsize_smp + 1;
            win = hamming(winsize_smp);
        else
            win = hamming(winsize_smp);
        end
        win = win./sum(win);
        cons = s(cons_start_sample:cons_end_sample);
        midpoint = round(length(cons)/2);
        winframe = cons(midpoint - (length(win)-1)/2:midpoint + (length(win)-1)/2).*win;

        % Apply FFT of NFFT points
        NFFT = 2048;
        FTmag = abs(fft(winframe, NFFT));
        FTh = 2*FTmag(1:NFFT/2+1);
        PowerSp = FTh.^2;
        PowerSpPMTM = pmtm(winframe, 4, NFFT);
        if normalize_flag
            PowerSp = PowerSp/sum(PowerSp);
            PowerSpPMTM = PowerSpPMTM./sum(PowerSpPMTM);
        end
        Fv = 0:fs/NFFT:fs/2;
        [~, ploc] = max(PowerSp);
        SpectralPeakLoc = ploc*fs/NFFT;
        [~, ploc] = max(PowerSpPMTM);
        SpectralPeakLocPMTM = ploc*fs/NFFT;
		% ========================================
		
		
		
        % =============================================
        % (c1) Spectral flatness and (c2) spectal slope
        % =============================================
        % (c1)
        PowerSp = PowerSpPMTM;
        SFl = exp(mean(log(PowerSp))) ./ (mean(PowerSp));
        % (c2)
        % compute means and terms
        term1 = Fv - mean(Fv);
        term2 = 10*log10(PowerSp) - mean(10*log10(PowerSp));
        numer = sum(term1(:).*term2(:));
        denom = sum(term1(:).^2);
        SSl = numer./denom;
		% ========================================
		
		
		
        % ======================================================
        % (d) Spectral Moments
        % ======================================================
        % conventional (using single FFT)
        % and multitaper-based method
        % ======================================================
        mt_win_ms = 25; % in ms
        mt_win_smp = round(mt_win_ms*10^(-3)*fs);
        if mod(mt_win_smp,2) == 0
            mt_win_smp = mt_win_smp + 1;
            mtwin = hamming(mt_win_smp);
        else
            mtwin = hamming(mt_win_smp);
        end
        mtwin = mtwin./sum(mtwin);
        Fv = 0:fs/NFFT:fs/2;
        
        % Beginning
        cons = s(cons_start_sample:cons_end_sample);
        winframe_begin = cons(1:length(mtwin)).*mtwin;
        FTmag_begin = abs(fft(winframe_begin, NFFT));
        FTh = 2*FTmag_begin(1:NFFT/2+1);
        Power_begin = FTh.^2;
		
        % MT method
        PowerPMTM(1,:) = pmtm(winframe_begin, 4, NFFT);
        if normalize_flag
            Power(1,:) = Power_begin/sum(Power_begin);
            PowerPMTM(1,:) = PowerPMTM(1,:)/sum(PowerPMTM(1,:));
        else
            Power(1,:) = Power_begin;
            PowerPMTM(1,:) = PowerPMTM(1,:);
        end
        
        
        % Middle
        cons = s(cons_start_sample:cons_end_sample);
        midpoint = round(length(cons)/2);
        winframe_mid = cons(midpoint - (length(mtwin)-1)/2:midpoint + (length(mtwin)-1)/2).*mtwin;
        FTmag_mid = abs(fft(winframe_mid, NFFT));
        FTh = 2*FTmag_mid(1:NFFT/2+1);
        Power_mid = FTh.^2;
		
        % MT method
        PowerPMTM(2,:) = pmtm(winframe_mid, 4, NFFT);
        if normalize_flag
            Power(2,:) = Power_mid/sum(Power_mid);
            PowerPMTM(2,:) = PowerPMTM(2,:)/sum(PowerPMTM(2,:));
        else
            Power(2,:) = Power_mid;
            PowerPMTM(2,:) = PowerPMTM(2,:);
        end
        
        
        % End
        cons = s(cons_start_sample:cons_end_sample);
        winframe_end = cons(end-length(mtwin)+1:end).*mtwin;
        FTmag_end = abs(fft(winframe_end, NFFT));
        FTh = 2*FTmag_end(1:NFFT/2+1);
        Power_end = FTh.^2;
		
        % MT method
        PowerPMTM(3,:) = pmtm(winframe_end, 4, NFFT);
        if normalize_flag
            Power(3,:) = Power_end/sum(Power_end);
            PowerPMTM(3,:) = PowerPMTM(3,:)/sum(PowerPMTM(3,:));
        else
            Power(3,:) = Power_end;
            PowerPMTM(3,:) = PowerPMTM(3,:);
        end

        if 1 % cutoff_flag always on
            % Threshold @ Fc
            Fc = 550;
            Fcbin = round(NFFT*Fc/fs);
			if normalize_flag
                for z = 1:3
				    PowerNorm(z, :) = Power(z, Fcbin:end)./sum(Power(z, Fcbin:end));
				    PowerPMTMNorm(z, :) = PowerPMTM(z, Fcbin:end)./sum(PowerPMTM(z, Fcbin:end));
                end
%             else %normally does not occur but it should!!!
% 				Power = Power(:, Fcbin:end);
% 				PowerPMTM = PowerPMTM(:, Fcbin:end);
            end
            Power = Power(:, Fcbin:end);
 			PowerPMTM = PowerPMTM(:, Fcbin:end);
            Fv = Fv(Fcbin:end);
        end
        
        if plot_flag
            subplot(2,1,1);
            plot(Fv, 10*log10(PowerPMTM)); title('Multitaper');
            legend('begin', 'mid', 'end');
            subplot(2,1,2); 
            plot(Fv, 10*log10(Power)); title('Single-frame');
            legend('begin', 'mid', 'end');
        end
        
        % Calculate moments - each line contains [beg, mid, end]
        for partOfFric=1:3
			% Standard FFT Power Spectrum Moments
            m1(partOfFric) = Fv*PowerNorm(partOfFric, :).';
            m2(partOfFric) = ((Fv-m1(partOfFric)).^2)*PowerNorm(partOfFric, :)';
            m3(partOfFric) = ((Fv-m1(partOfFric)).^3)*PowerNorm(partOfFric, :)';
            m4(partOfFric) = ((Fv-m1(partOfFric)).^4)*PowerNorm(partOfFric, :)';
            m3n(partOfFric) = m3(partOfFric)./((m2(partOfFric)).^(3/2));
            m4n(partOfFric) = m4(partOfFric)./(m2(partOfFric).^2) - 3; %normal distribution has kurtosis 0
            m2(partOfFric) = sqrt(m2(partOfFric)); % standard deviation, same units with mean
            
            % Multitaper Power Spectrum Moments
            m1MT(partOfFric) = Fv*PowerPMTMNorm(partOfFric,:).';
            m2MT(partOfFric) = ((Fv-m1MT(partOfFric)).^2)*PowerPMTMNorm(partOfFric, :).';
            m3MT(partOfFric) = ((Fv-m1MT(partOfFric)).^3)*PowerPMTMNorm(partOfFric, :).';
            m4MT(partOfFric) = ((Fv-m1MT(partOfFric)).^4)*PowerPMTMNorm(partOfFric, :).';
            m3nMT(partOfFric) = m3MT(partOfFric)./((m2MT(partOfFric)).^(3/2));
            m4nMT(partOfFric) = m4MT(partOfFric)./(m2MT(partOfFric).^2) - 3; %normal distribution has kurtosis 0
            m2MT(partOfFric) = sqrt(m2MT(partOfFric)); % standard deviation, same units with mean
        end
        clear PowerNorm
        clear PowerPMTMNorm
		% =============================================================================
		
		
        
        % Frequency bands (in bins)
        % LF -> 550 - 3000 Hz
        % MF -> 3000 - 7000 Hz
        % HF -> 7000 - 11025 Hz (fs/2)
        % ======================================================
        LF = round(NFFT*550/fs) - Fcbin + 1:round(NFFT*3000/fs) - Fcbin + 1;
        MF = LF(end)+1:round(NFFT*7000/fs) - Fcbin + 1;
        HF = MF(end)+1:NFFT/2 - Fcbin + 1;
        Df = fs/NFFT;
        
		
		% ======================================
        % Other parameters
		% =====================================
        for partOfFric=1:3
            % Peaks and amps
            % =====================================
            % (f) MidFreq, MidAmp - (Hz) & (dB)
			% =====================================
            [MidAmp(partOfFric), MidFreq(partOfFric)] = max(10*log10(Power(partOfFric, MF)));
            [MidAmpPMTM(partOfFric), MidFreqPMTM(partOfFric)] = max(10*log10(PowerPMTM(partOfFric, MF)));
            MidFreqPMTM(partOfFric) = MidFreqPMTM(partOfFric) + MF(1) - 1 + Fcbin - 1;
            MidFreq(partOfFric) = MidFreq(partOfFric) + MF(1) - 1 + Fcbin - 1;
			
			
			% =====================================
            % (g) LFmin, LFfeq - (Hz) & (dB)
			% =====================================
            [LFminAmp(partOfFric), LFminFreq(partOfFric)] = min(10*log10(Power(partOfFric, LF)));
            [LFminAmpPMTM(partOfFric), LFminFreqPMTM(partOfFric)] = min(10*log10(PowerPMTM(partOfFric, LF)));
            LFminFreqPMTM(partOfFric) = LFminFreqPMTM(partOfFric) + Fcbin - 1;
            LFminFreq(partOfFric) = LFminFreq(partOfFric) + Fcbin - 1;



			% ===========================================
            % (h) LFAmp (dB) (average amplitude in LF)
			% ===========================================
            LFAvgAmp(partOfFric) = mean(10*log10(Power(partOfFric, LF)));
            LFAvgAmpPMTM(partOfFric) = mean(10*log10(PowerPMTM(partOfFric, LF)));
            
			
			
			% =====================================
            % (i) HFfreq, HFAmp - (Hz) & (dB)
			% =====================================
            [HighAmp(partOfFric), HighFreq(partOfFric)] = max(10*log10(Power(partOfFric, HF)));
            [HighAmpPMTM(partOfFric), HighFreqPMTM(partOfFric)] = max(10*log10(PowerPMTM(partOfFric, HF)));
            HighFreq(partOfFric) = HighFreq(partOfFric) + HF(1) - 1 + Fcbin - 1;
            HighFreqPMTM(partOfFric) = HighFreqPMTM(partOfFric) + HF(1) - 1 + Fcbin - 1;



			% =====================================
            % Amplitude parameters
            % =====================================
            % (j) A_lo_avg
			% =====================================
            A_lo_avg(partOfFric) = MidAmp(partOfFric) - LFAvgAmp(partOfFric);
            A_lo_avgPMTM(partOfFric) = MidAmpPMTM(partOfFric) - LFAvgAmpPMTM(partOfFric);
            
			
			
			% =====================================
            % (k) A_lo_min
			% =====================================
            A_lo_min(partOfFric) = MidAmp(partOfFric) - LFminAmp(partOfFric);
            A_lo_minPMTM(partOfFric) = MidAmpPMTM(partOfFric) - LFminAmpPMTM(partOfFric);
            
			
			
			% =====================================
            % (l) A_midhi
			% =====================================
            A_midhi(partOfFric) = MidAmp(partOfFric) - HighAmp(partOfFric);
            A_midhiPMTM(partOfFric) = MidAmpPMTM(partOfFric) - HighAmpPMTM(partOfFric);
            
			
			
			% =====================================
            % (m) LD_midhi (dB)
			% =====================================
            LD_midhi(partOfFric) = Df*sum((Power(partOfFric, MF))) ...
                - Df*sum((Power(partOfFric, HF)));
            LD_midhiPMTM(partOfFric) = Df*sum((PowerPMTM(partOfFric, MF))) ...
                - Df*sum((PowerPMTM(partOfFric, HF)));
        end
        
		
		% =====================================
        % DeltaFq (Hz)
		% =====================================
        DeltaFq = mean(MidFreq(1:2)*fs/NFFT) - MidFreq(3)*fs/NFFT;
        DeltaFqPMTM = mean(MidFreqPMTM(1:2)*fs/NFFT) - MidFreqPMTM(3)*fs/NFFT;
        
		
		% ======================================================
        % (e) RMS amplitude (dB) (referenced on stressed vowel)
        % ======================================================
        RMSc = rms(cons);
        if firstVowelstress 
            Vow = s(prevVow_start_sample:prevVow_end_sample);
        else
            Vow = s(vowel_start_sample:vowel_end_sample);
        end
        [maxVowAmp, maxVowT] = max(Vow);
        xc = xcorr(Vow);
        xc = xc((length(xc)+1)/2:end);
        [amax, tmax] = findpeaks(xc);
        [ampmax, idxmax] = max(amax);
        % pitch period
        N0 = tmax(idxmax);
        halfWin = round(3*N0/2);
        % window around vowel max amplitude
        threePitchPidx = maxVowT-halfWin:maxVowT + halfWin;
        if sum(threePitchPidx <= 0)
            threePitchPidx = maxVowT-halfWin+sum(threePitchPidx < 0)+1:maxVowT + halfWin;
        end
        if threePitchPidx(end) > length(Vow)
            threePitchPidx = maxVowT-3*halfWin:maxVowT;
        end
        RMSv = rms(Vow(threePitchPidx));
        % RMS amplitude in dB
        RMSamp = 10*log10(RMSc/RMSv);
        
		
		
		
        % =====================================
        % Save stuff
		% =====================================
        savefile = strcat(A(i).name(1:end-9), '_data.mat');
        matname = fullfile(save_dir, savefile);
        save(matname, 'savefile', 's', 'fs', 'NormDur', 'Consnum', 'VsVden', 'NFFT',...
            'SpectralPeakLocPMTM', 'SFl', 'SSl', ...
            'm1MT', 'm2MT', 'm3nMT', 'm4nMT', ...
            'DeltaFq', 'DeltaFqPMTM', 'RMSamp', 'RMSc', 'RMSv', ...
            'MidFreqPMTM', 'LFminAmpPMTM', 'LFminFreqPMTM',...
            'LFAvgAmpPMTM', 'HighAmpPMTM', 'HighFreqPMTM', ...
            'A_lo_avgPMTM', 'A_lo_minPMTM', 'A_midhiPMTM', ...
            'LD_midhiPMTM', 'PowerPMTM');
        
        % Close the file and move on
        fclose(fid);
        % Clear variables
        clear Power PowerPMTM
    end    
    cd('..');
end