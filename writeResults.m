% Write results to a .txt file (you can import it to a .CSV later)
cd 'MATfiles';
curr_dir = pwd;
F = dir('*.mat');

fid = fopen('../Results.txt', 'w');

% Write the names of variables first
% =====================================
fprintf(fid, '%s,', 'Filename');
fprintf(fid, '%s, %s, %s,', 'NormDuration', 'NumDuration', 'VsVDuration');
fprintf(fid, '%s,', 'RMSamplitude');
fprintf(fid, '%s, %s,', 'RMScons', 'RMSvow');
fprintf(fid, '%s, %s, %s,', 'M1_begin', ...
    'M1_middle', 'M1_end');
fprintf(fid, '%s, %s, %s,', 'M2_begin', ...
    'M2_middle', 'M2_end');
fprintf(fid, '%s, %s, %s,', 'M3_begin', ...
    'M3_middle', 'M3_end');
fprintf(fid, '%s, %s, %s,', 'M4_begin', ...
    'M4_middle', 'M4_end');
fprintf(fid, '%s, %s,', 'DeltaFq', 'DeltaFqPMTM');
fprintf(fid, '%s, %s, %s,', 'MidFreq_begin', ...
    'MidFreq_middle', 'MidFreq_end');
fprintf(fid, '%s, %s, %s,', 'A_lo_min_begin', ...
    'A_lo_min_middle', 'A_lo_min_end');
fprintf(fid, '%s, %s, %s,', 'LDmidhi_begin', ...
    'LDmidhi_middle', 'LDmidhi_end');

fprintf(fid, '%s, %s, %s,', 'A_lo_avg_begin', ...
    'A_lo_avg_middle', 'A_lo_avg_end');
fprintf(fid, '%s, %s, %s,', 'A_midhi_begin', ...
    'A_midhi_middle', 'A_midhi_end');
fprintf(fid, '%s, %s, %s,', 'HighAmp_begin', ...
    'HighAmp_middle', 'HighAmp_end');
fprintf(fid, '%s, %s, %s,', 'HigFreq_begin', ...
    'HighFreq_middle', 'HighFreq_end');
fprintf(fid, '%s, %s, %s,', 'LFAvgAmp_begin', ...
    'LFAvgAmp_middle', 'LFAvgAmp_end');
fprintf(fid, '%s, %s, %s,', 'LFminAmp_begin', ...
    'LFminAmp_middle', 'LFminAmp_end');
fprintf(fid, '%s, %s, %s,', 'LFminFreq_begin', ...
    'LFminFreq_middle', 'LFminFreq_end');
fprintf(fid, '%s, %s, %s,', ...
    'SpecFlattness', 'SpecSlope', 'SpectralPeakLoc');
fprintf(fid, '%s\n', ...
    'FricDur');
    
	
	
	
% Process MAT files
% =====================================
for i = 1:length(F)
    load(F(i).name);
    fprintf(1, 'Processing: %s\n', F(i).name); 
    fprintf(fid, '%s,', F(i).name);
    fprintf(fid, '%f, %f, %f, %f,', NormDur, Consnum, VsVden, RMSamp);
    fprintf(fid, '%f, %f,', 10*log10(RMSc), 10*log10(RMSv));
    fprintf(fid, '%f, %f, %f,', m1MT(1), ...
         m1MT(2),  m1MT(3));
    fprintf(fid, '%f, %f, %f,', m2MT(1), ...
         m2MT(2),  m2MT(3));
    fprintf(fid, '%f, %f, %f,', m3nMT(1), ...
         m3nMT(2),  m3nMT(3));
    fprintf(fid, '%f, %f, %f,', m4nMT(1), ...
         m4nMT(2),  m4nMT(3));
    fprintf(fid, '%f, %f,', DeltaFq, DeltaFqPMTM);
    fprintf(fid, '%f, %f, %f,', MidFreqPMTM(1)*fs/2048, ...
        MidFreqPMTM(2)*fs/2048, MidFreqPMTM(3)*fs/2048);
    fprintf(fid, '%f, %f, %f,', A_lo_minPMTM(1), ...
        A_lo_minPMTM(2), A_lo_minPMTM(3));
    fprintf(fid, '%f, %f, %f,', LD_midhiPMTM(1), ...
        LD_midhiPMTM(2), LD_midhiPMTM(3));

    fprintf(fid, '%f, %f, %f,', A_lo_avgPMTM(1)', ...
         A_lo_avgPMTM(2),  A_lo_avgPMTM(3));
    fprintf(fid, '%f, %f, %f,', A_midhiPMTM(1)', ...
        A_midhiPMTM(2), A_midhiPMTM(3));
    fprintf(fid, '%f, %f, %f,', HighAmpPMTM(1), ...
        HighAmpPMTM(2), HighAmpPMTM(3));
    fprintf(fid, '%f, %f, %f,', HighFreqPMTM(1)*fs/2048, ...
        HighFreqPMTM(2)*fs/2048, HighFreqPMTM(3)*fs/2048);
    fprintf(fid, '%f, %f, %f,', LFAvgAmpPMTM(1), ...
        LFAvgAmpPMTM(2), LFAvgAmpPMTM(3));
    fprintf(fid, '%f, %f, %f,', LFminAmpPMTM(1), ...
        LFminAmpPMTM(2), LFminAmpPMTM(3));
    fprintf(fid, '%f, %f, %f,', LFminFreqPMTM(1)*fs/2048, ...
        LFminFreqPMTM(2)*fs/2048, LFminFreqPMTM(3)*fs/2048);
    fprintf(fid, '%f, %f, %f,', ...
        SFl, SSl, SpectralPeakLocPMTM);
    fprintf(fid, '%f\n', ...
        Consnum);

end       
    

fclose(fid);
cd ..;

return
