# Temporal, spectral and amplitude characteristics of the Greek fricative /s/ in hearing-impaired and normal-hearing speech

### Anna Sfakianaki, Katerina Nicolaidis, and George P. Kafentzis
---

### MATLAB source code 

The current script runs analyses on the intervocalic fricative /s/. In order to run the code you need a folder similar to the one provided (`SAMPLES`) that includes pairs of .WAV recordings and .TextGrid files that contain Praat-based annotations. This code has been tested on MATLAB 2018a.

### Annotation labels
Fricative /s/: s_c
Vowel: Vx_y, where x=1 when preceding /s/ and x=2 when following /s/, y=1 when stressed and y=2 when unstressed. Please refer to the .TextGrid files in `SAMPLES` folder for more information.

### Additional information
Vowel identity is extracted from folder name. Folders are named after the disyllable and include stress information, e.g. _pisi for [ 'pisi ] and pi_si for [ pi'si ].

A `SAMPLES` folder is provided with sample recordings from the first author.

## HOWTO:

1. Run `readextract.m`: a series of MAT files are constructed for each .WAV file.

2. Run `writeResults.m`: a .TXT file is constructed with the features discussed in the paper for each .WAV file.

3. The .TXT output can be imported to Excel and turned into a CSV for convenience.

## HOWTO for new data:

1. Create a folder similar to the one provided (`SAMPLES`).

2. Include a .WAV file recording and annotations in .TextGrid format.

3. Run the two .m scripts. Make sure you change path names (usually in the first 2-3 lines of each file).
