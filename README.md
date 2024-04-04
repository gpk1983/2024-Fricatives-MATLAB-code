# Temporal, spectral and amplitude characteristics of the Greek fricative /s/ in hearing-impaired and normal-hearing speech

### Anna Sfakianaki, Katerina Nicolaidis, and George P. Kafentzis
---

### MATLAB source code 

The current script runs analyses on the intervocalic fricative /s/. In order to run the code you need a folder similar to the one provided (`SAMPLES`) that includes pairs of .WAV recordings and .TextGrid files that contain Praat-based annotations. 

Annotation labels
Fricative /s/: s_c
Vowel: Vx_y, where x=1 when preceding /s/ and x=2 when following /s/, y=1 when stressed and y=2 when unstressed

Additional information
Vowel identity is extracted from folder name. Folders are named after the disyllable and include stress information, e.g. _pisi for [ 'pisi ] and pi_si for [ pi'si ].

A `SAMPLES` folder is provided with sample recordings from the first author.

## HOWTO:

1. Run `readextract.m`: a series of MAT files are provided for each .WAV file.

2. Run `writeResults.m`: a .TXT file is provided.

3. The output can be imported to Excel and turned into a CSV.
