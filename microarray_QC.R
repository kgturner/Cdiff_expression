#arrays with bubbles/dust/artifacts/identified spatial biases/outlier qc values/etc
#.xys files and experimental metrics reports generated with NimbleScan2.6
#artifacts trimmed using /home/morien/bin/sunflower_microarray/bad_features.kay_microarray.pl and remove_bad_features.sh

#sunflower expt QC

#large artifacts and weird metrics. probably no hope for these.
471145A03_532.tif	HDA JAS1 2	#top 1/4 of array ruined, throw out	#huge (dark) artifact/bubble? upper right, bad washing?	v high Uniformity CV (spatial bias in mean intensity)	v high Signal Range(indicative of spatial intensity bias)
471187A07_532.tif	IND JAS1 1	#throw out	#large bubble center left.	v high Signal Range(indicative of spatial intensity bias)	v high mean random (higher than normal non-specific binding).	v high Uniformity CV (spatial bias in mean intensity)
471312A03_532.tif	CRB JAS2 2	#throw out	#huge artifact below array. light affecting lower half of array.	v high Signal Range(indicative of spatial intensity bias)	v high Mean Random (higher than normal non-specific binding).	higher than normal Mean Empty (possible high background, improper washing, bad sample)

#weird metrics, not sure whether to throw out
471187A12_532.tif	HA89 DR2 1	#unsure whether to throw out 'rainbow' view looks ok i guess	#large v bright dust left of array. would this affect the image intensity? the image is darker than most.
471198A06_532.tif	HA89 DR2 3	#weird spot upper left, not sure if we should discard 	#v dark array with bright control spot grid pattern. low sample yield? might want to throw out.
470975A01_532.tif		#large washing artifact. metrics are fine. going to trim ~1/6 of the array. at this point, is it a throwaway?
102,56 259,256
#very dark corners, not sure whether to throw out. could save?
471186A11_532.tif	IND DR1 1	#upper left corner v dark
4,4 45,38
471145A02_532.tif	HDA JAS1 1	#dark spot lower left, dust lower left, lower right corner v dark
102,432 117,458; 43,496 65,514; 415,590 458,633

#smaller artifacts. don't affect array quality metrics. fix these with perl script.
470975A02_532.tif	HDA CTL 1	#dust center right
433,316 452,333
471139A09_532.tif	CMT CTL 3	#smudge (tiny bubble?) lower right quadrant
268,405 292,428
471148A06_532.tif	HA89 DR1 1	#dust upper right quadrant
267,188 284,198
471187A10_532.tif	GRE CTL 2	#small bubble lower left quadrant
25,481 76,496
471198A04_532.tif	KAN DR2 1	#small bubble upper left quadrant
92,57 97,63
471198A12_532.tif	COL DR2 3	#large bright dust upper left quadrant
158,67 166,92
471225A11_532.tif	HA89 JAS2 1	#dust lower right quadrant
385,429 405,454
471312A12_532.tif	UTA JAS1 3	#dust upper left quadrant
176,51 199,75

#weird metrics. would like to normalize and compare to the rest of the data before throwing them out.
471312A07_532.tif	IND JAS1 3	high Mean Empty (indicative of bad washing, bad sample)
471208A02_532.tif	KAN DR2 3	high Mean Empty (indicative of bad washing, bad sample)
471312A06_532.tif	CRB JAS1 2	high Mean Empty (indicative of bad washing, bad sample)


#in folder: /media/Childs/microarray/sunflower/
#"sunflower" expt pre-processing
a=read.delim("tif_aligned/experimental_metrics_report.txt", header=T, row.names=1, sep="\t")

#plot mean of 'empty' spots. good background measure
plot(sort(a[,10])) #a few v high. looks like 810 is the cutoff for the main group. above that there are a few outliers. 
#mean across all arrays is much higher than 'recommended' by nimblegen too (they say ~400 is normal. our mean is 600.61.)

> row.names(a)[which(a[,10] > 810)] 
#potential problem arrays/samples
[1] IND JAS1 3	894.422	471312A07_532.tif
[1] CRB JAS2 2	1211.461	471312A03_532.tif
[1] GRE JAS2 2	820.044	471312A01_532.tif
[1] CRB JAS1 2	961.321	471312A06_532.tif
[1] COL JAS2 1	848.016	471208A06_532.tif
[1] IND JAS1 1	815.578	471187A07_532.tif
[1] KAN DR2 3	927.899	471208A02_532.tif


#NIL expt QC
#THROW OUT
309556A01_532.tif	#washing artifacts lower half. v dark array. signal intensity gradient from L to R.	#horrible metrics	#THROW OUT

#undecided
309510A04_X4_532.tif	#spot near left edge	#high signal range	#high uniformity CV	#CONSIDER THROWING OUT
3,189 5,193
309510A03_X4_532.tif	#upper right & center darker than rest of array	#looks fine in nimblescan	#also high signal range	#also high uniformity CV	#CONSIDER THROWING OUT

#keep
309503A03_X4_532.tif	#dust top right
282,1 288,6
309503A04_X4_532.tif	#weird spot lower half
177,294 191,307
309510A02_X4_532.tif	#dust lower right quadrant, upper right quadrant
238,290 246,296; 277,55 285,64
309510A01_X4_532.tif	#dust upper right quadrant
211,132 214,140
#samples below look 'okay' in nimblescan view. will check metrics.	#metrics ok
309506A01_X4_532.tif	#dark array. low sample yield?
309506A02_X4_532.tif	#dark array. low sample yield?
309506A03_X4_532.tif	#dark array. low sample yield?
309506A04_X4_532.tif	#dark array. low sample yield?

#METRICS
#high_signal_range
230047A01_X4_532.tif	BC4S1_16:8_D+4_REP3_DW
230047A02_X4_532.tif	BC4S1_16:8_D+4_REP3_WW
309510A03_X4_532.tif	8:16_Mstem_25day_A3
309510A04_X4_532.tif	SD-LD_Mstem_25day_A1
309556A01_532.tif	8:16_Mstem_25day_B1

#high_uniformity_CV
309510A03_X4_532.tif	8:16_Mstem_25day_A3
309510A04_X4_532.tif	SD-LD_Mstem_25day_A1
309556A01_532.tif	8:16_Mstem_25day_B1

#v_high_mean_empty
309512A01_X4_532.tif	SD-LD_Mstem_25day_A2
309512A04_X4_532.tif	8:16_Mstem_25day_H2
309556A01_532.tif	8:16_Mstem_25day_B1	#throw away for sure
309556A02_532.tif	8:16_Mstem_25day_B3	#otherwise fine

#v_high_mean_experimental
309556A01_532.tif	8:16_Mstem_25day_B1	#doesn't look good. throw away
309556A02_532.tif	8:16_Mstem_25day_B3	#otherwise fine



#AMF expt QC

#definitely bad
309679A02_X4_Hannuus_532.tif	#large bubble upper left, dark line down center right, smudge center right. THROW OUT
226261A02_X4_U01_Heliantus_2_532_2.tif	#large washing artifact top of array. DO NOT USE.	#high mean empty	#high uniformity CV	#high signal range
226261A02_X4_U01_Heliantus_532_2.tif	#v dark array with R to L gradient.large washing artifact top of array. DO NOT USE. 	#high uniformity CV 	#high signal range
226261A02_X4_U01_Heliantus_2_532.tif	#large washing artifact top of array. DO NOT USE.	#high mean empty	#high uniformity CV	#high signal range
226261A02_X4_U01_Heliantus_532.tif	#v dark array with R to L gradient. large washing artifact top of array. DO NOT USE. 	#high uniformity CV 	#high signal range

#maybe bad
309514A01_X4_Hannuus_532.tif	#dark array. R to L gradient. LOOKS WEIRD. Background is too high. normalize but consider throwing out.
309514A02_X4_Hannuus_532.tif	#dust lower left quadrant. LOOKS WEIRD. Background is too high. normalize but consider throwing out.
123,352 128,361
309514A03_X4_Hannuus_532.tif	#LOOKS WEIRD. Background is too high. normalize but consider throwing out.
309514A04_X4_Hannuus_532.tif	#dust upper left quadrant. LOOKS WEIRD. Background is too high. normalize but consider throwing out.
91,168 94,176
309679A01_X4_Hannuus_532.tif	#large bubble upper right, bad spot on array, upper right, dust far right	#high signal range
170,1 207,32; 195,71 233,142; 146,129 174,155; 329,85 333,88
309706A03_X4_Hannuus_532.tif	#dark spot lower left	#high signal range
120,420 160,457
226261A02_X4_U01_Heliantus_3_532_2.tif	#v dark array with R to L gradient	#high uniformity CV	#high signal range
101,25 191,71
226261A02_X4_U01_Heliantus_3_532.tif	#v dark array with R to L gradient. better scan of the same image. can trim out bad spot.
101,25 191,71

#dust/scratches/bubbles but otherwise fine
225583A03_X4_U01_Heliantus_532.tif	#dust speck upper right corner
50,61 53,67
309518A04_X4_Hannuus_532.tif	#dust center right, specks under and to the left of the large piece
190,236 221,249; 204,319 207,321; 99,300 102,304
309654A02_X4_U01_Helianthus_532.tif	#large washing artifact right side	#metrics look fine
313,128 333,261
309659A03_X4_U01_Helianthus_532.tif	#dark spot below center
191,324 196,328
309670A03_X4_U01_Helianthus_532.tif	#small bubble center right
286,227 291,231
309679A03_X4_Hannuus_532.tif	#dark spot lower right
225,375 255,397
309709A04_X4_U01_Helianthus_532.tif	#dark spot lower right
177,440 201,459
309709A01_X4_U01_Helianthus_532.tif	#dark spot top center
167,3 201,25
309706A04_X4_Hannuus_532.tif	#dark spot lower right
178,426 202,452
309706A01_X4_Hannuus_532.tif	#washing artifact top center
170,3 218,37
309717A04_X4_U01_Helianthus_532.tif	#dark area lower left
66,425 71,427

#19933_20090408 and 19933_20090409 are dupes of the same scans, not just the same arrays. metrics are slightly off between sets, so they are new scans of the same array sets. 

#METRICS

#HIGH MEAN EMPTY
309514A01_X4_Hannuus_532.tif         309514A01_X4_Hannuus_532	STER CO WILD H106
309514A02_X4_Hannuus_532.tif         309514A02_X4_Hannuus_532	STER CO WILD H107
309514A03_X4_Hannuus_532.tif         309514A03_X4_Hannuus_532	STER CAW WEEDY H73
309514A04_X4_Hannuus_532.tif         309514A04_X4_Hannuus_532	STER UTW WEEDY H78
309679A02_X4_Hannuus_532.tif         309679A02_X4_Hannuus_532	STER UTW WEEDY H79
226261A02_X4_U01_Heliantus_2_532_2.tif 226261A02_X4_U01_Heliantus_2_532	STER CUL DOMESTICATED H13
226261A02_X4_U01_Heliantus_2_532.tif 226261A02_X4_U01_Heliantus_2_532	STER CUL DOMESTICATED H13

#HIGH MEAN EXPERIMENTAL	#none of these are extremely high, but they sit in a group separated from the rest of the arrays	#normalization may be ok for these
225575A04_X4_U01_Heliantus_532.tif 225575A04_X4_U01_Heliantus_532
225574A01_X4_U01_Heliantus_532_2.tif 225574A01_X4_U01_Heliantus_532
225574A02_X4_U01_Heliantus_532_2.tif 225574A02_X4_U01_Heliantus_532
225574A03_X4_U01_Heliantus_532_2.tif 225574A03_X4_U01_Heliantus_532
225574A04_X4_U01_Heliantus_532_2.tif 225574A04_X4_U01_Heliantus_532
225574A01_X4_U01_Heliantus_532.tif 225574A01_X4_U01_Heliantus_532
225574A02_X4_U01_Heliantus_532.tif 225574A02_X4_U01_Heliantus_532
225574A03_X4_U01_Heliantus_532.tif 225574A03_X4_U01_Heliantus_532
225574A04_X4_U01_Heliantus_532.tif 225574A04_X4_U01_Heliantus_532

#HIGH UNIFORMITY CV
225575A01_X4_U01_Heliantus_532.tif   225575A01_X4_U01_Heliantus_532
225575A02_X4_U01_Heliantus_532.tif   225575A02_X4_U01_Heliantus_532
309514A04_X4_Hannuus_532.tif         309514A04_X4_Hannuus_532
226261A01_X4_U01_Heliantus_532_2.tif   226261A01_X4_U01_Heliantus_532
226261A02_X4_U01_Heliantus_2_532_2.tif 226261A02_X4_U01_Heliantus_2_532
226261A02_X4_U01_Heliantus_3_532_2.tif 226261A02_X4_U01_Heliantus_3_532
226261A02_X4_U01_Heliantus_532_2.tif   226261A02_X4_U01_Heliantus_532
226261A03_X4_U01_Heliantus_532_2.tif   226261A03_X4_U01_Heliantus_532
226261A01_X4_U01_Heliantus_532.tif   226261A01_X4_U01_Heliantus_532
226261A02_X4_U01_Heliantus_2_532.tif 226261A02_X4_U01_Heliantus_2_532
226261A02_X4_U01_Heliantus_3_532.tif 226261A02_X4_U01_Heliantus_3_532
226261A02_X4_U01_Heliantus_532.tif   226261A02_X4_U01_Heliantus_532
226261A03_X4_U01_Heliantus_532.tif   226261A03_X4_U01_Heliantus_532

#HIGH UNIFORMITY MEAN	#none of these are extremely high, but they sit in a group separated from the rest of the arrays	#normalization may be ok for these
225575A04_X4_U01_Heliantus_532.tif 225575A04_X4_U01_Heliantus_532
225574A01_X4_U01_Heliantus_532_2.tif 225574A01_X4_U01_Heliantus_532
225574A02_X4_U01_Heliantus_532_2.tif 225574A02_X4_U01_Heliantus_532
225574A03_X4_U01_Heliantus_532_2.tif 225574A03_X4_U01_Heliantus_532
225574A04_X4_U01_Heliantus_532_2.tif 225574A04_X4_U01_Heliantus_532
225574A01_X4_U01_Heliantus_532.tif 225574A01_X4_U01_Heliantus_532
225574A02_X4_U01_Heliantus_532.tif 225574A02_X4_U01_Heliantus_532
225574A03_X4_U01_Heliantus_532.tif 225574A03_X4_U01_Heliantus_532
225574A04_X4_U01_Heliantus_532.tif 225574A04_X4_U01_Heliantus_532

#HIGH SIGNAL RANGE	#high uniformity cv club plus extras
225575A01_X4_U01_Heliantus_532.tif   225575A01_X4_U01_Heliantus_532
225575A02_X4_U01_Heliantus_532.tif   225575A02_X4_U01_Heliantus_532
309679A01_X4_Hannuus_532.tif         309679A01_X4_Hannuus_532
309679A02_X4_Hannuus_532.tif         309679A02_X4_Hannuus_532
309706A03_X4_Hannuus_532.tif         309706A03_X4_Hannuus_532
226261A01_X4_U01_Heliantus_532_2.tif   226261A01_X4_U01_Heliantus_532
226261A02_X4_U01_Heliantus_2_532_2.tif 226261A02_X4_U01_Heliantus_2_532
226261A02_X4_U01_Heliantus_3_532_2.tif 226261A02_X4_U01_Heliantus_3_532
226261A02_X4_U01_Heliantus_532_2.tif   226261A02_X4_U01_Heliantus_532
226261A01_X4_U01_Heliantus_532.tif   226261A01_X4_U01_Heliantus_532
226261A02_X4_U01_Heliantus_2_532.tif 226261A02_X4_U01_Heliantus_2_532
226261A02_X4_U01_Heliantus_3_532.tif 226261A02_X4_U01_Heliantus_3_532
226261A02_X4_U01_Heliantus_532.tif   226261A02_X4_U01_Heliantus_532

#outlier interquartile density
309514A01_X4_Hannuus_532.tif 309514A01_X4_Hannuus_532 STER CO WILD H106


#star thistle expt QC

#DEFINITELY THROW OUT
#220768A04_X4_U01_Centaurea_532.tif	#high background, looks pretty dirty...	dark spot lower left	#consider THROWING OUT	#V HIGH BACKGROUND (MEAN EMPTY)	#EXTREMELY HIGH MEAN RANDOM
78 364 101 381	
#266698A02_X4_U01_Centaurea_532.tif	#huge scratch down the center. THROW OUT.	#high signal range	#high uniformity CV
#299636A02_X4_U01_Centaurea_532.tif	#huge washing artifact	#THROW OUT	#high signal range
#299709A02_X4_U01_Centaurea_532.tif	#extremely low signal	#throw out	#very low uniformity mean	#VERY LOW EXPERIMENTAL MEAN

#QUESTIONABLE
220768A02_X4_U01_Centaurea_532.tif	#dark spot top center, weird streaking	#consider THROWING OUT
114 3 148 18
220768A03_X4_U01_Centaurea_532.tif	#weird streaking	#consider THROWING OUT
260942A01_X4_U01_Centaurea_532.tif	#upper right quadrant brighter	#dark spot top center	#V HIGH BACKGROUND (MEAN EMPTY)	#EXTREMELY HIGH MEAN RANDOM
167 3 184 11
260942A02_X4_U01_Centaurea_532.tif	#badly washed top left. nothing to remove though	#V HIGH BACKGROUND (MEAN EMPTY)
260942A03_X4_U01_Centaurea_532.tif	#badly washed bottom center
260942A04_X4_U01_Centaurea_532.tif	#dark spots bottom center	#V HIGH BACKGROUND (MEAN EMPTY)	#EXTREMELY HIGH MEAN RANDOM
217 432 232 450; 169 424 189 442
261033A01_X4_U01_Centaurea_532.tif	#very high background	#dark spot top center	#EXTREMELY HIGH MEAN RANDOM
188 3 210 18
261033A04_X4_bis_U01_Centaurea_532.tif	#gradient. bright bottom right
263785A03_X4_U01_Centaurea_532.tif	#low signal?
265944A02_X4_U01_Centaurea_2_532.tif	#EXTREMELY BRIGHT ARRAY	#dust bottom center	#EXTREMELY HIGH MEAN RANDOM
213 445 218 461
265944A03_X4_U01_Centaurea_2_532.tif	#EXTREMELY BRIGHT ARRAY	#high uniformity mean	#EXTREMELY HIGH MEAN RANDOM	#HIGH MEAN EXPERIMENTAL
265944A04_X4_U01_Centaurea_532.tif	#EXTREMELY BRIGHT ARRAY	#high uniformity mean	#EXTREMELY HIGH MEAN RANDOM	#HIGH MEAN EXPERIMENTAL
265945A01_X4_U01_Centaurea_4_532.tif	#EXTREMELY BRIGHT ARRAY	#high uniformity mean	#V HIGH BACKGROUND (MEAN EMPTY)	#EXTREMELY HIGH MEAN RANDOM	#HIGH MEAN EXPERIMENTAL
265945A02_X4_U01_Centaurea_4_532.tif	#EXTREMELY BRIGHT ARRAY	#bad washing top left	#high uniformity mean	#V HIGH BACKGROUND (MEAN EMPTY)	#EXTREMELY HIGH MEAN RANDOM	#HIGH MEAN EXPERIMENTAL
265945A04_X4_U01_Centaurea_3_532.tif	#EXTREMELY BRIGHT ARRAY	#high uniformity mean	#V HIGH BACKGROUND (MEAN EMPTY)	#EXTREMELY HIGH MEAN RANDOM	#HIGH MEAN EXPERIMENTAL
266698A01_X4_U01_Centaurea_532.tif	#dark array, dark spots in a few places. throw out?
266700A04_X4_U01_Centaurea_532_2.tif	#fuzzy spot

#BAD METRICS ONLY. CHECK AFTER NORMALIZATION
263785A03_bis_X4_U01_Centaurea_532.tif	#EXTREMELY HIGH MEAN RANDOM
263798A03_X4_U01_Centaurea_532.tif	#EXTREMELY HIGH MEAN RANDOM
265945A03_X4_U01_Centaurea_532.tif	#EXTREMELY HIGH MEAN RANDOM
261411A01_X4_U01_Centaurea_532.tif	#high uniformity mean	#HIGH MEAN EXPERIMENTAL
265944A01_X4_U01_Centaurea_2_532.tif	#EXTREMELY BRIGHT ARRAY. non-specific binding?	#EXTREMELY HIGH MEAN RANDOM
266700A01_X4_U01_Centaurea_532_2.tif	#bright array	#V HIGH BACKGROUND (MEAN EMPTY)	#EXTREMELY HIGH MEAN RANDOM
266700A01_X4_U01_Centaurea_532.tif	#bright array	#V HIGH BACKGROUND (MEAN EMPTY)	#EXTREMELY HIGH MEAN RANDOM
266700A03_X4_U01_Centaurea_532_2.tif	#bright array	#EXTREMELY HIGH MEAN RANDOM
266700A03_X4_U01_Centaurea_532.tif	#bright array	#EXTREMELY HIGH MEAN RANDOM
299709A01_X4_U01_Centaurea_532_2.tif	#gradient	#normalization should fix	#EXTREMELY HIGH MEAN RANDOM
299709A01_X4_U01_Centaurea_532.tif	#gradient	#normalization should fix	#EXTREMELY HIGH MEAN RANDOM
299709A04_X4_U01_Centaurea_532_2.tif	#gradient	#normalization should fix	#EXTREMELY HIGH MEAN RANDOM
299709A04_X4_U01_Centaurea_532.tif	#gradient	#normalization should fix	#EXTREMELY HIGH MEAN RANDOM


#BAD SPOTS ONLY
220768A01_X4_U01_Centaurea_532.tif	#dark spot top center, dark spot middle right
161 3 192 13; 196 296 203 301
260513A01_X4_U01_Centaurea_532.tif	#dark spot top center	#dust lower right
174 3 198 16; 287 281 304 313
260513A02_X4_U01_Centaurea_532.tif	#dark spot top left	#spot top left	#bright spot bottom center
124 5 141 21; 171 44 174 47; 155 457 161 461
260513A03_X4_U01_Centaurea_532.tif	#dark spot lower left
126 428 149 451
260513A04_X4_U01_Centaurea_532.tif	#dark spot bottom center
175 426 193 450
260948A01_X4_U01_Centaurea_532.tif	#dark spot top right
180 3 202 18
260961A01_X4_U01_Centaurea_532.tif	#dark spot top center
170 3 198 20
260961A02_X4_U01_Centaurea_532.tif	#dark spots top left
32 54 37 57; 121 3 146 17 
260966A01_X4_U01_Centaurea_532.tif	#dark spot top center
170 3 200 22
260967A02_X4_U01_Centaurea_532.tif	#dark spot top left	#dust top left
118 3 146 25; 140 24 156 40
260967A04_X4_U01_Centaurea_532.tif	#dark spot bottom center
165 426 183 446
261033A02_X4_U01_Centaurea_532.tif	#artifact upper left	#dark spot top center
1 31 14 47; 145 3 158 11
261033A03_X4_bis_U01_Centaurea_532.tif	#dark spot bottom right
132 432 157 456
261220A01_X4_U01_Centaurea_532.tif	#bright area top left
5 5 120 14
261415A01_X4_U01_Centaurea_532.tif	#dark spot top center
167 1 219 32
261415A02_X4_U01_Centaurea_532.tif	#bright spot upper left
106 74 156 91
261415A04_X4_U01_Centaurea_532.tif	#dark spots lower half
170 355 191 374; 187 425 210 447; 239 427 250 438
263797A03_X4_U01_Centaurea_532.tif	#dark spot bottom center
131 434 150 454
263797A04_X4_U01_Centaurea_532.tif	#dark spot bottom center
170 442 183 457
263798A01_X4_U01_Centaurea_532.tif	#bubbles upper left
70 37 76 43; 47 71 50 73
263798A02_X4_U01_Centaurea_532.tif	#dark spot top center
133 3 153 18
263798A04_X4_U01_Centaurea_532.tif	#dark spot bottom center
177 426 204 455
266698A04_X4_U01_Centaurea_532.tif	#dark spot bottom center
177 432 196 447
299628A02_X4_U01_Centaurea_532_2.tif	#fuzzy spot
243 149 253 159
299628A02_X4_U01_Centaurea_532.tif	#fuzzy spot
243 149 253 159
