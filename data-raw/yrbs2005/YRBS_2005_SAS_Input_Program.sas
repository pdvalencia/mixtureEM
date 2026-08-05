
/****************************************************************************************/
/*  This SAS program reads ASCII format (text format) 2005 YRBS data and creates a      */
/*  formatted and labeled SAS dataset.                                                  */
/*                                                                                      */
/*  Change the file location specifications from "c:\yrbs2005" to the location where 	*/
/*  you have stored the YRBS ASCII data file and the format library before you run this */
/*  program.  Change the location specification in three places - in the "filename"     */
/*  statement and in the two "libname" statements at the top of the program.            */
/*                                                                                      */
/*  Note: Run "YRBS_2005_SAS_Format_Program.sas" BEFORE you run                         */
/*  "YRBS_2005_SAS_Input_Program.sas" to create the 2005 YRBS dataset.                  */
/****************************************************************************************/

filename datain 'c:\yrbs2005\yrbs2005.dat';

libname dataout 'c:\yrbs2005';

libname library 'c:\yrbs2005';

data dataout.yrbs2005 ;
infile datain lrecl=400;
input
Q1 $ 17-17
Q2 $ 18-18
Q3 $ 19-19
Q4 $ 20-21
@22 Q5 4.2
@26 Q6 6.2
Q7 $ 32-32
Q8 $ 33-33
Q9 $ 34-34
Q10 $ 35-35
Q11 $ 36-36
Q12 $ 37-37
Q13 $ 38-38
Q14 $ 39-39
Q15 $ 40-40
Q16 $ 41-41
Q17 $ 42-42
Q18 $ 43-43
Q19 $ 44-44
Q20 $ 45-45
Q21 $ 46-46
Q22 $ 47-47
Q23 $ 48-48
Q24 $ 49-49
Q25 $ 50-50
Q26 $ 51-51
Q27 $ 52-52
Q28 $ 53-53
Q29 $ 54-54
Q30 $ 55-55
Q31 $ 56-56
Q32 $ 57-57
Q33 $ 58-58
Q34 $ 59-59
Q35 $ 60-60
Q36 $ 61-61
Q37 $ 62-62
Q38 $ 63-63
Q39 $ 64-64
Q40 $ 65-65
Q41 $ 66-66
Q42 $ 67-67
Q43 $ 68-68
Q44 $ 69-69
Q45 $ 70-70
Q46 $ 71-71
Q47 $ 72-72
Q48 $ 73-73
Q49 $ 74-74
Q50 $ 75-75
Q51 $ 76-76
Q52 $ 77-77
Q53 $ 78-78
Q54 $ 79-79
Q55 $ 80-80
Q56 $ 81-81
Q57 $ 82-82
Q58 $ 83-83
Q59 $ 84-84
Q60 $ 85-85
Q61 $ 86-86
Q62 $ 87-87
Q63 $ 88-88
Q64 $ 89-89
Q65 $ 90-90
Q66 $ 91-91
Q67 $ 92-92
Q68 $ 93-93
Q69 $ 94-94
Q70 $ 95-95
Q71 $ 96-96
Q72 $ 97-97
Q73 $ 98-98
Q74 $ 99-99
Q75 $ 100-100
Q76 $ 101-101
Q77 $ 102-102
Q78 $ 103-103
Q79 $ 104-104
Q80 $ 105-105
Q81 $ 106-106
Q82 $ 107-107
Q83 $ 108-108
Q84 $ 109-109
Q85 $ 110-110
Q86 $ 111-111
Q87 $ 112-112
Q88 $ 113-113
Q89 $ 114-114
Q90 $ 115-115
Q91 $ 116-116
Q92 $ 117-117
Q93 $ 118-118
Q94 $ 119-119
Q95 $ 120-120
Q96 $ 121-121
Q97 $ 122-122
QN7 $ 123-123
QN8 $ 124-124
QN9 $ 125-125
QN10 $ 126-126
QN11 $ 127-127
QN12 $ 128-128
QN13 $ 129-129
QN14 $ 130-130
QN15 $ 131-131
QN16 $ 132-132
QN17 $ 133-133
QN18 $ 134-134
QN19 $ 135-135
QN20 $ 136-136
QN21 $ 137-137
QN22 $ 138-138
QN23 $ 139-139
QN24 $ 140-140
QN25 $ 141-141
QN26 $ 142-142
QN27 $ 143-143
QN28 $ 144-144
QN29 $ 145-145
QN30 $ 146-146
QN31 $ 147-147
QN32 $ 148-148
QN33 $ 149-149
QN34 $ 150-150
QN35 $ 151-151
QN36 $ 152-152
QN37 $ 153-153
QN38 $ 154-154
QN39 $ 155-155
QN40 $ 156-156
QN41 $ 157-157
QN42 $ 158-158
QN43 $ 159-159
QN44 $ 160-160
QN45 $ 161-161
QN46 $ 162-162
QN47 $ 163-163
QN48 $ 164-164
QN49 $ 165-165
QN50 $ 166-166
QN51 $ 167-167
QN52 $ 168-168
QN53 $ 169-169
QN54 $ 170-170
QN55 $ 171-171
QN56 $ 172-172
QN57 $ 173-173
QN58 $ 174-174
QN59 $ 175-175
QN60 $ 176-176
QN61 $ 177-177
QN62 $ 178-178
QN63 $ 179-179
QN64 $ 180-180
QN65 $ 181-181
QN66 $ 182-182
QN67 $ 183-183
QN68 $ 184-184
QN69 $ 185-185
QN70 $ 186-186
QN71 $ 187-187
QN72 $ 188-188
QN73 $ 189-189
QN74 $ 190-190
QN75 $ 191-191
QN76 $ 192-192
QN77 $ 193-193
QN78 $ 194-194
QN79 $ 195-195
QN80 $ 196-196
QN81 $ 197-197
QN82 $ 198-198
QN83 $ 199-199
QN84 $ 200-200
QN85 $ 201-201
QN86 $ 202-202
QN87 $ 203-203
QN88 $ 204-204
QN89 $ 205-205
QN90 $ 206-206
QN91 $ 207-207
QN92 $ 208-208
QN93 $ 209-209
QN94 $ 210-210
QN95 $ 211-211
QN96 $ 212-212
QN97 $ 213-213
QNFRCIG $ 214-214
QNFRVG $ 215-215
QNDLYPE $ 216-216
QNROVWGT $ 217-217
QNOVWGT $ 218-218
QNANYTOB $ 219-219
QNMINPA $ 220-220
QNNOPA $ 221-221
QNASATCK $ 222-222
@350 Q4ORIG $char8.
@358 Weight 12.7
@370 PSU 5.0
@375 Stratum 4.0
BMIPct $ 379-383
ETHORIG $ 384-384
RACEORIG $ 385-389
;

************************************
*  Assign formats to SAS variables.*
************************************;
format 
Q1 $H1S. Q2 $H2S. Q3 $H3S. Q4 $H4S. Q5 4.2 Q6 6.2 
Q7 $H7S. Q8 $H8S. Q9 $H9S. Q10 $H10S. 
Q11 $H11S. Q12 $H12S. Q13 $H13S. Q14 $H14S. Q15 $H15S. 
Q16 $H16S. Q17 $H17S. Q18 $H18S. Q19 $H19S. Q20 $H20S. 
Q21 $H21S. Q22 $H22S. Q23 $H23S. Q24 $H24S. Q25 $H25S. 
Q26 $H26S. Q27 $H27S. Q28 $H28S. Q29 $H29S. Q30 $H30S. 
Q31 $H31S. Q32 $H32S. Q33 $H33S. Q34 $H34S. Q35 $H35S. 
Q36 $H36S. Q37 $H37S. Q38 $H38S. Q39 $H39S. Q40 $H40S. 
Q41 $H41S. Q42 $H42S. Q43 $H43S. Q44 $H44S. Q45 $H45S. 
Q46 $H46S. Q47 $H47S. Q48 $H48S. Q49 $H49S. Q50 $H50S. 
Q51 $H51S. Q52 $H52S. Q53 $H53S. Q54 $H54S. Q55 $H55S. 
Q56 $H56S. Q57 $H57S. Q58 $H58S. Q59 $H59S. Q60 $H60S. 
Q61 $H61S. Q62 $H62S. Q63 $H63S. Q64 $H64S. Q65 $H65S. 
Q66 $H66S. Q67 $H67S. Q68 $H68S. Q69 $H69S. Q70 $H70S. 
Q71 $H71S. Q72 $H72S. Q73 $H73S. Q74 $H74S. Q75 $H75S. 
Q76 $H76S. Q77 $H77S. Q78 $H78S. Q79 $H79S. Q80 $H80S. 
Q81 $H81S. Q82 $H82S. Q83 $H83S. Q84 $H84S. Q85 $H85S. 
Q86 $H86S. Q87 $H87S. Q88 $H88XX. Q89 $H89XX. Q90 $H90XX. 
Q91 $H91XX. Q92 $H92XX. Q93 $H93XX. Q94 $H94XX. Q95 $H95XX. 
Q96 $H96XX. Q97 $H97XX.;

***********************************
*  Assign labels to SAS variables.*
***********************************;
label 
Q1="How old are you"
Q2="What is your sex"
Q3="In what grade are you"
Q4="How do you describe yourself"
Q5="How tall are you"
Q6="How much do you weigh"
Q7="How do you describe your health"
Q8="How often wear bicycle helmet"
Q9="How often wore  a seat belt"
Q10="How often ride w/drinking driver 30 days"
Q11="How often drive while drinking 30 days"
Q12="Carried weapon 30 days"
Q13="Carried gun 30 days"
Q14="Carried weapon at school  30 days"
Q15="How many days feel unsafe@school 30 days"
Q16="How  many times threatened@school 12 mos"
Q17="Property stolen at school"
Q18="How many times in fight 12 mos"
Q19="How many times injured in fight 12 mos"
Q20="How many times in fight @ school  12 mos"
Q21="Did  boyfriend/girlfriend hit/slap 12 mo"
Q22="Have you been forced to have sex"
Q23="Ever feel sad or hopeless 12 mos"
Q24="Ever considered suicide 12 mos"
Q25="Ever make suicide plan 12 mos"
Q26="Ever attempt suicide 12 mos"
Q27="Ever injured from suicide attempt 12 mos"
Q28="Ever smoked"
Q29="How old when first smoked"
Q30="How many days smoked 30 days"
Q31="How many cigarettes/day 30 days"
Q32="How did you get cigarettes past 30 days"
Q33="How many days smoke @ school 30 days"
Q34="Have you ever smoked daily"
Q35="Tried to quit smoking past 12 months"
Q36="How many days use snuff past 30 days"
Q37="Days use snuff school property 30 days"
Q38="How many days smoke cigars 30 days"
Q39="How many days drink alcohol"
Q40="How old when first drank alcohol"
Q41="How many days drink alcohol 30 days"
Q42="How many days have 5+ drinks 30 days"
Q43="How many days drink @ school 30 days"
Q44="How many times smoke marijuana"
Q45="How old when first tried marijuana"
Q46="How many times use marijuana 30 days"
Q47="How many times marijuana@school 30 days"
Q48="How many times use cocaine"
Q49="How many times use cocaine 30 days"
Q50="How many times sniffed glue"
Q51="How many times used heroin"
Q52="How many times used methamphetamines"
Q53="Ecstacy one or more time"
Q54="How many times used steroids"
Q55="How many times injected drugs"
Q56="Offered drugs @ school 12 mos"
Q57="Ever had sex"
Q58="How old at first sex"
Q59="How many sex partners"
Q60="How many sex partners 3 mos"
Q61="Did you use alcohol/drugs @ last sex"
Q62="Did you use condom @ last sex"
Q63="What birth control @ last sex"
Q64="How do you describe your weight"
Q65="What are you trying to do about weight"
Q66="Did you exercise to lose weight 30 days"
Q67="Did you eat less to lose weight 30 days"
Q68="Did you fast to lose weight 30 days"
Q69="Did you take pill to lose weight 30 days"
Q70="Did you vomit to lose weight 30 days"
Q71="How many times fruit juice 7 days"
Q72="How many times fruit 7 days"
Q73="How many time green salad 7 days"
Q74="How many times potatoes 7 days"
Q75="How many times carrots 7 days"
Q76="How many times other vegetables 7 days"
Q77="How many glass of milk 7 days"
Q78="Did you do vigorous exercise 7 days"
Q79="Did you do moderate exercise 7 days"
Q80="Days active 60 min plus past 7 days"
Q81="How many hours watch TV"
Q82="How many days go to PE class"
Q83="How many minutes exercise in PE class"
Q84="On how many sports team 12 mos"
Q85="Ever taught about AIDS/HIV @ school"
Q86="Ever been told you have asthma"
Q87="Asthma attack in past 12 months"
Q88="How often wear motorcycle helmet"
Q89="Show age proof buying cigarettes 30 days"
Q90="Times used hallucinogens"
Q91="Computer use per day"
Q92="Injured while playing sports"
Q93="Ever tested for HIV"
Q94="Wear sunscreen when outside"
Q95="Protection from the sun"
Q96="Disability/health problem"
Q97="Days missed class w/out permission"
QN7="Described health as fair/poor"
QN8="Never/rarely wore bicycle helmet"
QN9="Never/rarely wore seat belt"
QN10="Rode 1+ times with drinking driver"
QN11="Drove 1+ times when drinking"
QN12="Carried weapon 1+ times past 30 days"
QN13="Carried gun 1+ past 30 days"
QN14="Carried weapon school 1+ past 30 days"
QN15="Missed school b/c unsafe 1+ 30 days"
QN16="Threatened at school 1+ times 12 mos"
QN17="Prop stolen at school 12 mos"
QN18="Fought 1+ times 12 mos"
QN19="Injured/treated 1+ times 12 mos"
QN20="Fought school 1+ times 12 mos"
QN21="Hit by bf/gf 12 mos"
QN22="Forced to have sex"
QN23="Sad 2 wks past 12 mos"
QN24="Considered suicide 12 mos"
QN25="Made suicide plan 12 mos"
QN26="Attempted suicide 1+ times 12 mos"
QN27="Suicide attempt w/injury 12 mos"
QN28="Ever tried cigarettes"
QN29="Smoked cigarette before 13"
QN30="Smoked 1+ past 30 days"
QN31="10+ cigarettes/day past 30 days"
QN32="Got cigarettes in store 30 days"
QN33="Smoked at school 1+ past 30 days"
QN34="Smoked daily for 30 days"
QN35="Among smokers, tried to quit smoking"
QN36="Used snuff/dip 1+ past 30 days"
QN37="Used snuff/dip at school 1+ 30 days"
QN38="Smoked cigars 1+ past 30 days"
QN39="Had 1 drink on 1+ days in life"
QN40="Had first drink before 13"
QN41="Had 1+ drinks past 30 days"
QN42="Five+ drinks 1+ past 30 days"
QN43="Had 1+ drinks at school 1+ 30 days"
QN44="Tried marijuana 1+ times in life"
QN45="Tried marijuana before 13"
QN46="Used marijuana 1+ times past 30 days"
QN47="Used marijuana school 1+ times 30 day"
QN48="Used cocaine 1+ times in life"
QN49="Used cocaine 1+ times past 30 days"
QN50="Sniffed glue 1+ times in life"
QN51="Used heroin 1+ times in life"
QN52="Used meth 1+ times in life"
QN53="Used ecstasy 1+ times in life"
QN54="Took steroids 1+ times in life"
QN55="Injected drugs 1+ times in life"
QN56="Offered/sold drugs at school 12 mos"
QN57="Had sex ever"
QN58="Had sex before 13"
QN59="Had sex with 4+ people in life"
QN60="Had sex with 1+ people 3 mos"
QN61="Of current sex, used alcohol last time"
QN62="Of current sex, used condom last time"
QN63="Of current sex, used birth ctl last sx"
QN64="Slightly/very overweight"
QN65="Trying to lose weight"
QN66="Exercised to lose weight past 30 days"
QN67="Ate less to lose weight past 30 days"
QN68="Fasted to lose weight past 30 days"
QN69="Took pills to lose weight past 30 days"
QN70="Vomited to lose weight past 30 days"
QN71="Drank fruit juice past 7 days"
QN72="Ate fruit past 7 days"
QN73="Ate green salad past 7 days"
QN74="Ate potatoes past 7 days"
QN75="Ate carrots past 7 days"
QN76="Ate vegetables past 7 days"
QN77="Drank 3+ glasses milk past 7 days"
QN78="Exercised vigorously past 7 days"
QN79="Exercised moderately past 7 days"
QN80="Active 60 min on 5+ past 7 days"
QN81="Watched 3+ hours of TV average day"
QN82="Got to PE class 1+ days average week"
QN83="Of enrolled in PE, exercised 20 min"
QN84="Played on 1+ sports teams 12 mos"
QN85="Taught about AIDS in school"
QN86="Told by doctor/nurse they had asthma"
QN87="Told had asthma and current asthma"
QN88="Never/rarely wore motorcycle helmet"
QN89="Not asked proof age buying cigarettes"
QN90="Used hallucinogenic drugs"
QN91="Played video games 3+ hours per day"
QN92="Injured playing sports seen by doctor"
QN93="Tested for HIV"
QN94="Most time/always wore sunscreen outside"
QN95="Most time/always protect from sun"
QN96="Has a disability/health problem"
QN97="Missed class w/out permission 30 days"
qnfrcig="Smoked on 20 past 30 days"
qnfrvg="Ate 5+ fruits/vegetables 7 days"
qndlype="Attended PE class daily"
qnrovwgt="At risk for becoming overweight"
qnovwgt="Overweight"
qnanytob="Used any tobacco past 30 days"
qnminpa="Did not exercise 20 min 3+ past 7 days"
qnnopa="No exercise"
qnasatck="Had asthma attack past 12 months"
weight="Analysis weight"
stratum="Stratum"
psu="PSU"
bmipct="Body Mass Index Percentage"
q4orig="Race/ethnicity as originally scanned"
ethorig="Ethnicity as originally scanned"
raceorig="Race as originally scanned"
;
run;


