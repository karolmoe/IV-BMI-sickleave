********************************************************************************
* Code to do IV & conventional analyses of BMI and risk of sick leave
* Adapted from code sent by David Carslake on 28/6/23
* Requires: `datafile1'; the main data file
* Preparation script for `datafile1'/"iv_BMI_sickleave.dta" is not provided
* Takes ~5 minutes to run
* Karoline Moe, october 2024
********************************************************************************

use "\\fil.nice.ntnu.no\nice\p499\Karoline\Arbeidsfiler\Do-filer\Artikkel 2- zBMI and risk of sick leave using offspring zBMI s an instrument\iv_BMI_sickleave.dta", clear

***************************************************************************************
*Table 1
{
*women
tab pmZbmi_cat, sum(part_age)
tab pmZbmi_cat, sum(bmi_use_Z6)
tab pmZbmi_cat, sum(female_bmi)
tab low_edu pmZbmi_cat, col 
tab pmZbmi_cat, sum(Occupation6_EGP)
tab Occupation_bw pmZbmi_cat, col
tab smoking pmZbmi_cat, col
tab alcohol pmZbmi_cat, col
tab phys_act pmZbmi_cat, col
*offspring of women
tab offspring_pmZbmi_cat, sum(part_age)
tab offspring_pmZbmi_cat, sum(bmi_use_Z6)
tab low_edu offspring_pmZbmi_cat, col 
tab offspring_pmZbmi_cat, sum(Occupation6_EGP)
tab Occupation_bw offspring_pmZbmi_cat, col
tab smoking offspring_pmZbmi_cat, col
tab alcohol offspring_pmZbmi_cat, col
tab phys_act offspring_pmZbmi_cat, col
*men
tab pfZbmi_cat, sum(part_age)
tab pfZbmi_cat, sum(bmi_use_Z6)
tab low_edu pfZbmi_cat, col
tab pfZbmi_cat, sum(Occupation6_EGP) 
tab Occupation_bw pfZbmi_cat, col
tab smoking pfZbmi_cat, col
tab alcohol pfZbmi_cat, col
tab phys_act pfZbmi_cat, col
*offspring of men
tab offspring_pfZbmi_cat, sum(part_age)
tab offspring_pfZbmi_cat, sum(bmi_use_Z6)
tab low_edu offspring_pfZbmi_cat, col 
tab offspring_pfZbmi_cat, sum(Occupation6_EGP)
tab Occupation_bw offspring_pfZbmi_cat, col
tab smoking offspring_pfZbmi_cat, col
tab alcohol offspring_pfZbmi_cat, col
tab phys_act offspring_pfZbmi_cat, col
}

***************************************************************************************
*Table 2, and IV estimates Table 3
{
*unadjusted
preserve
keep if Sex==1 //men
*keep if Sex==0 //women
local cov = 0
stset surv_time, failure(sf31_MSD==1) id(id) enter(time study_start) origin(BirthDate) scale(365.25)
*stset surv_time, failure(sf31_mental==1) id(id) enter(time study_start) origin(BirthDate) scale(365.25)
*stset surv_time, failure(sf31_allcause==1) id(id) enter(time study_start) origin(BirthDate) scale(365.25)
stcox offspring_z_bmi c.age_cspl*#i.wave
local B_num = _b[offspring_z_bmi]
local SE_num = _se[offspring_z_bmi]
scalar numerator = _b[offspring_z_bmi]
* regress BMIp against covariates only (not BMIo) and save R2:
regress z_bmi c.age_cspl*#i.wave
local nullR2=e(r2)
display `nullR2'
* regress BMIp against BMIo and covariates (the IV denominator):
regress z_bmi offspring_z_bmi c.age_cspl*#i.wave
* Get the partial R squared stat:
local partialR2=e(r2)-`nullR2'
display `partialR2'
* Test the null hypothesis for the instrument in the model above and return the F stat:
test offspring_z_bmi=0
local partialF=r(F) 
display `partialF'
*IV estimates Table 3 (needs to be exp)
local B_denom = _b[offspring_z_bmi]
local SE_denom = _se[offspring_z_bmi]
scalar denominator =_b[offspring_z_bmi]
display scalar(numerator)/scalar(denominator) 
scalar B_OAI = (`B_num')/(`B_denom')
scalar SE_OAI = sqrt(((`SE_num')^2/(`B_denom')^2) + ((`B_num')^2/(`B_denom')^4)*(`SE_denom')^2 - 2*((`B_num')/(`B_denom')^3)*`cov')
display B_OAI
display SE_OAI
scalar LCL_OAI = B_OAI-invnormal(0.975)*SE_OAI
scalar UCL_OAI = B_OAI+invnormal(0.975)*SE_OAI
display LCL_OAI
display UCL_OAI
restore


*adjusted 
preserve
keep if Sex==1
*keep if Sex==0 //women
local cov = 0
stset surv_time, failure(sf31_MSD==1) id(id) enter(time study_start) origin(BirthDate) scale(365.25)
*stset surv_time, failure(sf31_mental==1) id(id) enter(time study_start) origin(BirthDate) scale(365.25)
*stset surv_time, failure(sf31_allcause==1) id(id) enter(time study_start) origin(BirthDate) scale(365.25)
stcox offspring_z_bmi c.age_cspl*#i.wave i.smoking i.alcohol i.phys_act i.edu i.Occupation6_EGP
local B_num = _b[offspring_z_bmi]
local SE_num = _se[offspring_z_bmi]
scalar numerator = _b[offspring_z_bmi]
* regress BMIp against covariates only (not BMIo) and save R2:
regress z_bmi c.age_cspl*#i.wave i.smoking i.alcohol i.phys_act i.edu i.Occupation6_EGP
local nullR2=e(r2)
display `nullR2'
* regress BMIp against BMIo and covariates (the IV denominator):
regress z_bmi offspring_z_bmi c.age_cspl*#i.wave i.smoking i.alcohol i.phys_act i.edu i.Occupation6_EGP
* Get the partial R squared stat:
local partialR2=e(r2)-`nullR2'
display `partialR2'
* Test the null hypothesis for the instrument in the model above and return the F stat:
test offspring_z_bmi=0
local partialF=r(F) 
display `partialF'
*IV estimates Table 3 (needs to be exp)
local B_denom = _b[offspring_z_bmi]
local SE_denom = _se[offspring_z_bmi]
scalar denominator =_b[offspring_z_bmi]
display scalar(numerator)/scalar(denominator) 
scalar B_OAI = (`B_num')/(`B_denom')
scalar SE_OAI = sqrt(((`SE_num')^2/(`B_denom')^2) + ((`B_num')^2/(`B_denom')^4)*(`SE_denom')^2 - 2*((`B_num')/(`B_denom')^3)*`cov')
display B_OAI
display SE_OAI
scalar LCL_OAI = B_OAI-invnormal(0.975)*SE_OAI
scalar UCL_OAI = B_OAI+invnormal(0.975)*SE_OAI
display LCL_OAI
display UCL_OAI
restore
}

***************************************************************************************
*Table 3
{
*Own BMI
preserve
keep if Sex==1 //men
*keep if Sex==0 //women
stset surv_time, failure(sf31_MSD==1) id(id) enter(time study_start) origin(BirthDate) scale(365.25)
*stset surv_time, failure(sf31_mental==1) id(id) enter(time study_start) origin(BirthDate) scale(365.25)
*stset surv_time, failure(sf31_allcause==1) id(id) enter(time study_start) origin(BirthDate) scale(365.25)
stcox z_bmi
stcox z_bmi c.age_cspl*#i.wave
stcox z_bmi c.age_cspl*#i.wave i.smoking i.alcohol i.phys_act i.edu i.Occupation6_EGP
stptime, per(1000)
restore
*Offsprnig BMI
preserve
keep if Sex==1 //men
*keep Sex==0 //women
stset surv_time, failure(sf31_MSD==1) id(id) enter(time study_start) origin(BirthDate) scale(365.25)
*stset surv_time, failure(sf31_mental==1) id(id) enter(time study_start) origin(BirthDate) scale(365.25)
*stset surv_time, failure(sf31_allcause==1) id(id) enter(time study_start) origin(BirthDate) scale(365.25)
stcox offspring_z_bmi
stcox offspring_z_bmi c.age_cspl*#i.wave 
stcox offspring_z_bmi c.age_cspl*#i.wave i.smoking i.alcohol i.phys_act i.edu i.Occupation6_EGP
restore
}



*----------
* Settings:
*----------
* Is this the real data ("real") or mock data ("real")?
local dataset = "real"
* Define the data files and working directory on the basis of the response above:
if "`dataset'"=="real"{
	local datafile1 = "\\fil.nice.ntnu.no\nice\p499\Karoline\Arbeidsfiler\Do-filer\Artikkel 2- zBMI and risk of sick leave using offspring zBMI s an instrument\iv_BMI_sickleave.dta" 
	local working_dir = "\\fil.nice.ntnu.no\nice\p499\Karoline\Arbeidsfiler\Do-filer\Artikkel 2- zBMI and risk of sick leave using offspring zBMI s an instrument"
	
if "`dataset'"=="mock"{
	local datafile1 = "\\fil.nice.ntnu.no\nice\p499\Karoline\Arbeidsfiler\Do-filer\Artikkel 2- zBMI and risk of sick leave using offspring zBMI s an instrument\iv_BMI_sickleave.dta"
	local working_dir = "\\fil.nice.ntnu.no\nice\p499\Karoline\Arbeidsfiler\Do-filer\Artikkel 2- zBMI and risk of sick leave using offspring zBMI s an instrument"
	}
* List the outcomes (causes of death) which are to be analysed:
local outcomes = "sf31_allcause"
* List the variables to be summarised in Table 1:
local descriptive_vbls = "psmbmi69 psfbmi69 sbmi23 low_edu smoking alcohol phys_act Occupation_bw"
* List and define the adjustment sets:
local adjustments = "a b c d"
local adjustment_a = "c.BirthDate"
local adjustment_b = "c.BirthDate i.low_edu i.smoking i.alcohol i.phys_act i.Occupation_bw"
local adjustment_c = "c.BirthDate i.low_edu i.smoking i.alcohol i.phys_act i.Occupation_bw"
local adjustment_d = "c.BirthDate "
* Set a minimum number of deaths for analysis to proceed:
local min_deaths = 10
* For how many iterations should stcox attempt convergence?
set maxiter 100
* Set the random seed for reproducibility:
set seed 123456
* Run code as Stata 14:
version 14
* Install or update the metan package:
ssc install metan, replace

*---------------
* Preliminaries:
*---------------
clear all
cd "`working_dir'"
capture log close
log using "BMI_NCDS1958.log", replace
display c(current_date) ", " c(current_time)
use "`datafile1'", clear
save delme.dta, replace
* If necessary, simulate variable "serial" (ID number) for each parent (assuming no duplicates):
capture confirm variable serial
display _rc
if _rc!=0 generate serial = _n
* Set up places to store the output:
tempname memhold_F1
postfile `memhold_F1' str2 Dataset str1 BMI_person str8 Csex str1 Outcome_parent str15 Outcome str1 Adj BMIcat BMImean Beta LCL UCL Deaths N str80 Adjustment BMImean_event id bmi sf31_allcause Sex using "F1.dta", replace
tempname memhold_F2
postfile `memhold_F2' str20 Variable str1 Parent Beta_OAI SE_OAI Beta_CA SE_CA using "F2.dta", replace

*------------------------------
* Generate follow-up variables:
*------------------------------
* The origin is always the parent's DOB:
generate P_DOB = BirthDate
* Follow-up wrt offspring BMI can start at offspring DOB (actually 3-9/3/1958):
generate Start1a = offspring_BirthDate
* Follow-up wrt parental BMI can start at parental BMI measurement (actually :
generate Start2a = study_start
generate Start3a = study_start
* Start of follow-up is delayed if any covariates requiring a living parent were recorded later:
forvalues i = 1/3{
	* Although some covariates in sets b, c & d are recorded after 1958, I don't THINK they required the parent to be alive.
	* (This includes fsmok2 because missing data are included as a category):
	generate Start`i'b = Start`i'a
	generate Start`i'c = Start`i'a
	generate Start`i'd = Start`i'a
	* If they did, the start of follow up would be delayed as follows:
	*replace Start`i'd = date("01/07/1974", "DMY") if Start`i'd < date("01/07/1974", "DMY")
	}
* (The end of follow-up variable is fdate)
gen fdate= surv_time



*------------------------------
* Generate inclusion variables:
*------------------------------
*capital M is Sex (male) and small m is parent (mother)
gen sub="F" if Sex==0
replace sub="M" if Sex==1
gen pmbmi69=bmi if sub=="F"
gen pfbmi69=bmi if sub=="M"
rename offspring_bmi cbmi23
rename offspring_Sex csex 
gen pmschool=low_edu if sub=="F"
gen pfschool=low_edu if sub=="M"
*tabulate sub, missing
* INCLUSION VARIABLE Use1 requires offspring BMI:

generate Use1 = 1
replace Use1 = 0 if cbmi23>=.
tabulate sub Use1, missing
* INCLUSION VARIABLE Use2 requires parental BMI:

generate Use2 = 1
replace Use2 = 0 if (sub=="F" & pmbmi69>=.) | (sub=="M" & pfbmi69>=.)
tabulate sub Use2, missing
* INCLUSION VARIABLE Use3 requires offspring and parental BMI:
generate Use3 = 1
replace Use3 = 0 if Use1==0 | Use2==0
tabulate sub Use3, missing

forvalues i = 1/3{
	display "INCLUSION VARIABLE Use`i'a requires biological parents at follow-up and non-missing data on stratifying variables & covariates a:"
	generate Use`i'a = Use`i'
	display "INCLUSION VARIABLE Use`i'b requires Use`i'a==1, plus non-missing data on covariates b:"
	generate Use`i'b = Use`i'a
	display "INCLUSION VARIABLE Use`i'c requires Use`i'b==1, plus non-missing data on offspring smoking at 23:"
	generate Use`i'c = Use`i'b
	display "INCLUSION VARIABLE Use`i'd equals Use`i'b, because missing fsmoke is counted as a category:"
	generate Use`i'd = Use`i'b
	tabulate sub Use`i'd, missing
	drop Use`i'
	}

*--------------------------------------------
* Import and create some necessary variables:
*--------------------------------------------
* Make a single variable for parental BMI Z score in 1969:
codebook z_bmi 
rename z_bmi bmi69_Z
*Make z-scores separate for males and females
generate  psmbmi69 = bmi69_Z if sub=="F"
generate  psfbmi69 = bmi69_Z if sub=="M"
* Make a single variable for parental BMI Z score categories in 1969:
generate bmi69_Zg3 = pmZbmi_cat if sub=="F"
replace bmi69_Zg3 = pfZbmi_cat if sub=="M"
* Rename offspring_Zbmi_cat for consistency:
rename offspring_Zbmi_cat sbmi23g3
* Create squares of the exposure and instrument:
rename offspring_z_bmi sbmi23 
generate sbmi23_squared = sbmi23^2
generate bmi69_Z_squared = bmi69_Z^2
* Put BMI into plot categories:
generate oBMI_plotcat = floor(10*cbmi23)
recode oBMI_plotcat (min/184=1) (185/199=2) (200/249=3) (250/299=4) (300/349=5) (350/max=6) (.=.)
generate mBMI_plotcat = floor(10*pmbmi69)
recode mBMI_plotcat (min/184=1) (185/199=2) (200/249=3) (250/299=4) (300/349=5) (350/max=6) (.=.)
generate fBMI_plotcat = floor(10*pfbmi69)
recode fBMI_plotcat (min/184=1) (185/199=2) (200/249=3) (250/299=4) (300/349=5) (350/max=6) (.=.)
generate BMI_plotcat=mBMI_plotcat if sub=="F"
replace BMI_plotcat=fBMI_plotcat if sub=="M"
* Make a binary variable for parental smoking:
*smoking 
gen pmsmoking=smoking if sub=="F"
gen pfsmoking=smoking if sub=="M"
* Make a binary variable for parental alcohol:
*alcohol 
gen pmalcohol=alcohol if sub=="F"
gen pfalcohol=alcohol if sub=="M"
* Make a binary variable for parental physical activity:
*phys_act 
gen pmphys_act=phys_act if sub=="F"
gen pfphys_act=phys_act if sub=="M"
*Make binary variable for maternal and paternal HADS 
gen pmHADS_anxi = HADS_anxi if sub=="F" 
gen pfHADS_anxi = HADS_anxi if sub=="M"

gen pmHADS_dep = HADS_dep if sub=="F" 
gen pfHADS_dep = HADS_dep if sub=="M" 
* Make a variable for parental Occupation:
gen pmOccupation6_EGP=Occupation6_EGP if sub=="F"
gen pfOccupation6_EGP=Occupation6_EGP if sub=="M"

gen pmOccupation_bw=Occupation_bw if sub=="F"
gen pfOccupation_bw=Occupation_bw if sub=="M"

rename part_age page

foreach parent in "m" "f"{
		if "`parent'"=="m" local sub = "F"
		if "`parent'"=="f" local sub = "M"
}
foreach csex in "0" "1"{
		if "`csex'"=="0" local o = "F"
		if "`csex'"=="1" local o = "M"
}
**************************************************************************************************************************************************
*-----------------
* Do the analyses for Figure 4:
*-----------------
* Loop through datasets & parents:
foreach ds in "1a" "1b" "1c" "1d" "2a" "2b" "2c" "2d" "3a" "3b" "3c" "3d"{
	foreach parent in "m" "f"{
		if "`parent'"=="m" local sub = "F"
		if "`parent'"=="f" local sub = "M"
		foreach csex_condition in "csex==1" "csex==2" "csex<."{
						* Loop through BMI people:
			foreach bmiperson in "`parent'" "o"{
				if "`bmiperson'"=="o" local exposure = "sbmi23"
				if "`bmiperson'"=="o" local rawbmi = "cbmi23"
				if "`bmiperson'"=="`parent'" local exposure = "bmi69_Z"
				if "`bmiperson'"=="`parent'" local rawbmi = "p`parent'bmi69"
						*Do the analyses for F1:
						foreach outcome of local outcomes{
					stset fdate if sub=="`sub'" & Use`ds'==1 & `csex_condition', id(serial) failure(`outcome') enter(Start`ds') origin(P_DOB) scale(365.25)
					generate f_u = _t-_t0
					summarize f_u if sub=="`sub'" & Use`ds'==1 & `csex_condition', detail
					display "Median follow-up for parent sex `sub', dataset `ds', `csex_condition', bmiperson `bmiperson': "r(p50)
					drop f_u
						foreach adj of local adjustments{
							display "Doing F1 & FS1 analyses on outcome `outcome' in parent `parent' against BMI in bmiperson `bmiperson' with adjustment `adj', dataset `ds':"
							capture noisily stcox ib3.`bmiperson'BMI_plotcat `adjustment_`adj'' if sub=="`sub'" & Use`ds'==1 & `csex_condition', strata(csex) nohr

						if _rc==0{
							forvalues cat = 1/6{
								count if sub=="`sub'" & Use`ds'==1 & `csex_condition' & `bmiperson'BMI_plotcat<. & `bmiperson'BMI_plotcat==`cat' & `outcome'==1
								local D`cat' = r(N)
								count if sub=="`sub'" & Use`ds'==1 & `csex_condition' & `bmiperson'BMI_plotcat<. & `bmiperson'BMI_plotcat==`cat'
								local N`cat' = r(N)
								summarize `rawbmi' if sub=="`sub'" & Use`ds'==1 & `csex_condition' & `bmiperson'BMI_plotcat<. & `bmiperson'BMI_plotcat==`cat'
								local BMImean = r(mean)
								summarize `rawbmi' if sub=="`sub'" & Use`ds'==1 & `csex_condition' & `bmiperson'BMI_plotcat<. & `bmiperson'BMI_plotcat==`cat' & sf31_allcause==1
								local BMImean_event = r(mean)
								local Beta = .
								local LCL = .
								local UCL = .
								capture display _b[`cat'.`bmiperson'BMI_plotcat]
								if _rc==0{
									local Beta = _b[`cat'.`bmiperson'BMI_plotcat]
									local LCL = `Beta'+invnormal(0.025)*_se[`cat'.`bmiperson'BMI_plotcat]
									local UCL = `Beta'+invnormal(0.975)*_se[`cat'.`bmiperson'BMI_plotcat]
									}
								post `memhold_F1' ("`ds'") ("`bmiperson'") ("`csex_condition'") ("`parent'") ("`outcome'") ("`adj'") (`cat') (`BMImean') (`Beta') (`LCL') (`UCL') (`D`cat'') (`N`cat'') ("`adjustment_`adj''") (`BMImean_event') (id) (bmi) (sf31_allcause) (Sex) 
									}
								}
							}
						}
					}			
				}
			}
		}
postclose `memhold_F1'

		

*-----------------------------------------
* Calculate bias components for Figure 3, and estimate the effect of BMI on the vbl by the OAI method (Table 1):
*-----------------------------------------
{
tempname memhold_F2
postfile `memhold_F2' str20 Variable str1 Parent Beta_OAI SE_OAI Beta_CA SE_CA using "F2.dta", replace
display c(current_date)+": "+c(current_time)
recode alcohol (5 6 =1) (1 2 3 4=0)
recode phys_act (1 2 = 1) (3 4=0)
recode smoking (1 2 =0) (3 4=1)
local bcp_continuous = "page"
local bcp_binary = "low_edu smoking alcohol phys_act Occupation_bw"
local bcp_all = "`bcp_continuous' `bcp_binary'"
local cov = 0
* Calculate the bias components:
foreach vbl of local bcp_all{
	local type = cond(regexm(" `bcp_continuous' "," `vbl' ")==1,"continuous",cond(regexm(" `bcp_binary' "," `vbl' ")==1,"binary",""))
	assert "`type'"=="continuous" | "`type'"=="binary"
	local command = cond("`type'"=="continuous","regress",cond("`type'"=="binary","logistic","ERROR"))
	foreach parent in "m" "f"{
		if "`parent'"=="m" local sub = "F"
		if "`parent'"=="f" local sub = "M"
		* Check that the variable varies in that parent:
		quietly summarize `vbl' if sub=="`sub'"
		if r(sd)==0 | missing(r(sd)){
			display "Missing or invariant data for variable `vbl' and parent `parent'"
			post `memhold_F2' ("`vbl'") ("`parent'") (.) (.) (.) (.)	
			}
		else{
			* Estimate the effect of BMI on the vbl by the OAI method:
			`command' `vbl' sbmi23 if sub=="`sub'"
			local B_num = _b[sbmi23]
			local SE_num = _se[sbmi23]
			regress bmi69_Z sbmi23 if sub=="`sub'" 
			local B_denom = _b[sbmi23]
			local SE_denom = _se[sbmi23]
			local B_OAI = (`B_num')/(`B_denom')
			local SE_OAI = sqrt(((`SE_num')^2/(`B_denom')^2) + ((`B_num')^2/(`B_denom')^4)*(`SE_denom')^2 - 2*((`B_num')/(`B_denom')^3)*`cov')
			* Estimate the effect of BMI on the vbl by the CA method:
			`command' `vbl' bmi69_Z if sub=="`sub'" 
			local B_CA = _b[bmi69_Z]
			local SE_CA = _se[bmi69_Z]
			post `memhold_F2' ("`vbl'") ("`parent'") (`B_OAI') (`SE_OAI') (`B_CA') (`SE_CA')	
			}
		}
	}
}
postclose `memhold_F2'


*----------------------------
* Plot the data for Figure 4:
*----------------------------
use "\\fil.nice.ntnu.no\nice\p499\Karoline\Arbeidsfiler\Do-filer\Artikkel 2- zBMI and risk of sick leave using offspring zBMI s an instrument\iv_BMI_sickleave.dta", clear 
replace bmi=15 if bmi<15
replace bmi=40 if bmi>40
keep if sf31_allcause==1 & Sex==0 & bmi>=15 & bmi<=40
keep id bmi sf31_allcause Sex 
save "\\fil.nice.ntnu.no\nice\p499\Karoline\Arbeidsfiler\Do-filer\Artikkel 2- zBMI and risk of sick leave using offspring zBMI s an instrument\plot1_female.dta", replace
use "\\fil.nice.ntnu.no\nice\p499\Karoline\Arbeidsfiler\Do-filer\Artikkel 2- zBMI and risk of sick leave using offspring zBMI s an instrument\iv_BMI_sickleave.dta", clear 
replace bmi=15 if bmi<15
replace bmi=40 if bmi>40
keep if sf31_allcause==1 & Sex==1 & bmi>=15 & bmi<=40
keep id bmi sf31_allcause Sex 
save "\\fil.nice.ntnu.no\nice\p499\Karoline\Arbeidsfiler\Do-filer\Artikkel 2- zBMI and risk of sick leave using offspring zBMI s an instrument\plot1_male.dta", replace
local ds1 = "1b"
local ds2 = "2b"
local ds4 = "4"
local adj = "b"
local outcome = "sf31_allcause"
local csex_condition = "csex<."
use F1.dta, clear
append using "\\fil.nice.ntnu.no\nice\p499\Karoline\Arbeidsfiler\Do-filer\Artikkel 2- zBMI and risk of sick leave using offspring zBMI s an instrument\plot1_female.dta" 
append using "\\fil.nice.ntnu.no\nice\p499\Karoline\Arbeidsfiler\Do-filer\Artikkel 2- zBMI and risk of sick leave using offspring zBMI s an instrument\plot1_male.dta"
replace BMI_person="m" if BMI_person=="" & Sex==0
replace BMI_person="f" if BMI_person=="" & Sex==1
replace Dataset="4" if Dataset==""
replace Adj="b" if Adj==""
replace Outcome = "sf31_allcause" if Outcome==""
replace Outcome_parent="f" if Outcome_parent=="" & BMI_person=="f"
replace Outcome_parent="m" if Outcome_parent=="" & BMI_person=="m"
replace Csex= "csex<." if Csex==""
generate HR = exp(Beta)
generate LCL_HR = exp(LCL)
generate UCL_HR = exp(UCL)
summarize HR
summarize LCL_HR
summarize UCL_HR
gen LCL_HR2= cond(LCL_HR <0.50, 0.50, LCL_HR)
gen UCL_HR2= cond(UCL_HR >3, 3, UCL_HR)
*egen events= total(Deaths), by(BMImean)
*egen BMImean_events= BMImean if 
set graphics off
foreach parent in "f" "m"{
	if "`parent'"=="f"{
		local posh_parent = "Father" 
		local title = "B"
		}
	if "`parent'"=="m"{
		local posh_parent = "Mother"
		local title = "A"
		}
	#delimit ;
	graph twoway
	(rcap LCL_HR2 UCL_HR2 BMImean if BMI_person=="`parent'" & Dataset=="`ds2'", lcolor(gs0) lwidth(thin) msize(tiny))
	(line HR BMImean if BMI_person=="`parent'" & Dataset=="`ds2'", lcolor(gs0) lwidth(thin) lpattern(solid))
	(hist bmi if BMI_person=="`parent'" & Dataset=="`ds4'", freq discrete color(red%0) yaxis(2))
	(scatter HR BMImean if BMI_person=="`parent'" & Dataset=="`ds2'", mcolor(gs0) msymbol(O) msize(small))
	(rcap LCL_HR2 UCL_HR2 BMImean if BMI_person=="o" & Dataset=="`ds1'", lcolor(gs10) lwidth(thin) msize(tiny))
	(line HR BMImean if BMI_person=="o" & Dataset=="`ds1'", lcolor(gs10) lwidth(thin) lpattern(dash))
	(scatter HR BMImean if BMI_person=="o" & Dataset=="`ds1'", mcolor(gs10) msymbol(O) msize(small))
	if Adj=="`adj'" & Outcome=="`outcome'" & Outcome_parent=="`parent'" & Csex=="`csex_condition'"
	,
	xsize(6) ysize(3) yscale(log range("0.5 3.2")) ylabel("1 2 3", angle(0) labsize(small)) ylabel("20 40 60 80 100", axis(2) angle(0) labsize(small)) xlabel(, angle(0) labsize(small)) ymtick("0.6(0.2)3.2")
	title("`title'", position(11) span) xtitle(Own or offspring body mass index (kg m{superscript:-2}), size(small)) ytitle(Hazard ratio (95% CI), size(small)) 
	legend(position(1) order(3 6) label(3 "`posh_parent''s own BMI") label(6 "Offspring BMI") size(small))
	scheme(s2mono) graphregion(color(white)) name("`parent'", replace);
	#delimit cr
	}
set graphics on
graph combine m f, cols(1) xsize(6) ysize(6) imargin(vsmall) iscale(1) scheme(s2mono) graphregion(color(white)) 


*-------------------------
* Make the BCP (Figure 3):
*-------------------------
use F2.dta, clear
save delme.dta, replace
local LCL_truncation = -2
local UCL_truncation = 2
set graphics off
* Scale coefficients and SE by the absolute value of whichever of the two coefficients (OAI or CA) is larger:
generate Scale = max(abs(Beta_OAI), abs(Beta_CA))
replace Beta_OAI = Beta_OAI/Scale
replace SE_OAI = SE_OAI/Scale
replace Beta_CA = Beta_CA/Scale
replace SE_CA = SE_CA/Scale
* Calculate the CI for the Coefficients:
quietly generate LCL_OAI = Beta_OAI-invnormal(0.975)*SE_OAI
quietly generate UCL_OAI = Beta_OAI+invnormal(0.975)*SE_OAI
quietly generate LCL_CA = Beta_CA-invnormal(0.975)*SE_CA
quietly generate UCL_CA = Beta_CA+invnormal(0.975)*SE_CA
* Truncate CI at the chosen level:
generate LCL_OAI_truncated = .
replace LCL_OAI_truncated = `LCL_truncation' if LCL_OAI<`LCL_truncation'
replace LCL_OAI = `LCL_truncation' if LCL_OAI<`LCL_truncation'
generate UCL_OAI_truncated = .
replace UCL_OAI_truncated = `UCL_truncation' if UCL_OAI>`UCL_truncation' & UCL_OAI<.
replace UCL_OAI = `UCL_truncation' if UCL_OAI>`UCL_truncation' & UCL_OAI<.
generate LCL_CA_truncated = .
replace LCL_CA_truncated = `LCL_truncation' if LCL_CA<`LCL_truncation'
replace LCL_CA = `LCL_truncation' if LCL_CA<`LCL_truncation'
generate UCL_CA_truncated = .
replace UCL_CA_truncated = `UCL_truncation' if UCL_CA>`UCL_truncation' & UCL_CA<.
replace UCL_CA = `UCL_truncation' if UCL_CA>`UCL_truncation' & UCL_CA<.
* Determine the X axis range:
generate Extreme_X = max(abs(LCL_OAI), abs(LCL_CA), abs(UCL_OAI), abs(UCL_CA))
summarize Extreme_X
local xmax = r(max)
local xlabel = ""
if `xmax'>=1 local xlabel = "xlabel(-1 -0.5 0 0.5 1, labsize(small)) xline(0, lcolor(gs12))"
if `xmax'>=2 local xlabel = "xlabel(-2 -1 0 1 2, labsize(small)) xline(0, lcolor(gs12))"
if `xmax'>=3 local xlabel = "xlabel(-3 -2 -1 0 1 2 3, labsize(small)) xline(0, lcolor(gs12))"
if `xmax'>=4 local xlabel = "xlabel(-4 -2 0 2 4, labsize(small)) xline(0, lcolor(gs12))"
if `xmax'>=6 local xlabel = "xlabel(-6 -3 0 3 6, labsize(small)) xline(0, lcolor(gs12))"
if `xmax'>=8 local xlabel = "xlabel(-8 -4 0 4 8, labsize(small)) xline(0, lcolor(gs12))"
if `xmax'>=10 local xlabel = "xlabel(-10 -5 0 5 10, labsize(small)) xline(0, lcolor(gs12))"
drop Extreme_X
local xmin = -1*`xmax'
*Define "posh" covariate names for plotting:
generate Posh_Covariate = Variable
*replace Posh_Covariate = "Date of birth" if Variable=="P_DOB"
replace Posh_Covariate = "Education below university level" if Variable=="low_edu"
replace Posh_Covariate = "Participation age" if Variable=="page"
replace Posh_Covariate = "Frequent alcohol use" if Variable=="alcohol"
replace Posh_Covariate = "Current smoker" if Variable=="smoking"
replace Posh_Covariate = "Low physical activity" if Variable=="phys_act"
replace Posh_Covariate = "Occupation EGP class ≥III" if Variable=="Occupation_bw"
* Determine the order (based on the higher of M & F in abs(RB)) in which Covariates will be plotted:
generate RB = abs(Beta_OAI)/abs(Beta_CA)
bysort Posh_Covariate (Parent): egen Max_RB = max(RB)
sort Max_RB Variable Parent
generate Y_OAI = ceil(_n/2)-0.15
generate Y_CA = ceil(_n/2)+0.15
* Create plots for each parent:
local N_covariates = _N/2
foreach parent in "f" "m"{
	* Sort out some Y axis options:
	local ylabels = ""
	forvalues c = 1/`N_covariates'{
		levelsof Posh_Covariate if round(Y_OAI,1)==`c', local(posh_covariate) clean
		local ylabels = `"`ylabels' `c' "`posh_covariate'""'
		}
	local ylabels = strtrim(`"`ylabels'"')
	if "`parent'"=="m" local ylabels = `"ylabel(`ylabels', angle(0) labsize(vsmall))"'
	if "`parent'"=="f" local ylabels = `"ylabel(`ylabels', angle(0) labsize(vsmall) nolabels)"'
	local ysize = `N_covariates'
	* Sort out some X axis options:
	if "`parent'"=="m"{
		local xtitle = "Mothers"
		local xsize = "xsize(4)"
		}
	if "`parent'"=="f"{
		local xtitle = "Fathers"
		local xsize = "fxsize(31)"
		}
	* Code the individual plots:
	local code_Beta_OAI = "(scatter Y_OAI Beta_OAI, mcolor(gs0) mfcolor(gs16) msymbol(O))"
	local code_CI_OAI = "(rcap LCL_OAI UCL_OAI Y_OAI, horizontal lcolor(gs0) msize(zero))"
	local code_Truncation_OAI = "(scatter Y_OAI LCL_OAI_truncated, mcolor(gs0) mfcolor(gs0) msymbol(X)) (scatter Y_OAI UCL_OAI_truncated, mcolor(gs0) mfcolor(gs0) msymbol(X))"
	local code_Beta_CA = "(scatter Y_CA Beta_CA, mcolor(gs0) mfcolor(gs0) msymbol(O))"
	local code_CI_CA = "(rcap LCL_CA UCL_CA Y_CA, horizontal lcolor(gs0) msize(zero))"
	local code_Truncation_CA = "(scatter Y_CA LCL_CA_truncated, mcolor(gs0) mfcolor(gs0) msymbol(X)) (scatter Y_CA UCL_CA_truncated, mcolor(gs0) mfcolor(gs0) msymbol(X))"
	local x_opts = `"`xlabel' xscale(range(`xmin' `xmax')) xtitle("Bias Component" "(`xtitle')", margin(0 0 0 1) size(medsmall)) `xsize'"'
	local y_opts = `"`ylabels' ytitle("") ysize(`ysize')"'
	local other_opts = `"scheme(s2mono) graphregion(color(white)) legend(off) name("BCP_`parent'", replace)"'
	local if_cond = `"if Parent=="`parent'""'
	graph twoway `code_CI_OAI' `code_Beta_OAI' `code_Truncation_OAI' `code_CI_CA' `code_Beta_CA' `code_Truncation_CA' `if_cond', `x_opts' `y_opts' `other_opts'
	}
* Make a plot with just the legend:
local plot1 = `"(scatter Y_OAI Beta_OAI if Parent=="XXXX", mcolor(gs0) mfcolor(gs16) msymbol(O))"'
local plot2 = `"(scatter Y_CA Beta_CA if Parent=="XXXX", mcolor(gs0) mfcolor(gs0) msymbol(O))"'
local legend_opts = `"legend(label(1 "Instrumental variable") label(2 "Conventional analysis") cols(1) region(fcolor(none) lstyle(none) lwidth(none) margin(zero)) ring(0) bplacement(s) size(vsmall) rowgap(*0.5))"'
local other_opts = `"name(legend, replace) xtitle(off) ytitle(off) title() scheme(s2mono) aspectratio(0) xsize(5) ysize(1) xscale(off) yscale(off) plotregion(margin(1 1 1 1) color(none)) graphregion(style(none) margin(zero) color(white))"'
graph twoway `plot1' `plot2' `plot3', `legend_opts' `other_opts'
_gm_edit .legend.plotregion1.draw_view.set_false
_gm_edit .legend.ystretch.set fixed
_gm_edit .legend.title.draw_view.set_false
_gm_edit .legend.yaxis1.draw_view.set_false
_gm_edit .legend.xaxis1.draw_view.set_false
* Combine the graphs:	
set graphics on
graph combine legend BCP_m BCP_f, cols(2) ysize(7) xsize(6) holes(2) scheme(s1mono) l1title("") b1title("") iscale(1) title("") imargin(tiny) name(BCP,replace) saving("F22.gph",replace)
graph export "F22.eps", as(eps) preview(off) replace
graph export "BCP2.pdf", as(pdf) replace
clear all		
}
