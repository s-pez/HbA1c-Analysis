/******************************************************************************
 PROJECT: Evaluating the HbA1c Screening Cutoff for Type 2 Diabetes
          A Sensitivity/Specificity, ROC, and Disparities Analysis

 PURPOSE: Simulates adult population health data calibrated to published NHANES 
          / CDC National Diabetes Statistics Report figures (overall diabetes 
          prevalence, mean HbA1c, and known racial/ethnic prevalence disparities), 
          then evaluates HbA1c as a screening test against a diagnostic reference
          standard - the same cutoff-tradeoff and ROC methodology used for
          a public health screening program deciding where to draw a line
          between "screen positive" and "screen negative."
          
          American Diabetes Association cutoffs were themselves
          set using ROC/Youden analysis against fasting plasma glucose 
          and oral glucose tolerance test reference standards. It is also a 
          well-documented case where a screening test's accuracy is not uniform
          across racial/ethnic groups, which makes it a strong vehicle for a 
          disparities-in-test-performance analysis, not just a 
          disparities-in-outcome analysis.

 DATA:    This uses simulated data calibrated to published
          national estimates (overall adult diabetes prevalence ~11-12%;
          mean HbA1c in the non-diabetic population ~5.3-5.5%; higher
          prevalence among Non-Hispanic Black and Mexican American adults
          per CDC National Diabetes Statistics Report).

 SOFTWARE: SAS 9.4
******************************************************************************/
/* REPLACE WITH YOUR FILE PATH:
*/
ODS RTF FILE = '/path/HbA1c_Project_Tables.rtf' STYLE= statistical;

/* REPLACE WITH YOUR FILE PATH:
*/
%let outpath = /path/HbA1c Project;

/*
 SECTION 1: SIMULATE NHANES-CALIBRATED ADULT POPULATION DATA
*/
data diabetes_screen;
    call streaminit(20260721);
    length race_eth $20;
    array racecat[4] $20 _temporary_
        ("Non-Hispanic White","Non-Hispanic Black","Mexican American","Other/Multiracial");
    array raceprob[4] _temporary_ (0.62, 0.14, 0.12, 0.12);
    array racediab_or[4] _temporary_ (1.00, 1.55, 1.45, 1.10);  /* relative risk multipliers, per published disparity pattern */

    do subj_id = 1 to 6000;

        /* --- assign race/ethnicity category --- */
        u = rand('UNIFORM');
        cum = 0;
        do r = 1 to 4;
            cum + raceprob[r];
            if u <= cum then do;
                race_eth = racecat[r];
                race_or  = racediab_or[r];
                leave;
            end;
        end;

        /* --- covariates --- */
        Age = round(min(79, max(20, rand('NORMAL', 47, 15))));
        BMI = round(min(55, max(16, rand('NORMAL', 29, 6))), 0.1);

        /* --- latent diabetes risk (logistic function calibrated to ~11.5% overall prevalence) --- */
        logit_p = -4.2 + 0.045*Age + 0.075*(BMI-25) + log(race_or);
        p_diab  = 1/(1+exp(-logit_p));
        True_Diabetes = (rand('UNIFORM') < p_diab);

        /* --- HbA1c: correlated with true status, age, BMI, plus test noise --- */
        if True_Diabetes = 0 then
            HbA1c = max(4.0, rand('NORMAL', 5.35 + 0.004*Age + 0.01*(BMI-25), 0.35));
        else
            HbA1c = max(5.0, rand('NORMAL', 7.6 + 0.01*Age, 1.3));
        HbA1c = round(HbA1c, 0.1);

        /* --- Fasting Plasma Glucose (FPG): correlated but imperfectly with HbA1c --- */
        if True_Diabetes = 0 then
            FPG = max(65, rand('NORMAL', 95 + 0.3*Age, 9));
        else
            FPG = max(90, rand('NORMAL', 155 + 0.5*Age, 35));
        FPG = round(FPG);

        /* --- screen result flags at candidate HbA1c cutoffs --- */
        Flag_57 = (HbA1c >= 5.7);
        Flag_60 = (HbA1c >= 6.0);
        Flag_63 = (HbA1c >= 6.3);
        Flag_65 = (HbA1c >= 6.5);   /* current ADA diagnostic cutoff */

        output;
    end;
    drop u cum r race_or logit_p p_diab;
run;

title "Simulated Adult Screening Data - Structure Check";
proc contents data=diabetes_screen varnum; run;
proc freq data=diabetes_screen; tables race_eth True_Diabetes; run;
title;


/*
 SECTION 2: DESCRIPTIVE / DISTRIBUTION ANALYSIS
*/
title "HbA1c Distribution by True Diabetes Status";
proc univariate data=diabetes_screen;
    class True_Diabetes;
    var HbA1c;
    histogram HbA1c / normal kernel;
    inset mean std min max / position=ne;
run;
title;

ods graphics on;
title "Overlapping HbA1c Distributions: Diabetic vs. Non-Diabetic Adults";
proc sgplot data=diabetes_screen;
    histogram HbA1c / group=True_Diabetes transparency=0.4 binwidth=0.2;
    density HbA1c / group=True_Diabetes type=kernel;
    xaxis label="HbA1c (%)" max=14;
    yaxis label="Density";
    keylegend / title="True Diabetes Status (0=No, 1=Yes)";
run;
title;

title "HbA1c vs. FPG (illustrates imperfect agreement between the two tests)";
proc sgplot data=diabetes_screen;
    scatter x=FPG y=HbA1c / group=True_Diabetes transparency=0.6;
    xaxis label="Fasting Plasma Glucose (mg/dL)";
    yaxis label="HbA1c (%)";
run;
title;


/*
 SECTION 3: SENSITIVITY / SPECIFICITY / PPV / NPV ACROSS CANDIDATE CUTOFFS
*/
%macro cutoff_stats(flagvar, cutlabel);
    title "Screening Performance at HbA1c Cutoff = &cutlabel%";
    proc freq data=diabetes_screen;
        tables True_Diabetes*&flagvar / nocol nopercent sensspec;
    run;
    title;
%mend;

%cutoff_stats(Flag_57, 5.7);
%cutoff_stats(Flag_60, 6.0);
%cutoff_stats(Flag_63, 6.3);
%cutoff_stats(Flag_65, 6.5);   /* current ADA cutoff */

%macro build_summary;
    %do i = 1 %to 4;
        %let v = %scan(Flag_57 Flag_60 Flag_63 Flag_65, &i);
        %let c = %scan(5.7 6.0 6.3 6.5, &i, %str( ));

        proc freq data=diabetes_screen noprint;
            tables True_Diabetes*&v / out=_ct&i (drop=percent);
        run;

        data _sum&i;
            length Cutoff 8;
            Cutoff = &c;
            set _ct&i end=last;
            retain TP FP TN FN 0;
            if True_Diabetes=1 and &v=1 then TP=count;
            if True_Diabetes=0 and &v=1 then FP=count;
            if True_Diabetes=0 and &v=0 then TN=count;
            if True_Diabetes=1 and &v=0 then FN=count;
            if last then do;
                Sensitivity = TP / (TP+FN);
                Specificity = TN / (TN+FP);
                PPV         = TP / (TP+FP);
                NPV         = TN / (TN+FN);
                Screen_Pos_per_1000 = ((TP+FP)/(TP+FP+TN+FN))*1000;
                output;
            end;
            keep Cutoff TP FP TN FN Sensitivity Specificity PPV NPV Screen_Pos_per_1000;
        run;
    %end;
    data cutoff_summary;
        set _sum1-_sum4;
    run;
%mend;
%build_summary;

title "Cutoff Comparison Summary: Sensitivity/Specificity/PPV Tradeoff";
proc print data=cutoff_summary noobs;
    format Sensitivity Specificity PPV NPV percent8.1;
    var Cutoff TP FP TN FN Sensitivity Specificity PPV NPV Screen_Pos_per_1000;
run;
title;


/*
 SECTION 4: ROC ANALYSIS
*/
ods graphics on;
title "ROC Curve: HbA1c as a Continuous Predictor of True Diabetes Status";
proc logistic data=diabetes_screen plots(only)=roc;
    model True_Diabetes(event='1') = HbA1c / outroc=roc_data;
run;
title;

data roc_youden;
    set roc_data;
    youden = _SENSIT_ - _1MSPEC_;
run;
proc sort data=roc_youden; by descending youden; run;

title "Top 5 Candidate Cutoffs by Youden's Index (from ROC curve)";
proc print data=roc_youden(obs=5) noobs;
    var _PROB_ _SENSIT_ _1MSPEC_ youden;
    label _PROB_="Predicted Prob. Cutpoint" _SENSIT_="Sensitivity" _1MSPEC_="1-Specificity";
run;
title;


/*
 SECTION 5: AGE/BMI-ADJUSTED MODEL
*/
title "Logistic Model: True Diabetes ~ HbA1c + Age + BMI";
proc logistic data=diabetes_screen plots(only)=roc;
    model True_Diabetes(event='1') = HbA1c Age BMI / outroc=roc_adj;
run;
title;

title "Model Comparison: Unadjusted vs. Age/BMI-Adjusted ROC (AUC)";
proc logistic data=diabetes_screen;
    model True_Diabetes(event='1') = HbA1c Age BMI;
    roc 'Unadjusted' HbA1c;
    roc 'Age+BMI Adjusted' HbA1c Age BMI;
    roccontrast reference('Unadjusted') / estimate e;
run;
title;


/*
 SECTION 6: DISPARITIES ANALYSIS - PREVALENCE AND TEST PERFORMANCE BY
 RACE/ETHNICITY
*/
title "Diabetes Prevalence by Race/Ethnicity";
proc freq data=diabetes_screen;
    tables race_eth*True_Diabetes / nocol nopercent chisq;
run;
title;

title "HbA1c=6.5% Cutoff Performance by Race/Ethnicity (screening equity check)";
proc sort data=diabetes_screen; by race_eth; run;
proc freq data=diabetes_screen;
    by race_eth;
    tables True_Diabetes*Flag_65 / nocol nopercent sensspec;
run;
title;


/*
 SECTION 7: EXPORT SUMMARY REPORT
*/

/* REPLACE WITH YOUR FILE PATH:
*/
ods pdf file="/path/HbA1c_Cutoff_Analysis_Report.pdf" style=journal;
/* 
*/

title "HbA1c Diabetes Screening Cutoff Review - Summary Report";
proc print data=cutoff_summary noobs label;
    format Sensitivity Specificity PPV NPV percent8.1;
    var Cutoff TP FP TN FN Sensitivity Specificity PPV NPV Screen_Pos_per_1000;
    label Cutoff="HbA1c Cutoff (%)"
          Screen_Pos_per_1000="Screen-Positive / 1,000 Screened";
run;
title;

ods pdf close;

proc datasets library=work nolist;
    delete _ct1-_ct4 _sum1-_sum4;
quit;

ODS RTF CLOSE;
