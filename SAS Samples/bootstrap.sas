options realmemsize = max;
options compress = yes;
options obs = max;

/* estimate differences in prevalence, incidence, persistence between ACOs and non-ACOs for 500 bootstrap samples */
%macro bootstrap(obs=, rep=, build=, bootsamp=, matrix=, ratio=, bootse=);

libname data "../../decomp_incidence_cure/temp";
libname match "/n/data1/project/uhsp_tgm/tbn8/analysis/decomp_matching/output";
libname temp "../temp";
libname out "../output";

%include "/n/data1/project/uhsp_tgm/tbn8/lib/macros/macarray.sas";

%macro all_hccs;
  HCC1    HCC2    HCC6    HCC8    HCC9    HCC10   HCC11   HCC12
  HCC17   HCC18   HCC19   HCC21   HCC22   HCC23   HCC27   HCC28
  HCC29   HCC33   HCC34   HCC35   HCC39   HCC40   HCC46   HCC47
  HCC48                   HCC54   HCC55   HCC57   HCC58   HCC70
  HCC71   HCC72   HCC73   HCC74   HCC75   HCC76   HCC77   HCC78
  HCC79   HCC80   HCC82   HCC83   HCC84   HCC85   HCC86   HCC87
  HCC88   HCC96   HCC99   HCC100  HCC103  HCC104  HCC106  HCC107
  HCC108  HCC110  HCC111  HCC112  HCC114  HCC115  HCC122  HCC124
  HCC134  HCC135  HCC136  HCC137
  HCC157  HCC158                  HCC161  HCC162  HCC166  HCC167
  HCC169  HCC170  HCC173  HCC176  HCC186  HCC188  HCC189
%mend all_hccs;


%macarray(group= hccs, number = 79, list = %all_hccs);

%let startyr = 2017;
%let endyr = 2018;

%if &build. = 1 %then %do;
  data temp.original_matched_sample;
    merge match.matched_ids_&endyr.(in=in_id) data.hcc_flags_all_&startyr._&endyr.(in=in_flag);
    by bene_id;
    if in_flag and in_id;

    array _xxx_ _numeric_;
    do i=1 to dim( _xxx_);
      if missing(_xxx_[i]) then _xxx_[i]=0;
    end;

    %do i = 1 %to &nhccs.;
      %let hcc = &&hccs&i..;
      rename &hcc._neither=neither_&hcc.;
      rename &hcc._yr&startyr.=first_&hcc.;
      rename &hcc._yr&endyr.=second_&hcc.;
      rename &hcc._both=both_&hcc.;
    %end;
  run;

%end; /*end of build*/

/* make bootstrap samples */
%if &bootsamp. = 1 %then %do;

  proc sort data=temp.original_matched_sample out=original_sample(obs=&obs.); by in_aco; run;
  proc surveyselect data=original_sample(obs=&obs.) out=temp.bootsample
      seed = 1347 method = urs
      samprate = 1 outhits rep = &rep.;
    strata in_aco;
  run;

  proc sql noprint;
    select count (*) into :nonaco from temp.bootsample where in_aco=0 & replicate=2;
    select count (*) into :aco from temp.bootsample where in_aco=1 & replicate=2;
  quit;
  %put &=nonaco;
  %put &=aco;

  proc sort data=temp.bootsample tagsort; by replicate; run;
%end;

/* count benes by HCC */
%if &matrix. = 1 %then %do;
  ods graphics off;  ods exclude all;  ods noresults;

  proc means data = temp.bootsample stackodsoutput sum;
    var neither: first: second: both:;
    ods output summary = matrix;
    class in_aco;
    by replicate;
  run;
  ods graphics on;  ods exclude none;  ods results;

  data matrix(drop=variable);
    set matrix(drop=_control_ nobs);
    type=scan(variable,1,'_');
    hcc=scan(variable,2,'_');
  run;

  proc sort data=matrix; by replicate in_aco hcc; run;
  proc transpose data=matrix out=temp.matrix(drop=_name_) prefix=sum_;
    by replicate in_aco hcc;
    id type;
    var sum;
  run;

%end; /*end of matrix*/

/* calculate prevalence/incidence/persistence */
%if &ratio. = 1 %then %do;
  /*get number of observations from merged dataset*/
  data _null_;
    if 0 then set temp.original_matched_sample nobs = n;
     call symputx ('totpop',n/2);
    stop;
  run;
  %put &=totpop.;

  data temp.ratios;
    set temp.matrix;

    prev&startyr. = (sum_first + sum_both)/%eval(&totpop.);
    prev&endyr. = (sum_second + sum_both)/%eval(&totpop.);
    incidence = sum_second/(sum_second + sum_neither);
    persistence = sum_both/(sum_first + sum_both);
    ss_prev = incidence/(1 + incidence - persistence);

  run;

%end; /*end of ratio*/

/* compute SEs and CIs */
%if &bootse.=1 %then %do;
  data temp.ratio_diffs;
    merge temp.ratios (drop=sum: prev&startyr. where=(in_aco=1) rename=(prev&endyr.=prev1 incidence=inc1 persistence=pers1 ss_prev=ss_prev1))
      temp.ratios (drop=sum: prev&startyr. where=(in_aco=0) rename=(prev&endyr.=prev0 incidence=inc0 persistence=pers0 ss_prev=ss_prev0));
    by replicate hcc;
    diff_prev = prev1 - prev0;
    diff_inc = inc1 - inc0;
    diff_pers = pers1 - pers0;
    diff_ss_prev = ss_prev1 - ss_prev0;
  run;

  proc sort data=temp.ratio_diffs; by hcc; run;

  ods graphics off;  ods exclude all;  ods noresults;
  proc means
    data=temp.ratio_diffs mean std;
    var diff:;
    by hcc;
    ods output summary=bootse(drop=vname:);
  run;

  proc univariate data=temp.ratio_diffs noprint;
    var diff_prev diff_inc diff_pers diff_ss_prev;
    output out=bootci pctlpts=2.5 97.5
      pctlpre=diff_prev_ diff_inc_ diff_pers_ diff_ss_prev_
      pctlname=lci uci;
    by hcc;
  run;

  data bootse;
    merge bootse(keep=hcc diff_inc:) bootci(keep=hcc diff_inc:)
          bootse(keep=hcc diff_pers:) bootci(keep=hcc diff_pers:)
          bootse(keep=hcc diff_prev:) bootci(keep=hcc diff_prev:)
          bootse(keep=hcc diff_ss_prev:) bootci(keep=hcc diff_ss_prev:);
    by hcc;
    hccnum = input(substr(hcc,4),3.);
  run;

  proc sort data=bootse out=out.bootse(drop=hccnum); by hccnum; run;

  proc export data = out.bootse
              dbms = xlsx
              outfile = "../output/diff_bootse.xlsx"
              replace;
            run;
%end;

%mend bootstrap;

%bootstrap(obs=max, rep=500, build=0, bootsamp=0, matrix=0, ratio=0, bootse=0);
