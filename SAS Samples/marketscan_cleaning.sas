/*********Clean Marketscan data + restrict sample*******/

/*set options*/
options realmemsize = max;
options compress = yes;
options obs = max;
options helpbrowser = sas;

libname outaf "/data/markscan_authorized_users/nham/clean_data/adult/full";
libname outcf "/data/markscan_authorized_users/nham/clean_data/child/full";
libname outas "/data/markscan_authorized_users/nham/clean_data/adult/sample";
libname outcs "/data/markscan_authorized_users/nham/clean_data/child/sample";

%let startyear = 2015;
%let endyear = 2018;

/*macro to loop through list and assign each element to a global variable*/
%macro macarray(group, number, list);
  %global n&group.;
  %let n&group.= &number.;
  %put n&group. = &&n&group.;

  %do a=1 %to &&n&group.;
     %global &group.&a.;
     %let &group.&a. = %scan(&list, &a., ' ');
     %put &group.&a. = &&&group.&a. ;
  %end;
%mend macarray;

/*assign file names to global variables*/
%macarray(group = ages, number = 2, list = a c);
%macarray(group = files, number = 6, list = a i s f o d);

/*main macro to clean data*/
%macro clean(obs=, make_ids=, inputs=);

  %do year = &startyear. %to &endyear.; %let yr = %substr(&year.,3,2);
    libname raw&yr. "/data/markscan_authorized/data/commercial/&year.";
  %end;

  /*1. make sample restrictions and create ID files*/
  %if &make_ids = 1 %then %do;
    %do year = &startyear. %to &endyear.; %let yr = %substr(&year.,3,2);

    	/*codes specifying how datasets are named, e.g.: ccaea171 is ccae admissions dataset for year 2017, num=1 because 2017>2019*/
      %if &year. < 2009 %then %do;
        %let num = 2;
      %end;
      %if &year. >= 2009 %then %do;
        %let num = 1;
      %end;

      /* 1.1. Clean enrollment file: keep if enrollment-months is not missing, patient has drug benefit and mental health coverage, region not missing*/
      data enrol_&yr._a enrol_&yr._c;
        set raw&yr..ccaea&yr.&num. (obs = &obs. where = (enrmon ^= . and rx = '1' and mhsacovg = '1' and region ^= '5' ));
        by enrolid;

        /*number of months in HMO/non-HMO plan*/
        nonhmo_mths = 0;
        hmo_mths = 0;

        %do t = 1 %to 12;
          mth&t. = 0; /*flag for hmo*/
          if plntyp&t. in (3,5,6,8,9) then do; /*for epo, pos, ppo, cdhp, hdhp*/
              nonhmo_mths + 1;
              mth&t. = 1; /*flag for non-hmo*/
          end;
          if plntyp&t. in (1,2,4,7) then hmo_mths + 1; /*for basic/major medical, comprehensive, hmo, pos with capitation -> hmo*/
        %end; *end of 12;

        if nonhmo_mths ^= 12 then partial = 1; else partial = 0; /*flag for partial non-hmo beneficiaries*/
        total_mths = nonhmo_mths + hmo_mths; /*total number of months enrolled in any plan*/
        if total_mths ^= 12 then real_partial = 1; else real_partial = 0; /*flag for partial beneficiaries, any plan type*/

  	  if nonhmo_mths > 0; /*keep only ppl with at least 1 non-hmo month*/

        if age >= 18 then output enrol_&yr._a; /*output people >=18yo to adult dataset*/
        else if age < 18 then output enrol_&yr._c; /*output people <18yo to child dataset*/

      run; /*end of enrollment file*/

      /*1.2. Drop benes who ever had a negative, missing or capitated claim in inpatient, outpatient or drug files*/

  	  /*get IDs with negative claims from I file*/
      data drop_ip_&yr. (keep = enrolid);
        set raw&yr..ccaei&yr.&num. (obs = &obs. keep = enrolid hosppay);
        by enrolid;
        if first.enrolid then ev_neg = 0;
        retain ev_neg;
        if  hosppay < 0 then ev_neg = 1;
        if last.enrolid and ev_neg = 1;
      run;

  	  /*get IDs with capitated claims from S file*/
      data drop_sv_&yr. (keep = enrolid);
        set raw&yr..ccaes&yr.&num. (obs = &obs. keep = enrolid cap_svc);
        by enrolid;
        if first.enrolid then do;
            ev_cap = 0;
        end;
        retain ev_cap;
        if cap_svc = "y" then ev_cap = 1;
        if last.enrolid and ev_cap = 1;
      run;

  	  /*get IDs with negative or capitated claims from O file*/
      data drop_op_&yr. (keep = enrolid);
        set raw&yr..ccaeo&yr.&num. (obs = &obs. keep = enrolid pay cap_svc);
        by enrolid;
        if first.enrolid then do;
            ev_neg = 0;
            ev_cap = 0;
        end;
        retain ev_neg ev_cap;
        if  pay < 0 then ev_neg = 1;
        if cap_svc = "y" then ev_cap = 1;
        if last.enrolid and (ev_neg = 1 or ev_cap = 1);
      run;

  	  /*get IDs with negative or capitated claims from D file*/
      data drop_rx_&yr. (keep = enrolid);
        set raw&yr..ccaed&yr.&num. (obs = &obs. keep = enrolid pay cap_svc);
        by enrolid;
        if first.enrolid then do;
            ev_neg = 0;
            ev_cap = 0;
        end;
        retain ev_neg ev_cap;
        if  pay < 0 then ev_neg = 1;
        if cap_svc = "y" then ev_cap = 1;
        if last.enrolid and (ev_neg = 1 or ev_cap = 1);
      run;

  	/*merge with enrollment file to drop IDs with negative/capitated claims*/
      %do v = 1 %to &nages.; %let age = &&ages&v..;
        data out&age.f.id_cl_&yr. (keep = enrolid mth: hmo_mths nonhmo_mths partial real_partial);
          merge enrol_&yr._&age. (in = in_enr) drop_ip_&yr. (in = in_ip) drop_sv_&yr. (in = in_sv) drop_op_&yr. (in = in_op)  drop_rx_&yr. (in = in_rx);
            by enrolid;
          if in_enr and (not in_ip and not in_sv and not in_op and not in_rx);
        run;

        /*reshape ID file from wide to long*/
        data out&age.f.id_cl_&yr._long (keep = enrolid month);
          set %do t = 1 %to 12; out&age.f.id_cl_&yr. (in = in_&t. where = (mth&t. = 1)) %end;;
          month = 0;
          %do b=1 %to 12;
            if in_&b. then month = &b.;
          %end;
        run;

      %end; /*end of age loop*/
    %end; /*end of year loop*/
  %end; /*end of make_ids*/

  /*2. Merge ID files to raw files*/
  %if &inputs = 1 %then %do;

    %do year = &startyear. %to &endyear.; %let yr = %substr(&year.,3,2);

      %if &year. < 2009 %then %do;
        %let num = 2;
      %end;
      %if &year. >= 2009 %then %do;
        %let num = 1;
      %end;

      %do v = 1 %to &nages.; %let age = &&ages&v..;

        proc sort data = out&age.f.id_cl_&yr._long (obs = &obs. keep = enrolid month) out = id_cl_&yr._long;
          by enrolid month;
        run;

      %do b = 1 %to &nfiles.; %let file = &&files&b..;
        %if &file. ^= a %then %do; /*for s f o d i tables*/
            data ccae&file.&yr.&num.;
              set raw&yr..ccae&file.&yr.&num. (obs = &obs.);
                 by enrolid;

                 %if &file. ^= i %then %do; /*for s f o d tables*/
                   month = month(svcdate);
                 %end;
                 %else %if &file. = i %then %do;
                   month = month(admdate);
                 %end;
            run;

            proc sort data =  ccae&file.&yr.&num. out = ccae&file.&yr.&num.;
              by enrolid month;
            run;

            data out&age.f.&file.&yr.;
                merge ccae&file.&yr.&num. (in = in_raw obs = &obs.) id_cl_&yr._long (in = in_ids obs = &obs.);
                  by enrolid month;
                if in_ids and in_raw;
              run;
        %end;

        %else %if &file. = a %then %do;
          data out&age.f.&file.&yr.;
              merge raw&yr..ccae&file.&yr.&num. (in = in_raw obs = &obs.) out&age.f.id_cl_&yr. (in = in_ids obs = &obs. keep = enrolid mth:);
                by enrolid;
              if in_ids;
            run;
        %end;

      %end; /*end of file loop*/

    %end; /*end of age loop*/

  %end; /*end of year loop*/

%end; /*end of inputs*/

%mend clean;

%*execute;
%clean(obs = max, make_ids = 1, inputs = 1);
