***********************************************************************
	PROGRAMMER:      Tracie Shing
  
    DATE:            09/02/2021
----------------------------------------------------------------------
	NAME:         	 CLF_IMPUTATION_MACRO.sas

	DESCRIPTION: 	 Macro for imputation of binary disease status of
                     tooth sites based upon the conditional linear 
					 family (CLF) of distributions for simulating
                     correlated variables


	LANGUAGE/VER:    SAS - Ver 9.4

	HISTORY (PROG):  tshing (09/02/2021)

	PROGRAM NOTES:	
----------------------------------------------------------------------
REQUIRED INPUT:

	seed 		= number to set a seed
	nimp 		= number of imputations
	betaest 	= dataset of beta parameter estimates from marginal
					mean model (ex. output from GEECORR macro)
	betacov 	= dataset of covariance of beta parameters (ex. output 
					from GEECORR macro)
	alphaest 	= dataset of alpha parameter estimates from marginal 
					pairwise correlation model (ex. output from GEECORR macro)
	alphacov 	= dataset of covariance of alpha parameters (ex. output 
					from GEECORR macro),
	idvar 		= cluster(subject)-level id variable
	unitvar		= unit-level id variable (within cluster)
	xvar 		= variables names of marginal mean model covariates
	yvar 		= variable name of outcome
	zvar 		= variable name of pairwise correlation model covariates
	type		= FULL/PARTIAL (default=FULL)
					-- FULL = imputes up to 168 tooth sites, regardless of
						of whether or not the tooth was missing (for pi1/pi2)
					-- PARTIAL = imputes only up to the number of sites on
						non-missing teeth in the mouth (for piA)
	nvar 		= variable name of total cluster size
					(ex. 168 for full mouth cluster)
	xdata 		= dataset of outcome to be imputed and marginal mean covariates, 
					sorted with all observed values first within each 
					cluster and missing values for all values that need to be imputed
					-- contains idvar, unitvar, xvar, and yvar
	zdata 		= dataset of pairwise correlation covariates. should 
					correspond to xdata where records are sorted with all 
					observed values first within each cluster and missing values 
					for all values that need to be imputed.
					ex. if the full cluster size is 168 but only 162 records are observed, 
					the zdata should still be 168C2=14028 records but the matrix is 
					ordered such that the pairwise correlations for the missing 
					values are last
					-- contains idvar, zvar
	cdata 		= dataset of the total cluster size needed and the number of 
					observed values for each cluster(subject)
					-- contains idvar, nvar, nobsvar
	testxdata 	= dataset of test xdata corresponding to the largest cluster size
					-- contains idvar, unitvar, and xvar
	testzdata 	= dataset of test zdata corresponding to the largest cluster size, 
					should be the same for all clusters
					-- contains zvar
	testcdata 	= dataset of test cdatda corresponding to the largest cluster size 
					-- contains idvar, nvar
	outdata 	= name output dataset containing nimp-imputations 

OPTIONAL INPUT:

	processing 	= YES/NO (default=NO)
					if YES, then creates output dataset containing idvar unitvar xvar yvar and 
					nimp-columns of imputed/observed values in variables &yvar._imp1-&yvar._imp&nimp.
					(one column=one imputation)
					-- observed values are ordered first and imputed values are ordered last (matching xdata)						
						*for any impution during which observations could not be imputed due to
							range violations or the variance is not positive definite, all values
							for the imputation are set to missing to indicate a failed imputation
					-- conditional mean for imputed values in variables (yvar)_clf1-(yvar)_clf(nimp)
						*values of . correspond to observations that were not imputed
						*for any impution during which observations could not be imputed due to
							range violations or the variance is not positive definite, all values
							for the imputation are set to missing to indicate a failed imputation
					-- error flags to indicate if the conditional mean for imputed values are
						outside the range [0,1] in variables (yvar)_clferr1-(yvar)_clferr(nimp)
						*values of . correspond to observations that were not imputed
					if NO, then creates dataset containing idvar unitvar xvar yvar and 
					nimp-columns of imputed/observed values in variables &yvar._imp1-&yvar._imp&nimp.
					(one column=one imputation)
					-- observed values are ordered first and imputed values are ordered last (matching xdata)
						*values of -1 correspond observations that could not be imputed due to
							range violations in the cluster or the variance of the cluster is not positive definite
					-- conditional mean for imputed values in variables (yvar)_clf1-(yvar)_clf(nimp)
						*values of -1 correspond to observations that were not imputed
	impbounds = YES/NO (default=NO)
					if YES, then values where the conditional mean is <0 or >1 are set to 0 and 1 respectively
					if NO, then values are 0 and 1 with site-level probability 

	OUTPUT: Dataset called &outdata containing
			-- &idvar &unitvar &xvar &yvar
			-- &yvar._imp1-&yvar._imp&nimp. (nimp-columns of imputed/observed values in variables)
			-- &yvar._clf1-&yvar._clf&nimp. (nimp-columns of conditional means from CLF)
***********************************************************************;

%macro clf_imp (seed=,nimp=, 
				betaest=,betacov=,alphaest=,alphacov=,
				idvar=,unitvar=,xvar=,yvar=,zvar=,nvar=,
				xdata=,zdata=,cdata=,
				testxdata=,testzdata=,testcdata=,outdata=,
				processing=,impbounds=);
proc IML;
	*load modules;
	************************************************************;
	*helper functions;
	************************************************************;
	start BEGINEND(first, last, n);
		last = cusum(n);
		first = last - n + 1;
		finish BEGINEND;
	/* returns: <-
		first = a vector of first row # of cluster
		last = a vector of last row # of cluster */
	************************************************************;
	start binvar(p);
		return(p#(1-p));
	finish;
	/* returns: <- the binary variance function */
	************************************************************;
	start logit(p);
	  return (log(p/(1-p)));
	finish;
	/* returns: <- logit(p) */
	************************************************************;
	start alogit(t);
	  return (1/(1+exp(-t)));
	finish;
	/* returns: <- antilogit(t) */
	************************************************************;
	start ch2(n);
	  return ((n#(n-1))/2);
	finish;
	/* returns: <- n choose 2 */
	************************************************************;
	start ch2inv(m);
	  n = (1+sqrt(1+8*m))/2;
	  return (round(n));
	finish;
	/* returns: <- n such that m == (n choose 2) */
	************************************************************;
	start cor2var (r, u);
	/*
	returns variance matrix of binary variables with mean
	vector u[] and corr matrix r[,].
	*/
	  s = sqrt(u#(1-u));
	  return ( s # r # t(s) );
	finish;
	************************************************************;
	start var2cor (v);
	/*
	returns corr matrix
	*/
	  s = 1/sqrt(vecdiag(v));
	  return ( s # v # t(s) );
	finish;
	************************************************************;
	start allreg (v);
	/*
	v[1:n, 1:n] is the input covariance matrix of Y[1:n].
	v[,] is assumed +ve definite symmetric, not checked.
	Returns  b[1:n,1:n] such that b[1:t-1, t] are the coefficients
	of y[1:t-1] in the cond. mean of Y[t], for t=2:n.
	Diagonals and lower half of b[,] are zero-filled.
	*/

	  n  = nrow(v);
	  u  = root(v);   * upper Choleski root;
	  u0 = u;
	  d  = do(1, n * n, n + 1);    * diagonals;
	  u0[d] = j(n, 1, 0);    * set diagonals = 0;
	  b  = trisolv(1, u, u0);

	  return (b);
	finish;
	************************************************************;
	*functions for CLF checking;
	***********************************************************;
	start is_pos_def (A);
	    * A = symmetric matrix;
	    return (min(eigval(A)) > 0);
	finish;
	/* returns 1 if A is positive definite
	           0 otherwise */
	************************************************************;
	start chkbinc  (r, mu, i, j);
	/*
	Find all pairwise range violations.

	Input: r, mu
	Output: i, j, return value

	r[1:n, 1:n] =  the correlation matrix (only upper half is checked)
	mu[1:n] =  the mean vector

	To check that the pairwise correlations are compatible with the
	means:

	  err =  chkbinc (r, mu, i, j);

	err <- number of range violation detected (0 means no violations).
	(i,j) <- (row,column) locations of all errors found.
	If no violations are found, i and j are set to "." (missing).
	i=row, j=column, i<j only (upper half).
	*/
	  n = nrow(mu);
	  v = cor2var (r, mu);
	  nc2 = (n # (n-1)) / 2;
	  ierror = j(nc2, 1, 0);
	  jerror = j(nc2, 1, 0);
	  kk = 0;        * error counter;
	  do ii = 1 to n-1;
	    ui = mu[ii];
	    do jj = ii+1 to n;
	      uj  = mu[jj];
	      uij = ui*uj + v[ii,jj];
	      ok  = (uij <= min(ui, uj)) & (uij >= max(0, ui+uj-1));
	      if (^ok) then do;
	        kk = kk + 1;
	        ierror[kk] = ii;
	        jerror[kk] = jj;
	      end;
	    end;
	  end;

	  if (kk > 0) then do;
	    i = ierror[1:kk];
	    j = jerror[1:kk];
	  end;
	  else do; i=.; j=.; end;

	  return (kk);
	finish;
	************************************************************;
	start get_bnds (/*check_u, check_u2, b2,*/ nu_min, nu_max, i, u, b);
	/*
	output: nu_min <- min of the i-th cond. mean.
	output: nu_max <- max of the i-th cond. mean.

	input: u[1:n] = mean vector,
	input: b[1:(i-1)] = reg coeffs
	range of i is 2:n
	*/
	  y = (b > 0);
	  nu_max = u[i] + (y-u[1:(i-1)])` * b;
	  y = 1 - y;
	  nu_min = u[i] + (y-u[1:(i-1)])` * b;

	  /*for checking uncomment:
	  b2= b;
	  check_u = u[i];
	  check_u2 = u[1:(i-1)];*/
	finish;
	************************************************************;
	start blrchk1 (u, b);
	/* check for CLF compatibiltiy */
	/* u[1:n] := mean, b[1:n, 1:n] := matrix computed by allreg() */
	/* no printed output */
	/* returns unit number if there is an error (cond. means out of range) */
	/* returns 0 otherwise  */

	  n = nrow(b);
	  rc = 0;       * return code;
	  do i = 3 to n  while (rc=0); * n=2 is guarnteed CLF compatible;
	    run get_bnds (nu_min, nu_max, i, u, b[1:(i-1), i] );
		if ((nu_max > 1) | (nu_min < 0)) then rc = i; *this returns the unit at which the cond. mean is out of range;
	  end;
	  return (rc);
	finish;
	************************************************************;
	*CLF module;
	***********************************************************;
	start mbsclfimp1y (yimp, ciimp, u, B, y, ntot);
	/*
	Modified Multivariate Bernoulli simulation from the CLF for imputation

	  yimp = mbsclfimp1y (u, B, y, nobs, ntot);

		u = mean vector
		B = matrix computed by B = allreg(V)
		y = vector with observed values and values to be imputed
		nobs = number of observed values
		ntot = total size of cluster
		
		returns:
		yimp <- one simulated column vector, contains both observed and imputed values
		ciimp <- conditional means and error flags
					i=number of error,
					ci=conditional mean, 
					1 for Conditional mean outside [0,1], 0 otherwise
	*/
	  nend = ntot;
	  nstart=2;
	  yimp = y;
	  ciimp = J(ntot,1,-1);

	  if yimp[1] = -1 then do;
	  	yimp[1] = rand("Bernoulli",u[1]);
	  end;

	  r = (yimp - u)`;                            * residuals (row);

	do i = nstart to nend;
		i1 = i - 1;

		if yimp[i] = -1 then do;
			ci = u[i] + (r[1, 1:i1]) * b[1:i1, i];     * cond. mean;

			if (ci < 0 | ci > 1) then do;
				*print "ERROR: Conditional mean outside [0,1] in mbsclfimp1y:" i ci;
				%if %upcase(&impbounds)=YES %then %do;
					if ci<0 then yimp[i] = 0;
					if ci>1 then yimp[i] = 1;
				%end;
				%if %upcase(&impbounds)^=YES %then %do;
					yimp[i] = rand("Bernoulli",u[i]);
				%end;
				r[i] = yimp[i] - u[i]; * residual;
			end;
			else do;
				yimp[i] = rand("Bernoulli",ci);
				r[i] = yimp[i] - u[i]; * residual;
			end;
			ciimp[i] = ci;
		end;
	end;

	finish;
	***********************************************************;
	*end load modules;
	***********************************************************;

	call randseed(&seed);
	call streaminit(&seed);

	*read in beta, alpha, and empirical variance;
	USE &betaest;
	READ ALL VAR {estimate} INTO betahat;
	close &betaest;

	USE &betacov;
	READ ALL VAR {&xvar} INTO betacov;
	close &betacov;

	USE &alphaest;
	READ ALL VAR {estimate} INTO alphahat;
	close &alphaest;

	USE &alphacov;
	READ ALL VAR {&zvar} INTO alphacov;
	close &alphacov;


	*Step 1. Random draw of beta and alpha parameters from MVN
	(each imputation replicate will have a different mean and correlation);

	betastar = J(&nimp,nrow(betahat),0); *initialize matrix for drawn beta for all imputations;
	alphastar = J(&nimp,nrow(alphahat),0); *initialize matrix for drawn alpha for all imputations;


	*Step 2. Test CLF compatability for hypothetical 168 tooth sites 
	for each age group;

	*read in data for hypothetical subjects;
	use &testcdata;
	read all var{ntot} into ntot; *size of cluster;
	read all var{seqn} into clusterid; *list of clusterids;
	close &testcdata;

	USE &testxdata;
	READ ALL VAR {&xvar} INTO testx; *xdata;
	close &testxdata;

	USE &testzdata;
	READ ALL VAR {&zvar} INTO testz; *zdata - same for all participants;
	close &testzdata;

	run BEGINEND(clusterstart, clusterend, ntot); *create starting and end row number for each cluster (for iterating);

	*random draw of beta, alpha for all imputations;
	do impnum = 1 to &nimp;
	print "IMPUTATION #:" impnum;

		counter=0; *redraw counter;
		do j = 1 to 50 until(compat_min=1);

			counter=counter+1;
			print "RANDOM DRAW NUMBER:" counter;
			
			*random draw from MVN;
			betastar_tmp = RandNormal(1,betahat,betacov);
			alphastar_tmp = RandNormal(1,alphahat,alphacov);
			print betastar_tmp;
			print alphastar_tmp;

			*testing CLF compatibility;
			compat_min=0; *tracker for minimum compatibility (no range violations and var is pos def);
			do i = 1 to nrow(ntot);
				c = clusterid[i];
				print "Cluster ID:" c;

				testx_c = testx[clusterstart[i]:clusterend[i],];
				testz_c = testz[1:14028,];

				*****Specify mu and R based on drawn betastar,alphastar;
				xb_c = testx_c*t(betastar_tmp);
				mu_c = alogit(xb_c); *mu should be vector of 168 x 1;

				rho_c = testz_c*t(alphastar_tmp); *rho should be vector of 14028 x 1;

				*****Create correlation matrix from correlation vector rho_c;
				R_dim = ch2inv(nrow(rho_c));                * dimension of full matrix - should be 168 x 168;
				R_c = J(R_dim,R_dim,0);          			* initialize zero matrix;
				upperTri = loc(row(R_c) < col(R_c)); 		* upper tri indices in row major order;
				R_c[upperTri] = rho_c`;       				* copy elements;
				R_c = R_c + R_c`;            				* make symmetric;
				diag = loc(row(R_c) = col(R_c));     		* diagonal elements;
				R_c[diag] = 1;           					* put 1 on diagonal;
				free R_dim upperTri diag;

				*****Create covariance matrix and intermediate variance matrix;
				V_c = cor2var (R_c, mu_c);    * v <- cov matrix;

				*****Checks;
				*check correlation ranges;
				_i=.; _j=.;     * SAS/IML requires defined variables;
				err1 =  chkbinc (R_c, mu_c, _i, _j);
				if (err1=1) then print "Pairwise range violations at" _i _j;
				if (err1=1) then CONTINUE;

				*check if variance is positive definite;
				var_pos_def = is_pos_def(V_c);
				if (var_pos_def^=1) then print "Variance is NOT positive definite";
				if (var_pos_def^=1) then CONTINUE;

				corr_pos_def = is_pos_def(R_c);
				if (corr_pos_def^=1) then print "Correlation is NOT positive definite";
				if (corr_pos_def^=1) then CONTINUE;

				compat_min=1;
			end;
		end;

		betastar[impnum,] = betastar_tmp;
		alphastar[impnum,] = alphastar_tmp;
	end;
	*clean up matrices for actual CLF imputation;
	free c testx_c testz_c xb_c mu_c rho_c R_c V_c B_c
		var_pos_def corr_pos_def err1 err2
		i j _i _j counter compat
		impnum betastar_tmp alphastar_tmp
		ntot clusterid testx testx_id testz clusterstart clusterend;

	print "PARAMETER VALUES FOR &nimp IMPUTATIONS", "(each row is 1 imputation)", betastar, alphastar;

	*Step 3. CLF data augmentation;

	*read in data for subjects;
	use &cdata;
	read all var{&nvar} into ntot; *total size of cluster;
	read all var{&idvar} into clusterid; *list of clusterids;
	close &cdata;

	USE &xdata;
	READ ALL VAR {&xvar} INTO xdata; *xdata;
	READ ALL VAR {&yvar} INTO ydata; *ydata;
	close &xdata;
		*set missing values of ydata to -1;
		idx = loc( ydata = . );
		ydata[ idx ] = -1;
		free idx;

	USE &zdata;
	READ ALL VAR {&zvar} INTO zdata;
	close &zdata;

	run BEGINEND(clusterstart, clusterend, ntot); *create starting and end row number for xydata of each cluster (for iterating);
	run BEGINEND(zclusterstart, zclusterend, ch2(ntot)); *create starting and end row number for zdata of each cluster (for iterating);

	*perform data augmentation for all imputation replicates;
	yimp = J(ntot[+,],&nimp,-1); *initialize matrix for all imputation output;
	ciimp = J(ntot[+,],&nimp,-1); *initialize matrix for all imputation output;

	do impnum = 1 to &nimp;
	print "IMPUTATION NUMBER:" impnum;
	betastar_imp = betastar[impnum,];
	alphastar_imp = alphastar[impnum,];
	print betastar_imp;
	print alphastar_imp;

		do i = 1 to nrow(clusterid);
			c = clusterid[i];

			*****gather data for cluster;
			ntot_c = ntot[i];

			y_c = ydata[clusterstart[i]:clusterend[i],];

				x_c = xdata[clusterstart[i]:clusterend[i],];
				z_c = zdata[zclusterstart[i]:zclusterend[i],];

				*****Specify mu and R based on drawn betastar,alphastar;
				xb_c = x_c*t(betastar_imp);
				mu_c = alogit(xb_c); *mu should be vector of 168 x 1;

				rho_c = z_c*t(alphastar_imp); *rho should be vector of 14028 x 1;

				*****Create correlation matrix from correlation vector rho_c;
				R_dim = ch2inv(nrow(rho_c));                * dimension of full matrix - should be 168 x 168;
				R_c = J(R_dim,R_dim,0);          			* initialize zero matrix;
				upperTri = loc(row(R_c) < col(R_c)); 		* upper tri indices in row major order;
				R_c[upperTri] = rho_c`;       				* copy elements;
				R_c = R_c + R_c`;            				* make symmetric;
				diag = loc(row(R_c) = col(R_c));     		* diagonal elements;
				R_c[diag] = 1;           					* put 1 on diagonal;
				free R_dim upperTri diag;

				*****Create covariance matrix;
				V_c = cor2var (R_c, mu_c);    * v <- cov matrix;

				*****Checks;
				*check correlation ranges;
				_i=.; _j=.;     * SAS/IML requires defined variables;
				err1 =  chkbinc (R_c, mu_c, _i, _j);

				*check if variance is positive definite;
				var_pos_def = is_pos_def(V_c);

				corr_pos_def = is_pos_def(R_c);
				
				if err1=1 | var_pos_def ^=1 | corr_pos_def ^=1 then do;
					print "Errors for Cluster ID:" c;
					if (err1) then print "Pairwise range violations at" _i _j;
					if (var_pos_def ^=1) then print "Variance is NOT positive definite";
					if (corr_pos_def ^=1) then print "Correlation is NOT positive definite";

					yimp[clusterstart[i]:clusterend[i],impnum] = y_c;
					ciimp[clusterstart[i]:clusterend[i],impnum] = J(nrow(y_c),1,-1);
				end;
				else do;
					*****Create intermediate variance matrix - only if variance is positive definite;
					B_c = allreg(V_c); *intermediate variance matrix needed for mbsclf;

					*****simulation (imputation);
					run mbsclfimp1y(yimp_c, ciimp_c, mu_c, B_c, y_c, ntot_c);

					yimp[clusterstart[i]:clusterend[i],impnum] = yimp_c;
					ciimp[clusterstart[i]:clusterend[i],impnum] = ciimp_c;
				end;
		end;

		free c ntot_c x_c y_c z_c xb_c mu_c rho_c R_c B_c V_c yimp_c ciimp_c
				i _i _j err1 err2 corr_pos_def var_pos_def;
	end;
	free betastar_imp alphastar_imp impnum;

	**** output;
	USE &xdata;
	READ ALL VAR {&idvar} INTO xydata_id; *xdata including the clusterid;
	READ ALL VAR {&unitvar &yvar} INTO xydata; *xdata including the clusterid;
	close &xdata;

	xyimp = yimp || ciimp;

	impnames = "&yvar._imp1":"&yvar._imp&nimp";
	cinames = "&yvar._clf1":"&yvar._clf&nimp";
	names = {&idvar &unitvar &yvar}||impnames||cinames;

	create outdata_tmp from xydata_id xydata xyimp [colname=names];
		append from xydata_id xydata xyimp;
	close outdata_tmp;
quit;

%if %upcase(&processing)=YES %then %do;
	*programming to check/process output data;
	title "Checking completeness of imputations";
	title2 "if any imputations contain -1 then imputation failed";
	proc tabulate data=outdata_tmp missing;
		class &yvar._imp1-&yvar._imp&nimp;
		tables &yvar._imp1-&yvar._imp&nimp, n pctn;
	run;
	title2;
	title;

	***set all values for failed imputations to missing;
	proc means data=outdata_tmp noprint;
		var &yvar._imp1-&yvar._imp&nimp;
		output out=check min=&yvar._imp1-&yvar._imp&nimp;
	run;
	proc transpose data=check out=checkt prefix=MIN;
		var &yvar._imp1-&yvar._imp&nimp;
	run;
	data checkt;
		set checkt;
		name2 = tranwrd(_name_,'imp','clf');
	run;
	proc sql noprint;
		select strip(_name_) into :impfailvar separated by ' '
		from checkt
		where min1=-1;
		select strip(name2) into :impfailclfvar separated by ' '
		from checkt
		where min1=-1;
	quit;

	data &outdata;
		set outdata_tmp;
		
		%if &impfailvar ^= %then %do;
		*set all values for failed imputations to missing;
		array imperr {*} &impfailvar;
		array imperrclf {*} &impfailclfvar;

		do i=1 to dim(imperr);
			imperr[i]=.;
			imperrclf[i]=.;
		end;
		%end;

		*set conditional mean to missing if conditional mean is -1 
			(i.e. the value was not imputed)
		*create flag to see which values were imputed to 0 or 1
			because conditional mean<0 or >1;
		array ci{*} &yvar._clf1-&yvar._clf&nimp;
		array cierr{*} &yvar._clferr1-&yvar._clferr&nimp;
		do j=1 to dim(ci);
			if ci[j]=-1 then ci[j]=.;
			else do;
				if ci[j]<0 or ci[j]>1 then cierr[j]=1; else cierr[j]=0;
			end;
		end;
		drop i j;
	run;

	title "Checking recoding of failed imputations";
	proc tabulate data=&outdata missing;
		class &yvar._imp1-&yvar._imp&nimp;
		tables &yvar._imp1-&yvar._imp&nimp, n pctn;
	run;
	title;

	title "Checking percent of records where conditional means outside of [0,1]";
	proc tabulate data=&outdata missing;
		class &yvar._clferr1-&yvar._clferr&nimp;
		tables &yvar._clferr1-&yvar._clferr&nimp, n pctn;
	run;
	title;

	title "Checking conditional means distributions";
	proc means data=&outdata n mean std p1 p5 p10 p25 median p75 p90 p95 p99;
		var &yvar._clf1-&yvar._clf&nimp;
	run;
	title;

	*clean up tmp datasets;
	proc datasets library=work noprint;
		delete check checkt outdata_tmp;
	run;
	quit;
%end;
%else %do;
	data &outdata;
		set outdata_tmp;
	run;

	*clean up tmp datasets;
	proc datasets library=work noprint;
		delete outdata_tmp;
	run;
	quit;

%end;
%mend;
