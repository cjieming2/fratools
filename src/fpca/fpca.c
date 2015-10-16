#include <stdio.h>
#include <nicklib.h>
#include <stdlib.h>
#include "eigsubs.h"
#include <unistd.h>
        
int NSAMPLES, nSNP;

int main(int argc, char **argv)
{
    int k, n, m, nn, col, i, j, extensionLength, val1, val2, N, x, y, rowvalid, colvalid, strLen, pafFile, tgFile;
    int printCovarianceMatrix = 0;
    int rowCenter = 1;
    int pcNo = 20;
    int individualNormalization = 0;
    int populationNormalization = 0;
    char Xchar;
    char *str;
    char **samples;
    char **snps;
    double *X, *XTX, *syyArray, rowsum, rowmean, rowmeanbayes, colsum, colmean, colmeanbayes, tempdouble, sxx, syy, sxy;
    double *eval, *evec, sum;
    FILE *fp, *fpcor, *fpout, *fpeval, *fpcov;
    char *INFILE = NULL;
    char *PCFILE = NULL;
    char *EVALFILE = NULL;
    char *COVFILE = NULL;
    char *CORFILE = NULL;
    
    if(argc==1)
    {
        printf("usage: fpca [options] <paf-file|tg-file>\n");
        printf("\n");
        printf("       -c       columnwise centering (rowwise centering by default)\n");
        printf("       -i       normalization by rate of genetic drift : sqrt(p*(1-p)) - applicable for individuals\n");
        printf("       -p       normalization by : sqrt(pbar*(1-pbar)) - applicable for population\n");
        printf("                default normalization is a centering of the data\n");
        printf("       -v       print out covariance matrix\n");
        printf("       -e       number of principal components to print (default 20)\n");
        printf("       paf-file population allele frequency file\n");
        printf("       tg-file  SNPs x Samples genotype file\n");
        printf("\n");
        printf("       example: fpca -p pscalare.paf\n");
        printf("                fpca -i pscalare.tg\n");
        printf("\n");
        exit(1);
    }
    	
    /* process flags */
    while((i = getopt(argc,argv,"ivpe:")) != -1)
    {
        switch(i)
        {
            case 'v':          
                printCovarianceMatrix = 1;
                break; 
            case 'c':
            	/*DISABLED!!!*/ 
                rowCenter = 1;
                break;
            case 'i':          
                individualNormalization = 1;
                break;              
            case 'p':          
                populationNormalization = 1;
                break; 
            case 'e':
                pcNo = atoi(optarg); 
                break;
            case '?':
            	fprintf(stderr, "Unrecognized option: -%c\n", optopt);
            	exit(1);
        }
    }

	if (optind != argc-1)
	{
		fprintf(stderr, "1 non-option argument expected: paf-file or tg-file\n");
		exit(1);
	}

	INFILE = argv[optind];
	strLen = strlen(INFILE);
	    
    if (strLen >= 4 &&
    	(INFILE[strLen-3] == '.' &&
    	 INFILE[strLen-2] == 't' && 
    	 INFILE[strLen-1] == 'g'))
    {
    	tgFile = 1;
    	pafFile = 0; 
    	extensionLength = 2;
    }
    else if (strLen >= 5 &&
        	 (INFILE[strLen-4] == '.' &&
        	  INFILE[strLen-3] == 'p' &&
        	  INFILE[strLen-2] == 'a' && 
        	  INFILE[strLen-1] == 'f'))
    {   	    	   	    	
    	tgFile = 0;
    	pafFile = 1; 
		extensionLength = 3;
    }
    else
    {
    	fprintf(stderr,"%s not a tgFile or pafFile\n", INFILE);
        exit(1);
    }
	
	if (tgFile)
	{
		if((str = (char *) malloc((strLen+2)*sizeof(*str))) == NULL)
	    { fprintf(stderr,"CM\n");  exit(1); }
	}
	else
	{
		if((str = (char *) malloc((strLen+1)*sizeof(*str))) == NULL)
	    { fprintf(stderr,"CM\n");  exit(1); }		
	}
    	
    strncpy(str, INFILE, strLen-extensionLength);
    str[strLen-extensionLength] = '\0';    		
	PCFILE = strcat(str, "pca"); 	    	
      
	/* open output files */
    if( (fpout = fopen(PCFILE, "w")) == NULL)
    {
        fprintf(stderr,"Could not open pca file %s\n", PCFILE);  exit(1);
    }
   
	if (tgFile)
	{
		if((str = (char *) malloc((strLen+3)*sizeof(*str))) == NULL)
    	{ fprintf(stderr,"CM\n");  exit(1); }
    }
    else
    {
		if((str = (char *) malloc((strLen+2)*sizeof(*str))) == NULL)
    	{ fprintf(stderr,"CM\n");  exit(1); }
    }
	
    strncpy(str, INFILE, strLen-extensionLength);
    str[strLen-extensionLength] = '\0';
	EVALFILE = strcat(str, "eval");
	
	if( (fpeval = fopen(EVALFILE, "w")) == NULL)
    {
        fprintf(stderr,"Could not open eval file %s\n", EVALFILE);  exit(1);
    }

    strncpy(str, INFILE, strLen-extensionLength);
    str[strLen-extensionLength] = '\0';
	CORFILE = strcat(str, "cor"); 	    	

	if(printCovarianceMatrix)
	{
    	if((str = (char *) malloc((strLen+1)*sizeof(*str))) == NULL)
	    { fprintf(stderr,"CM\n");  exit(1); }
	    
	    strncpy(str, INFILE, strLen-extensionLength);
        str[strLen-extensionLength] = '\0';
    	COVFILE = strcat(str, "cov");
	   	if( (fpcov = fopen(COVFILE, "w")) == NULL)
	    {
	        fprintf(stderr,"Could not open covariance file %s\n", COVFILE);  exit(1);
	    }
	}
    
	/* allocate string length */    
    if((str = (char *) malloc(100*sizeof(*str))) == NULL)
    { 
    	fprintf(stderr,"CM\n");  
    	exit(1); 
    }
    
    if( (fp = fopen(INFILE, "r")) == NULL)
    {
        fprintf(stderr,"Could not open input file %s\n", INFILE);  
        exit(1);
    }
    
    if((samples = (char **) malloc(sizeof(*samples))) == NULL)
    { 
    	fprintf(stderr, "CM\n");  
    	exit(1); 
   	}

    if((snps = (char **) malloc(sizeof(*snps))) == NULL)
    { 
    	fprintf(stderr, "CM\n");  
    	exit(1); 
   	}
   	   	
   	/*count number of SNPs and samples*/
   	fprintf(stderr, "Examining matrix\n");
   	
    strLen = 0;
    n = 0;
    m = 0;
    while(1)
    {   	
        fscanf(fp,"%c",&Xchar);
        
        if(feof(fp))
        {
        	break;
        }
        else if (m==0)
        {
	        if(Xchar=='\t')
	        {
	            str[strLen] = '\0';
	            strLen = 0;
	            
	            /*ignores snp-id or sample-id*/
	            if(n!=0)
	            {
	            	if((samples = (char **) realloc(samples, (n+1)*sizeof(*samples))) == NULL)
				    { 
				    	fprintf(stderr, "CM\n");  
				    	exit(1); 
				   	}			    
				   
	            	samples[n-1] = strdup(str);
	            }
	            n++;
	        }
	        else if (Xchar=='\n')
	        {
	            str[strLen] = '\0';
	            strLen = 0;
	            
            	if((samples = (char **) realloc(samples, (n+1)*sizeof(*samples))) == NULL)
			    { 
			    	fprintf(stderr, "CM\n");  
			    	exit(1); 
			   	}			    
			   
            	samples[n-1] = strdup(str);
                		            
	            n++;
	            m++;
	            continue;
	        }
	        else
	        {
	            str[strLen++] = Xchar;
	        }
	    }
	    else
	    {
	        if(col==0 && Xchar=='\t')
	        {
	            str[strLen] = '\0';
            
            	if((snps = (char **) realloc(snps, (m+1)*sizeof(*snps))) == NULL)
			    { 
			    	fprintf(stderr, "CM\n");  
			    	exit(1); 
			   	}
			   
            	snps[m-1] = strdup(str);
	            
	            col = 1;
	            strLen = 0;
	        }
	        else if (col==0)
	        {
	            str[strLen++] = Xchar;
	        }
	        else if (Xchar == '\n')
	        {
	        	m++;
	        	col = 0;
	        	
	            continue;
	        }
	    }
    }
    NSAMPLES = n-1;
    pcNo = pcNo>NSAMPLES ? NSAMPLES : pcNo;
    nSNP = m-1;
    
    fprintf(stderr, "  No. of columns = %d\n", NSAMPLES);
    fprintf(stderr, "  No. of rows = %d\n", nSNP);

    fclose(fp);
    
    /* malloc */
    if((eval = (double *) malloc(NSAMPLES*sizeof(*eval))) == NULL)
    { fprintf(stderr,"CM\n");  exit(1); }
    if((evec = (double *) malloc(NSAMPLES*NSAMPLES*sizeof(*evec))) == NULL)
    { fprintf(stderr,"CM\n");  exit(1); }
    if((X = (double *) malloc(nSNP*NSAMPLES*sizeof(*X))) == NULL)
    { fprintf(stderr,"CM\n");  exit(1); }
    if((XTX = (double *) malloc(NSAMPLES*NSAMPLES*sizeof(*X))) == NULL)
    { fprintf(stderr,"CM\n");  exit(1); }

    /* initialize XTX */
    for(n=0; n<NSAMPLES; n++)
    {
        for(nn=0; nn<NSAMPLES; nn++)
            XTX[NSAMPLES*n+nn] = 0.0;
    }

    /* get data */
    if( (fp = fopen(INFILE, "r")) == NULL)
    {
        fprintf(stderr,"Could not open input file %s\n", INFILE);  exit(1);
    }

 	fprintf(stderr, "Reading matrix");
 
	m = -1;
    while(m<nSNP)
    {
		/*reads 1 line*/
		if(m==-1)
		{
			/*ignores header*/
			while (1)
			{
				fscanf(fp, "%c", &Xchar);
				if (Xchar=='\n')
				{
					break;
				}
			}			
		}
		else
		{
			/*ignore first column*/
			while (1)
			{
				fscanf(fp, "%c", &Xchar);
				if (Xchar=='\t')
				{
					break;
				}
			}
			
			/*reads in remaining row*/
			n=0;
			strLen=0;
			
			while (1)
			{
				fscanf(fp, "%c", &Xchar);
				if (Xchar=='\t')
				{
					str[strLen] = '\0';
					strLen = 0;
					X[m*NSAMPLES+n] = atof(str);
					if (X[m*NSAMPLES+n]==-1)
					{
						X[m*NSAMPLES+n] = -100;
					}
					
					n++;
					
					continue;
				}
				else if (Xchar=='\n')
				{
					str[strLen] = '\0';
					strLen = 0;
					X[m*NSAMPLES+n] = atof(str);
					if (X[m*NSAMPLES+n]==-1)
					{
						X[m*NSAMPLES+n] = -100;
					}
					
					n++;
					
					break;
				}
				else
				{
					str[strLen++] = Xchar;
				}
			}
		}

        m++;
    }
    
    fclose(fp);
    fprintf(stderr, " ... completed\n");
    
    printf("Normalization\n");
	printf("  1)Centering\n");
	if (individualNormalization)
	{
		printf("  2)Genetic drift rate\n");
	}
	else if (populationNormalization)
	{
		printf("  2)sqrt(pbar*(1-pbar))\n");
	}
		
    fprintf(stderr, "Constructing covariance matrix\n");
	if (rowCenter)
    {
    	for (m=0; m<nSNP ;m++)
    	{
	        /* mean-adjust this SNP */
	        rowvalid = 0;
	        rowsum = 0.0;
	        for(n=0; n<NSAMPLES; n++)
	        {
	        	if(X[m*NSAMPLES+n] >= -99.0)
        		{
	                rowvalid++;
	                rowsum += X[m*NSAMPLES+n];
	        	}
	        }
	    
	        if (individualNormalization)
            {	                	
            	rowsum /= 2;
            	rowmean = (rowsum)/((double)(rowvalid));
            	rowmeanbayes = (rowsum+0.5)/((double)(1+rowvalid));
        	}
	        else if (populationNormalization)
            {	                	
            	rowmean = (rowsum)/((double)(rowvalid));
            	rowmeanbayes = (rowsum+0.5)/((double)(1+rowvalid));
        	}
        	else
        	{
        		rowmean = (rowsum)/((double)(rowvalid));
        	}
	    
	        for(n=0; n<NSAMPLES; n++)
	        {	        	
	            if(X[m*NSAMPLES+n] >= -99.0)
	            {
	                if (individualNormalization)
	                {	                	
	                	X[m*NSAMPLES+n] = X[m*NSAMPLES+n]/2 - rowmean;
	                	X[m*NSAMPLES+n] /= sqrt(rowmeanbayes*(1.0-rowmeanbayes));
	            	}
	                else if (populationNormalization)
		            {	                	
	                	X[m*NSAMPLES+n] -= rowmean;
	                	X[m*NSAMPLES+n] /= sqrt(rowmeanbayes*(1.0-rowmeanbayes));
	                }
	            	else
	            	{
	            		X[m*NSAMPLES+n] -= rowmean;
	            	}
	            }
	            else
	            {
	                X[m*NSAMPLES+n] = 0.0;
	        	}
	        	
	        	/*printf("%f\n", X[m*NSAMPLES+n]);*/
	        }
	        
	        /* update XTX */
	        for(n=0; n<NSAMPLES; n++)
	        {
	            for(nn=n; nn<NSAMPLES; nn++)
	            {
	                XTX[NSAMPLES*n+nn] += X[m*NSAMPLES+n]*X[m*NSAMPLES+nn];
	        	}
	        }
    	}
	}
	/*column centre, NOT UPDATED*/
	else
	{	    
	    /*Mean adjust samples*/
	    for(n=0; n<NSAMPLES; n++)
	    {
	    	colvalid = 0;
	        colsum = 0.0;
	            
	        for(m=0; m<nSNP; m++)
	        {
	        	if(X[m*NSAMPLES+n] >= -99.0)
        		{
		        	colvalid++;
		            colsum += X[m*NSAMPLES+n];
	        	}
	        }
	        
	        colmean = (colsum)/((double)(colvalid));
	        colmeanbayes = (colsum+0.5)/((double)(1+colvalid));
	        
	        for(m=0; m<nSNP; m++)
	        {
	        	if(X[m*NSAMPLES+n] >= -99.0)
	            {
		        	X[m*NSAMPLES+n] -= colmean;
		        	if (individualNormalization)
		            {
		               	X[m*NSAMPLES+n] /= sqrt(colmeanbayes*(1.0-colmeanbayes));
		            }
	        	}
	        	else
	            {
	                X[m*NSAMPLES+n] = 0.0;
	        	}
	        }     
	    } 
		
	    /* update XTX */
	    for(m=0; m<nSNP; m++)
	    {
		    for(n=0; n<NSAMPLES; n++)
		    {
		        for(nn=n; nn<NSAMPLES; nn++)
		            XTX[NSAMPLES*n+nn] += X[m*NSAMPLES+n]*X[m*NSAMPLES+nn];
		    }
		}
	}    
    
    /* complete XTX */
    for(n=0; n<NSAMPLES; n++)
    {
        for(nn=n; nn<NSAMPLES; nn++)
        {
            XTX[NSAMPLES*n+nn] /= ((double)nSNP);
    	}
    }
    for(n=0; n<NSAMPLES; n++)
    {
        for(nn=0; nn<n; nn++)
        {
            XTX[NSAMPLES*n+nn] = XTX[NSAMPLES*nn+n];
    	}
    }
    
    /* singular value decomposition */
    fprintf(stderr, "Calculating eigen vectors and values\n");
    eigvecs(XTX, eval, evec, NSAMPLES); /* eigenvector k is evec[k*NSAMPLES+n] */       
    
    if (printCovarianceMatrix)
    {
	    /* print Covariance Matrix */
	    fprintf(stderr, "Printing covariance matrix\n");
	    for(n=0; n<NSAMPLES; n++) 
	    {
	    	for(nn=0; nn<NSAMPLES-1; nn++) 
		    {
		    	fprintf(fpcov,"%.06f\t", XTX[NSAMPLES*n+nn]);
		    }
		    
		    fprintf(fpcov,"%.06f\n", XTX[NSAMPLES*n+nn]);
	    }
	    fclose(fpcov);
	}

    fprintf(stderr, "Printing eigen vectors and values\n");
    /* print eval */
    sum = 0;
    for(k=0; k<NSAMPLES; k++) 
    {
    	sum += eval[k];
	}
	fprintf(fpeval,"PC\teigenvalue\tpercentage-of-variance\n",eval[k],eval[k]/sum);
    for(k=0; k<NSAMPLES; k++) 
    {
    	fprintf(fpeval,"PC%d\t%.06f\t%0.06f\n", k+1, eval[k],eval[k]/sum);
	}
	fclose(fpeval);

	/* print principal components */
	for(n=0; n<NSAMPLES; n++)
	{
		if(n==0)
		{
			if (tgFile)
			{
				fprintf(fpout, "sample-id");
			}
			else if (pafFile)
			{
				fprintf(fpout, "population-id");
			}
	    			
			for(k=0; k<pcNo; k++)
	    	{
   				fprintf(fpout, "\tPC%d", k+1);
		    }
		    
			fprintf(fpout, "\n");
		}
		
		fprintf(fpout, "%s", samples[n]);
		
	    for(k=0; k<pcNo; k++)
    	{
			fprintf(fpout, "\t%.04f", evec[k*NSAMPLES+n]);
	    }
	    
	    fprintf(fpout, "\n");
    }
    
    fclose(fpout);
    
    /* allocate memory to syyArray */
	if((syyArray = (double *) malloc(pcNo * sizeof(*syyArray))) == NULL)
    { fprintf(stderr,"CM\n");  exit(1); }
                        
    fprintf(stderr, "Printing SNP correlations\n");
	/* print SNP correlations */
	
	/* open snp-correlation file */
    if( (fpcor = fopen(CORFILE, "w")) == NULL)
    {
        fprintf(stderr,"Could not open cor file %s\n", CORFILE);  exit(1);
    }
        
	for(m=0; m<nSNP; m++)
	{   
		if(m==0)
		{
		    /* print header */
		    fprintf(fpcor, "snp-id");
		    
			for(k=0; k<pcNo; k++)
	    	{
				fprintf(fpcor, "\tPC%d", k+1);
		    }

			fprintf(fpcor, "\n");

		    for(k=0; k<pcNo; k++)
        	{   
        	    syy = 0;
        	    
        	    for(n=0; n<NSAMPLES; n++)
        	    {
        	        syy += evec[k*NSAMPLES+n] * evec[k*NSAMPLES+n];
        	    }
        	    
        	    syyArray[k] = syy;
    	    }
    	    
		}
        
		fprintf(fpcor, "%s", snps[m]);

        sxx = 0;
	    for(n=0; n<NSAMPLES; n++)
	    {
            sxx += X[m*NSAMPLES+n] * X[m*NSAMPLES+n];
	    }
		
		for(k=0; k<pcNo; k++)
    	{
    	    sxy = 0;
    	 	for(n=0; n<NSAMPLES; n++)
	        {   
    	        sxy += X[m*NSAMPLES+n] * evec[k*NSAMPLES+n];
            }
            
            if (sxx==0 || syyArray[k]==0)
            {
			    fprintf(fpcor, "\t0");
	        }
	        else
	        {
	            fprintf(fpcor, "\t%.04f", (sxy/sqrt(sxx*syyArray[k])));
	        }
	    }
	    
	    fprintf(fpcor, "\n");
    }   

    fclose(fpcor);
}
