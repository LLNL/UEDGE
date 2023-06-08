
#define nispmx 31
	! maximum number of ion species
#define ngspmx 6
	! maximum number of gas species; also must set in com.v
#define nmcmx 12
        ! maximum number of EIRENE test species in data file 'fort.44'
#define ndomainmx 32
         ! maximum number of domains for domain decomposition
#define nxptmx 2
         ! maximum number of x-points in R-Z domain
#define ndcsmx nispmx*nxptmx
         ! data dimension for csfaclb and csfacrb
#define nispmxngspmx nispmx*ngspmx 
         ! tot numb ion*gas species
#define nstramx  10 
         ! maximum number of strata for MC neutrals code
