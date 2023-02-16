
#define NTEMP   48
#define NDEN    11
#define NCASET  40
#define NCASENO 40
#define NCASENT 40
#define NZSORMX 10
      ! maximum number of impurity sources
#define MXMISO 5
       ! maximum number of charged isotopes; also must set in com.v
#define MXNZCH 26
      ! maximum number of charge states for any isotope
#define MXMINZ MXMISO*MXNZCH
#define KXA 3
#define MXMISO1 MXMISO+1
#define KNX 3*KXA*KXA*MXNZCH
#define KMXZ KXA*MXMINZ
#define NBA 5
      ! used in fmombal

