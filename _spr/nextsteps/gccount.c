#include <R.h> 
#include <Rdefines.h>

SEXP gccount(SEXP inseq) {
  int i, l;
  SEXP ans, dnaseq;    
  PROTECT(dnaseq = STRING_ELT(inseq, 0)); 
  l = LENGTH(dnaseq); 
  printf("length %d\n",l);
  PROTECT(ans = NEW_NUMERIC(4));

  for (i = 0; i < 4; i++) 
    REAL(ans)[i] = 0;

  for (i = 0; i < l; i++) {
    char p = CHAR(dnaseq)[i];
    if (p=='A') 
      REAL(ans)[0]++;
    else if (p=='C') 
      REAL(ans)[1]++;
    else if (p=='G') 
      REAL(ans)[2]++;
    else if (p=='T') 
      REAL(ans)[3]++;
    else 
      error("Wrong alphabet");
  }
  UNPROTECT(2);
  return(ans);
}
