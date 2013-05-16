#include <R.h> 
#include <Rdefines.h>

SEXP gccount(SEXP inseq) {
  int i, l; 
  char p;
  SEXP ans, dnaseq; 

  PROTECT(dnaseq = STRING_ELT(inseq, 0)); // a CHARSXP
  l = length(dnaseq);
  
  PROTECT(ans = allocVector(INTSXP, 4));
  memset(INTEGER(ans), 0, 4 * sizeof(int));

  for (i = 0; i < l; i++) {
    p = CHAR(dnaseq)[i];
    if (p == 'A')
      INTEGER(ans)[0]++;
    else if (p == 'C')
      INTEGER(ans)[1]++;
    else if (p == 'G')
      INTEGER(ans)[2]++;
    else if (p == 'T')
      INTEGER(ans)[3]++;
    else
      error("Wrong alphabet");
  }
  UNPROTECT(2);
  return(ans);
}

