#include <Rcpp.h>

using namespace Rcpp;

RcppExport SEXP gccount2(SEXP inseq)
{
  Rcpp::IntegerVector ans(4);
  Rcpp::CharacterVector dnaseq(inseq);
  std::string s = Rcpp::as<std::string>(dnaseq[0]);

  for (int i = 0; i < s.size(); i++) {
    char p = s[i];
    if (p=='A') 
      ans[0]++;
    else if (p=='C') 
      ans[1]++;
    else if (p=='G') 
      ans[2]++;
    else if (p=='T') 
      ans[3]++;
    else 
      Rf_error("Wrong alphabet");
  }
  
  return(ans);
}
