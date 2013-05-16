library("Rcpp")

cppFunction(code = "
NumericVector gccount2(CharacterVector inseq) {
  Rcpp::CharacterVector dnaseq(inseq);
  Rcpp::NumericVector ans(4);
  std::string s = Rcpp::as<std::string>(dnaseq);

  for (int i = 0; i < s.size(); i++) {
    char p = s[i];
    if (p == \'A\') 
      ans[0]++;
    else if (p == \'C\') 
      ans[1]++;
    else if (p == \'G\') 
      ans[2]++;
    else if (p == \'T\') 
      ans[3]++;
    else 
      Rf_error(\"Wrong alphabet\");
  }  
  return(ans);
}
")
