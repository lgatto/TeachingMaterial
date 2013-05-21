library("Rcpp")

cppFunction("
IntegerVector ingccount2(CharacterVector inseq) {
  IntegerVector ans(4);
  std::string s = Rcpp::as<std::string>(inseq[0]);
  int n = inseq[0].size();
  for (int i = 0; i < n; i++) {
    if (s[i] == 'A') 
      ans[0]++;
    else if (s[i] == 'C') 
      ans[1]++;
    else if (s[i] == 'G') 
      ans[2]++;
    else if (s[i] == 'T') 
      ans[3]++;
    else 
      Rf_error(\"Wrong alphabet\");
  }
  return wrap(ans);
}
")

