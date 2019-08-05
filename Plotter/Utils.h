#ifndef __UTILS_H__
#define __UTILS_H__
vector<TString> Split(TString s,TString del){
  TObjArray* array=s.Tokenize(del);
  vector<TString> out;
  for(const auto& obj:*array){
    out.push_back(((TObjString*)obj)->String());
  }
  array->Delete();
  return out;
}

TString Dirname(TString s){
  if(s.Last('/')!=-1) return s(0,s.Last('/'));
  else return ".";
}
TString Basename(TString s){
  if(s.Last('/')!=-1) return s(s.Last('/')+1,s.Length());
  else return s;
}
TString Replace(TString str,TRegexp reg,TString repl){
  int extent;
  int start=str.Index(reg,&extent);
  return str.Replace(start,extent,repl);
}

#endif
