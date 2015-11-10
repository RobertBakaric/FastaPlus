/*
 * XNU.hpp
 * 
 * Copyright 2015 Robert Bakaric <rbakaric@irb.hr>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 * 
 */

/* NOTE:
  * 
  * This library has been written according to xnu(.c) program:
  * 
  *    Jean Michel Claverie & David J. States (1993)
  *    Computers and Chemistry 17: 191-201.
  * 
  * In order to preserve the logic, underlying  data structures and code, snippetss  
  * used here are based on the exact code fragments form xnu.c (and associated) 
  * source files. Therefore  labeles are marked using same or similar variable names. 
  * 
  * */



#include <cctype>
#include <cmath>


#include <Filters/XNUData.hpp>

#include <vector>
#include <cstring>
#include <string>

namespace fastaplus {
   
/**
 * @brief XNU filter class.
 */
template <typename Tint>
class XNU : protected XnuScores{
   
   Tint ascend;
   Tint descend;
   double K;
   Tint H;

   Tint ncut;
   Tint mcut;

   double pcut;
   Tint scut;

   char subchar;
   Tint  repeats;

   
   double Lambda;
   vector<vector<Tint>> Mtx;
/*!
 * Function returns the index value associated to the AA alphabet
 * @param ch [char] // AA character
 * @param Alph [const string&] // "ARNDCQEGHILKMFPSTWYVBZX*-"
 */
   Tint AlphaToNum(char ch,  string Alph);
/*!
 * Function sets the expected parameters given a chosen matrix
 * @param f [vector<double>&] 
 */
   double EInfo(vector<double>& f);
   
   public:
   
/*!  
 * Default constructor
 */
 
   XNU();
         
/*! 
 * Constructor overload
 * @param Param [unordered_map<string,string>&] 
 */
 
   XNU(unordered_map<string,string>& Param);
   
/*! 
 * Destructor 
 */
   ~XNU();
/*! Function executing filtering procedure 
 * @param str [const string&] 
 */
  string Filter(const string& str); 
   
};

template <typename Tint>
XNU<Tint>::XNU():subchar('X'),scut(0), pcut(0.01), repeats(0), mcut(1), ncut(4), descend(1), ascend(1),K(0.2), XnuScores(){
   Lambda = this->Lambda60;
   Mtx    = this->Pam60;
   H = EInfo(this->Dayhoff);   
};

template <typename Tint>
XNU<Tint>::XNU(unordered_map<string,string>& Arg):XnuScores(){
   subchar   = (Arg.find("subchar") !=Arg.end()) ? Arg["subchar"][0]                     : 'X';
   scut      = (Arg.find("scut") !=Arg.end())    ? StringToNumeric<Tint>(Arg["scut"])    : 0;
   pcut      = (Arg.find("pcut") !=Arg.end())    ? StringToNumeric<double>(Arg["pcut"]) : 0.01;
   repeats   = (Arg.find("repeats") !=Arg.end()) ? StringToNumeric<Tint>(Arg["repeats"]) : 0;
   mcut      = (Arg.find("mcut") !=Arg.end())    ? StringToNumeric<Tint>(Arg["mcut"])    : 1;
   ncut      = (Arg.find("ncut") !=Arg.end())    ? StringToNumeric<Tint>(Arg["ncut"])    : 4;
   descend   = (Arg.find("descend") !=Arg.end()) ? StringToNumeric<Tint>(Arg["descend"]) : 1;
   ascend    = (Arg.find("ascend") !=Arg.end())  ? StringToNumeric<Tint>(Arg["ascend"])  : 1;
   
   if (Arg.find("pam") !=Arg.end()){
      if(Arg["pam"].compare("PAM60") == 0){
         Lambda = this->Lambda60;
         Mtx    = this->Pam60;
      }else if(Arg["pam"].compare("PAM12") == 0){
         Lambda = this->Lambda120;
         Mtx    = this->Pam120;
      }else if(Arg["pam"].compare("PAM250") == 0){
         Lambda = this->Lambda250;
         Mtx    = this->Pam250;
      }
   }else{
      Lambda = this->Lambda60;
      Mtx    = this->Pam60;
   }
   
   
   K = 0.2;
   H = EInfo(this->Dayhoff);
   

}

template <typename Tint>
XNU<Tint>::~XNU(){};

template <typename Tint>
double XNU<Tint>::EInfo(vector<double>& f){

	double sum,tot,fij,eij;

	sum = tot = 0.0;

	for (Tint i=0; i<Mtx.size()-5; i++)
	   for (Tint j=0; j<Mtx.size()-5; j++) {

   		fij = f[i]*f[j];
	   	tot += fij;

		   eij = Mtx[i][j]*fij*exp(Lambda*Mtx[i][j]);
		   sum += eij;
	   }

	return(Lambda*sum/tot);
}

template <typename Tint>
Tint  XNU<Tint>::AlphaToNum(char ch,  string Alph){
   /* no lower characters !!!! */
        char *p;
        if (isupper(Alph[0]) && islower(ch)) ch = toupper(ch);
        if (islower(Alph[0]) && isupper(ch)) ch = tolower(ch);
        p = strchr(&Alph[0],ch);
        return (p != NULL) ? (p - &Alph[0]) : (22);
}



template <typename Tint>
string XNU<Tint>::Filter(const string& s){
   
	Tint off = 0,sum = 0,beg = 0,end= 0,top= 0,noff=0;
	Tint topcut=0,fallcut=0;
   string str = s;
	double s0;
   
   vector<Tint> iseq(str.size()+1,0);
   vector<Tint> hit(str.size()+1,0);
   
	for (Tint i=0; i<str.size(); i++){
      iseq[i] = AlphaToNum(str[i], this->Alphabet);
   }
	noff = str.size()-1;
	if (ncut>0) noff=ncut;
   
	if (scut!=0) {
		topcut = scut;
	} else {
		s0 = 0 - log( pcut*H / (noff*K) ) / Lambda;
		if (s0>0) topcut = floor(s0 + log(s0)/Lambda + 0.5);
		else topcut = 0;
	}
	fallcut = static_cast<Tint>(log(K/0.001)/Lambda);

	for (Tint off=mcut; off<=noff; off++) {

      sum=top=0;
      beg=off;
      end=0;

		for (Tint i=off; i<str.size(); i++) {
			sum += Mtx[iseq[i]][iseq[i-off]];
			if (sum>top) {
				top=sum;
				end=i;
			}
			if (top>=topcut && top-sum>fallcut) {
				for (Tint k=beg; k<=end; k++) {
					if (ascend) hit[k] = 1;
					if (descend) hit[k-off] = 1;
				}
				sum=top=0;
				beg=end=i+1;
			} else if (top-sum>fallcut) {
				sum=top=0;
				beg=end=i+1;
			}
			if (sum<0) {
				beg=end=i+1;
				sum=top=0;
			}
		}
		if (top>=topcut) {
			for (Tint k=beg; k<=end; k++) {
				if (ascend) hit[k] = 1;
				if (descend) hit[k-off] = 1;
			}
		}
	}

	for (Tint i=0; i<str.size(); i++) {
		str[i] = toupper(str[i]);
		if (hit[i] ^ repeats) {
			if (subchar==0) {
				str[i] = tolower(str[i]);
			} else {
				str[i] = subchar;
			}
		}
	}
   
   return str;
   
}

}
