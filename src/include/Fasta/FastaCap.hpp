/*
 * FastaCap.hpp
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



#include <iostream>
#include <algorithm>
#include <string>
#include <vector>
#include <unordered_map>

using namespace std;

/** @file FastaCap.hpp
 * Classes and functions for handling information within Fasta header 
 */
 
/**
 * @brief The namespace of FastaPlus container.
 */

namespace fastaplus {

 /** @brief FastaCap class handles the information located in the header line of a fasta record
 */
class FastaCap {
   
   unordered_map<string,string> SiToTi;
   unordered_map<string,string> SiToFeat;
   unordered_map<string,string> SiToSs;
   unordered_map<string,string> SsToSi;
   unordered_map<string,vector<string>> TiToSi;

/*!
 * ParseCap function extracts the information from a given string.
 * The extraction is based upon a given deliminator
 * @param Head [const string&]
 * @param del [char]
 * @par Example:
 * @code  ParseCap(">si|76612|ti|8363|ss|0|[tab]Additional information ",'|');
 * @endcode
 */
   vector <string> ParseCap(const string& Head, char del);
   
   
public:


   FastaCap();
/*!
 * FastaCap class constructor for adding multiple fasta records. By
 * default constructor assumes >si|***|ti|***|ss|***|[tab]Add... format. 
 * @param Caps [vector<string>&]
 * @par Example:
 * @code
 * vector<string> Caps ={">si|***|ti|***|ss|***|[tab]Add...",">si|***|ti|***|ss|***|[tab]Add...",">si|***|ti|***|ss|***|[tab]Add..."}
 * FastaCap fastaCap(Caps);
 * @endcode
 */
 
   FastaCap(vector<string>& Caps);

/*!
 * FastaCap class constructor for adding single fasta record. By
 * default constructor assumes >si|***|ti|***|ss|***|[tab]Add... format. 
 * @param Cap [const string&]
 * @par Example:
 * @code  FastaCap fastaCap(">si|***|ti|***|ss|***|[tab]Add...");
 * @endcode
 */
 
   FastaCap(const string& Cap);
   
/*!
 * FastaCap class destructor
 */
   ~FastaCap();

/*!
 * LoadCap function for adding single fasta record. By
 * default constructor assumes >si|***|ti|***|ss|***|[tab]Add... format. 
 * @param Cap [const string&]
 * @par Example:
 * @code  LoadCap(">si|***|ti|***|ss|***|[tab]Add...");
 * @endcode
 */
   void LoadCap(const string& Cap);
   
/*!
 * LoadCap function for adding multiple fasta records. By
 * default constructor assumes >si|***|ti|***|ss|***|[tab]Add... format. 
 * @param Caps [vector<string>&]
 * @par Example:
 * @code
 * vector<string> Caps ={">si|***|ti|***|ss|***|[tab]Add...",">si|***|ti|***|ss|***|[tab]Add...",">si|***|ti|***|ss|***|[tab]Add..."}
 * FastaCap fastaCap(Caps);
 * @endcode
 */
 
   void LoadCap(vector<string>& Caps);

/*!
 * GetCapSiForTi function returns all si identifiers for a given ti. 
 * @param Ti [const string&]
 * @par Example:
 * @code vector<string> vec = GetCapSiForTi("12345");
 * @endcode
 */ 
 
   vector<string>& GetCapSiForTi(const string& Ti);
   
/*!
 * GetCapSiForTi function overload returns all si identifiers for a given ti ones. 
 * @param Tis [vector<string>&]
 * @par Example:
 * @code vector<string> tiss = {"12345", "67890"}
 * vector<vector <string>> vec = GetCapSiForTi(tis);
 * @endcode
 */ 
 
   vector<vector <string>> GetCapSiForTi(vector<string>& Tis);

/*!
 * GetCapSiForSs function returns all si identifiers for a given ti. 
 * @param Ss [const string&]
 * @par Example:
 * @code vector<string> ss = GetCapSiForTi("12345");
 * @endcode
 */ 
 
   string& GetCapSiForSs(const string& Ss);
   
/*!
 * GetCapSiForSs function overload returns all si identifiers for a given ss ones. 
 * @param Sss [const string&]
 * @par Example:
 * @code vector<string> sss = {"12345", "67890"}
 * vector<vector <string>> vec = GetCapSiForTi(sss);
 * @endcode
 */ 
   vector <string> GetCapSiForSs(vector<string>& Sss);

/*!
 * GetCapSiForCap function extracts ss identifier from a given fasta indexed header. 
 * @param Cap [const string&]
 * @par Example:
 * @code  
 * string si= GetCapSiForCap(">si|12345|ti|***|ss|***|\tAdd...");
 * cout << si << endl;// prints: 12345
 * @endcode
 */   
   string GetCapSiForCap(const string& Cap);
/*!
 * GetCapSiForCap function extracts ss identifier from a given fasta indexed header. 
 * @param Caps [vector<string>&]
 * @par Example:
 * @code
 * vector<string> vec = {">si|12345|ti|***|ss|***|[tab]Add...",">si|67890|ti|***|ss|***|[tab]Add..."}
 * vector<string> si= GetCapSiForCap(vec);
 * // si contains: 12345, 67890
 * @endcode
 */ 
   vector<string> GetCapSiForCap(vector<string>& Caps);

/*!
 * GetCapTiForSi function returns ti identifier for a given si identifier. 
 * @param Si [string&]
 */ 
   string& GetCapTiForSi(const string& Si);
/*!
 * GetCapTiForSi function returns a set of ti identifiers for a given si identifier. 
 * @param Sis [vector<string>&]
 */ 
   vector <string> GetCapTiForSi(vector<string>& Sis);

/*!
 * GetCapSsForSi function returns si identifier for a given si identifier. 
 * @param Si [string&]
 */ 
   string& GetCapSsForSi(const string& Si);
/*!
 * GetCapSsForSi function returns a set of ss identifiers for a given si identifier. 
 * @param Sis [vector<string>&]
 */ 
   vector <string> GetCapSsForSi(vector<string>& Sis);

/*!
 * GetCapMetaForSi function returns associated meta information for a given si identifier. 
 * @param Si [string&]
 */ 
   string& GetCapMetaForSi(const string& Si);
/*!
 * GetCapMetaForSi function returns a set of associated meta information for a given set of si identifiers.
 * @param Sis [vector<string>&]
 */ 
   vector <string> GetCapMetaForSi(vector<string>& Sis);

/*!
 * GetCapAll function returns all header identifiers.
 */ 
   vector <string> GetCapAll();

/*!
 * The function eraces the containor.
 */
   void Clear();
   
};
   
FastaCap::FastaCap(){}

FastaCap::FastaCap(const string& Cap){
   LoadCap(Cap);
}

FastaCap::FastaCap(vector<string>& Caps){
   LoadCap(Caps);
}

FastaCap::~FastaCap(){
   Clear();
}

void FastaCap::LoadCap(const string& Cap){
   vector<string> First = ParseCap(Cap,'\t');
   vector<string> Second = ParseCap(First[0], '|');
   SiToTi[Second[1]]   = Second[3];
   SiToFeat[Second[1]] = First[1];
   SiToSs[Second[1]]   = Second[5];
   SsToSi[Second[5]]   = Second[1];
   TiToSi[Second[3]].push_back(Second[1]);
}

void FastaCap::LoadCap(vector<string>& Caps){
   for(long i = 0; i< Caps.size(); i++)
     LoadCap(Caps[i]);
}

vector<string> FastaCap::GetCapSiForCap(vector<string>& CapRaw){
   vector<string> res;
    for(long i = 0; i< CapRaw.size(); i++)
       res.push_back(GetCapSiForCap(CapRaw[i]));
   return res;
}
string FastaCap::GetCapSiForCap(const string& CapRaw){
   vector<string> First = ParseCap(CapRaw,'\t');
   vector<string> Second = ParseCap(First[0], '|');
   return Second[1];
}

vector <string> FastaCap::ParseCap(const string& Line, char del){
   stringstream ss(Line);
   string item;
   vector<string> tokens;
   while (getline(ss, item, del)) {
        tokens.push_back(item);
   }
   return tokens;
}

vector <string>& FastaCap::GetCapSiForTi(const string& Ti){
   return TiToSi[Ti];
}

vector<vector <string>> FastaCap::GetCapSiForTi(vector <string>& Tis){
   vector<vector <string>> res;
   for(long i=0; i< Tis.size(); i++)
      res.push_back(TiToSi[Tis[i]]);
   return res;
}

string& FastaCap::GetCapSiForSs(const string& Ss){
   return SsToSi[Ss];
}

vector <string> FastaCap::GetCapSiForSs(vector <string>& Sss){
   vector <string> res;
   for(long i=0; i< Sss.size(); i++)
      res.push_back(SsToSi[Sss[i]]);
   return res;
}

string& FastaCap::GetCapTiForSi(const string& Si){
   return SiToTi[Si];  
}

vector <string> FastaCap::GetCapTiForSi(vector<string>& Sis){
   vector <string> res;
   for(long i=0; i< Sis.size(); i++)
      res.push_back(GetCapTiForSi(Sis[i]));
   return res;
}
string& FastaCap::GetCapSsForSi(const string& Si){
   return SiToSs[Si];
}
vector <string> FastaCap::GetCapSsForSi(vector<string>& Sis){
   vector <string> res;
   for(long i=0; i< Sis.size(); i++)
      res.push_back(GetCapSsForSi(Sis[i]));
   return res;
}
  
string& FastaCap::GetCapMetaForSi(const string& Si){
   return SiToFeat[Si];
}

vector <string> FastaCap::GetCapMetaForSi(vector<string>& Sis){
   vector <string> res;
   for(long i=0; i< Sis.size(); i++)
      res.push_back(GetCapMetaForSi(Sis[i]));
   return res;
}

vector <string> FastaCap::GetCapAll(){
   vector <string> res;
   for ( auto it = SiToTi.begin(); it != SiToTi.end(); ++it )
      res.push_back(it->first);
   return res;
}

void FastaCap::Clear(){
   SiToTi.clear();
   SiToFeat.clear();
   SiToSs.clear();
   SsToSi.clear();
   TiToSi.clear();
}

}
