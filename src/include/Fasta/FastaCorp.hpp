/*
 * FastaCorp.hpp
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
#include <unordered_set>

using namespace std;

/** @file FastaCorp.hpp
 * Classes and functions for handling Fasta formated sequence 
 */
 
/**
 * @brief The namespace of FastaPlus container.
 */

namespace fastaplus {

 /** @brief FastaCorp class processes the sequence of a given Fasta record
 */

class FastaCorp {
   
   unordered_map<string,string> Corpus;
   vector<string>              Identifiers;
   unordered_set<string>       CheckId;
   
/*!
 * ToUpperCase function converts lower case letters to their upper case form
 * @param Str [string&]
 * @par Example:
 * @code
 * string x = "AAAhagA";
 * ToUpperCase(x);
 * cout << x << endl; // result: AAAHAGAA
 * @endcode
 */
   inline void ToUpperCase(string& Str);
/*!
 * RemoveSpaces function removes spaces from strings
 * @param Str [string&]
 * @par Example:
 * @code
 * string x = "AAA   A";
 * RemoveSpaces(x);
 * cout << x << endl; // result: AAAAA
 * @endcode
 */
   inline void RemoveSpaces(string& Str);
/*!
 * MaskDubious function replaces all non alphabet characters with X characters
 * @param  Str [string&]
 * @par Example:
 * @code
 * string x = "AAA 8.A";
 * MaskDubious(x);
 * cout << x << endl; // result: AAAXXXAA
 * @endcode
 */
   inline void MaskDubious(string& Str);
   
public:

   FastaCorp();
/*!
 * FastaCorp class constructor for adding multiple fasta records. By 
 * default constructor assumes strings are not clean.
 * @param Corp [unordered_map<string,string>&] 
 * @par Example:
 * @code
 * FastaCorp fastaCorp(map);
 *   // where map is a construct of a string identifier and its content
 *   //  id => ATVYYWQEGGGESS...
 * @endcode
 */
   FastaCorp(unordered_map<string,string>& Corp);
/*!
 * FastaCorp class constructor overload for adding a single fasta record.
 * By default constructor assumes the string is not clean.
 * @param Id [const string&]
 * @param Corp [const string&]
 * @par Example:
 * @code  FastaCorp fastaCorp("ID", "ATVYYWQEGGGESS...");
 * @endcode
 */
   FastaCorp(const string& Id, const string& Corp);

/*!
 * FastaCorp class destructor.
 */
   ~FastaCorp();

/*!
 * LoadClean function loads strings into a container as they are
 * @param Corp [unordered_map<string,string>& ]
 */
   void LoadCleanCorp(unordered_map<string,string>& Corp);
/*!
 * LoadClean function overload loads a single string into a container as is
 * @param Id [const string& ]
 * @param Corp [const string& ]
 */
   void LoadCleanCorp(const string& Id, const string& Corp);
/*!
 * Load function cleans strings before loading them into a container
 * @param Corp [unordered_map<string,string>& ]
 */
   void LoadCorp(unordered_map<string,string>& Corp);
/*!
 * Load function overload cleans a string before loading it into a container
 * @param Id [const string& ]
 * @param Corp [const string& ]
 */
   void LoadCorp(const string& Id, const string& Corp);
/*!
 * GetCorpOnly function returns a specific string
 * @param Id [const string&]
 */
   unordered_map<string,string> GetCorpOnly(const string& Id);
/*!
 * GetCorpOnly function overload returns a set of specified strings
 * @param Ids [vector<string>&]
 */
   unordered_map<string,string> GetCorpOnly(vector<string>& Ids);
/*!
 * GetCorpAllExcept function returns all strings except the one scpecified
 * @param Id [const string&]
 */
   unordered_map<string,string> GetCorpAllExcept(const string& Id);
/*!
 * GetCorpAllExcept function overload returns all strings except the ones 
 * within a given vector
 * @param Ids [vector<string>&]
 */
   unordered_map<string,string> GetCorpAllExcept(vector<string>& Ids);
/*!
 * GetCorpAllExcept function returns all strings within a containor
 */
   unordered_map<string,string> GetCorpAll();
/*!
 * The function clears the containor.
 */
   void Clear();
};

FastaCorp::FastaCorp(){}

FastaCorp::FastaCorp(unordered_map<string,string>& Corp){
   LoadCorp(Corp);
}


FastaCorp::~FastaCorp(){
   Clear();
}


void FastaCorp::Clear(){
   Corpus.clear();
   vector<string>().swap(Identifiers);
}


FastaCorp::FastaCorp(const string& Id,const string & Corp){
   LoadCorp(Id,Corp);
}


void FastaCorp::LoadCorp(unordered_map<string,string>& Corp){
   for ( auto it = Corp.begin(); it != Corp.end(); ++it )
      LoadCorp(it->first, it->second);
}


void FastaCorp::LoadCorp(const string& Id,const string & Corp){
   string s = Corp;
   ToUpperCase(s);
   RemoveSpaces(s);
   MaskDubious(s);
   Corpus[Id] = s;
   if(CheckId.find(Id) != CheckId.end()){
      Identifiers.push_back(Id);
      CheckId.insert(Id);
   }
}



void FastaCorp::LoadCleanCorp(unordered_map<string,string>& Corp){
   for ( auto it = Corp.begin(); it != Corp.end(); ++it )
      LoadCleanCorp(it->first, it->second);

}

void FastaCorp::LoadCleanCorp(const string& Id,const string & Corp){
   Corpus[Id] = Corp;
   if(CheckId.find(Id) != CheckId.end()){
      Identifiers.push_back(Id);
      CheckId.insert(Id);
   }
}


unordered_map<string,string> FastaCorp::GetCorpAll(){
   return Corpus;
}


unordered_map<string,string> FastaCorp::GetCorpOnly(const string& Id){
   unordered_map<string,string> str;
   str[Id] = Corpus[Id];
   return str;
}


unordered_map<string,string> FastaCorp::GetCorpOnly(vector<string>& Ids){
   unordered_map<string,string> str;
   for(long i =0; i< Ids.size(); i++)
      str[Ids[i]] = Corpus[Ids[i]];
   return str;
}


unordered_map<string,string> FastaCorp::GetCorpAllExcept(const string& Id){
   unordered_map<string,string> str;
   for ( auto it = Corpus.begin(); it != Corpus.end(); ++it )
      if(Id.compare(it->first) != 0) 
         str[it->first] = it->second;
   return str;
}


unordered_map<string,string> FastaCorp::GetCorpAllExcept(vector<string>& Ids){

   vector<string> get;
   copy_if(Identifiers.begin(), Identifiers.end(), back_inserter(get),
     [&Ids](const string& arg) { return (find(Ids.begin(), Ids.end(), arg) == Ids.end());});
 
   return GetCorpOnly(get);
}


inline void FastaCorp::ToUpperCase(string& Str){
   transform(Str.begin(), Str.end(),Str.begin(), ::toupper);
}


inline void FastaCorp::RemoveSpaces(string& Str){
   Str.erase(remove_if(Str.begin(), Str.end(), ::isspace), Str.end());
}


inline void FastaCorp::MaskDubious(string& Str){
   replace_if(Str.begin(), Str.end(), [](char c) { return !isalpha(c); }, 'X' );
}
}
