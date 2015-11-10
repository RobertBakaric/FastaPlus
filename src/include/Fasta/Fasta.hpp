/*
 * Fasta.hpp
 * 
 * Copyright 2016 Robert Bakaric <rbakaric@irb.hr>
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
#include <fstream>
#include <sstream>
#include <dirent.h>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <unistd.h>
#include <stdexcept>
#include <string>
#include <vector>
#include <unordered_map>
#include <Fasta/FastaCap.hpp>
#include <Fasta/FastaCorp.hpp>



/** @file Fasta.hpp
 * Utility classes and functions for Fasta files/records handling
 */



using namespace std;


/**
 * @brief The namespace of FastaPlus container.
 */

namespace fastaplus {


/**
 * @brief Class for handling fasta records.
 */
template <typename Tint>
class Fasta :  public FastaCorp,  public  FastaCap{
   
   Tint sid;  // remember ss is 345 -> start position  and 0-> whole seq
   Tint NumOfSeq;
   Tint TotSize;
/*!
 * CapToIndex function creates an index structure given a raw headder fasta line.
 * @param Cap [const string&]
 * @param Taxid [const string&]
 * @par Example:
 * @code  string s = CapToIndex(">ENS937474 additional information", "9606"); 
 * @return cout << s << endl; // si|00300432122/0|ti|300|ss|0|[tab]ENS937474 additional information
 * @endcode
 */
   string CapToIndex(const string& Cap,const string& TaxId,const string& Ss);
/*!
 * LoadFasta function loads (multi)fasta records from a file.
 * @param File [const string&]
 * @param Taxid [const string&]
 */
   void LoadFasta(const string& File, const string& TaxId);
/*!
 * Dmp function saves (multi)fasta records to a given file.
 * @param File [const string&]
 * @param Seqs [unordered_map<string, string>&]
 */
   void Dmp(const string& File, unordered_map<string, string> Seqs);

public:

/*!
 * Fasta class constructor 
 */
   Fasta();
/*!
 * Fasta class constructor overload .\n
 * Constructor assumes raw header line in each record
 * @param File [const string&]
 * @param TaxId [const string&]
 */
   Fasta(const string& File, const string& TaxId);
/*!
 * Fasta class constructor overload.\n
 * Constructor assumes formated header line in each record.
 * [Ex: >si|***|ti|***|ss|***|[tab]Add... ]
 * @param File [const string&]
 */
   Fasta(const string& File);
/*!
 * Fasta class desctructor 
 */
   ~Fasta();
/*!
 * Fasta file loader. \n
 * The loader assumes formatted fasta header: [Ex: >si|***|ti|***|ss|***|[tab]Add...]
 *  @param File [const string&]
 */
   void LoadFastaFile(const string& File);
/*!
 * Fasta file loader. \n
 * The loader assumes raw fasta header.
 * @param File [const string&]
 * @param TaxId [const string&]
 */
   void LoadFastaFile(const string& File, const string& TaxId);
/*!
 * Fasta record loader. \n
 * The loader assumes formatted fasta header: [Ex: >si|***|ti|***|ss|***|[tab]Add...]
 * @param Cap [const string&]
 * @param Corp [const string&]
 */
   void LoadFastaRec(const string& Cap, const string& Corp);
/*!
 * Fasta record loader. \n
 * The loader assumes raw fasta header with only 3 parameters specified.
 * @param Cap [const string&]
 * @param Corp [const string&]
 * @param TaxId [const string&]
 */
   string LoadFastaRec(const string& Cap, const string& Corp,const string& TaxId);
/*!
 * Fasta record loader. \n
 * The loader assumes raw fasta header with all 4 parameters specified.
 * @param Cap [const string&]
 * @param Corp [const string&]
 * @param TaxId [const string&]
 * @param Ss [const string&]
 */
   string LoadFastaRec(const string& Cap, const string& Corp, const string& TaxId, const string& Ss);
/*!
 * Fasta record loader. \n
 * The loader allows loading from a map with the assumption that fasta headers are formatted: [Ex: >si|***|ti|***|ss|***|[tab]Add...]
 * @param Records [unordered_map<string, string>&]
 */
   void LoadFastaRec(unordered_map<string, string>& Records);
/*!
 * Fasta record loader. \n
 * The loader allows loading from a map with the assumption that fasta headers are not indexed.
 * @param Records [unordered_map<string, string>&]
 * @param TaxId [const string&]
 */
   vector<string> LoadFastaRec(unordered_map<string, string>& Records, const string& TaxId);

/*!
 * GetCorp retrieves all sequences assigned to a given taxonomy identifier. \n
 * @param TaxId [const string&]
 */
   unordered_map<string,string> GetCorp(const string& TaxId);
/*!
 * Fasta record dumper. \n
 * The dumper retrieves all fasta records from the database assigned to a specified taxonomy identifier and saves them to a given file location.
 * @param File [const string&]
 * @param TaxId [const string&]
 */
   void DmpFastaAll(const string& File,const string& TaxId);
/*!
 * Fasta record dumper. \n
 * The dumper retrieves all fasta records from the database and saves them to a given file location.
 * @param File [const string&]
 */
   void DmpFastaAll(const string& File);
/*!
 * Fasta record dumper. \n
 * The dumper retrieves all fasta records from the database except the one in Cap and saves them to a given file location.
 * @param File [const string&]
 * @param Cap [const string&]
 */
   void DmpFastaAllExcept(const string& File, const string& Cap);
/*!
 * Fasta record dumper. \n
 * The dumper retrieves all fasta records from the database except those in Caps and saves them to a given file location.
 * @param File [const string&]
 * @param Caps [const vector<string>&]
 */
   void DmpFastaAllExcept(const string& File, const vector<string>& Caps);
/*!
 * Fasta record dumper. \n
 * The dumper retrieves only the fasta record from the database specified by Cap and saves them to a given file location.
 * @param File [const string&]
 * @param Cap [const string&]
 */
   void DmpFastaOnly(const string& File, const string& Cap);
/*!
 * Fasta record dumper. \n
 * The dumper retrieves only those fasta record from the database specified in Caps and saves them to a given file location.
 * @param File [const string&]
 * @param Caps [const vector<string>&]
 */
   void DmpFastaOnly(const string& File, const vector<string>& Caps);
   

   /*!
 * Fasta record getter. \n
 * The getter retrieves all fasta records from the database assigned to a specified taxonomy identifier and saves them to a given file location.
 * @param TaxId [const string&]
 */
   unordered_map<string,string> GetFastaAll(const string& TaxId);
/*!
 * Fasta record getter. \n
 * The getter retrieves all fasta records from the database and saves them to a given file location.
 */
   unordered_map<string,string> GetFastaAll();
/*!
 * Fasta record getter. \n
 * The getter retrieves all fasta records from the database except the one in Cap and saves them to a given file location.
 * @param Cap [const string&]
 */
   unordered_map<string,string> GetFastaAllExcept(const string& Cap);
/*!
 * Fasta record getter. \n
 * The getter retrieves all fasta records from the database except those in Caps and saves them to a given file location.
 * @param Caps [const vector<string>&]
 */
   unordered_map<string,string> GetFastaAllExcept(const vector<string>& Caps);
/*!
 * Fasta record getter. \n
 * The getter retrieves only the fasta record from the database specified by Cap and saves them to a given file location.
 * @param Cap [const string&]
 */
   unordered_map<string,string> GetFastaOnly(const string& Cap);
/*!
 * Fasta record getter. \n
 * The getter retrieves only those fasta record from the database specified in Caps and saves them to a given file location.
 * @param Caps [const vector<string>&]
 */
   unordered_map<string,string> GetFastaOnly(const vector<string>& Caps);
   
   
/*!
 * The function clears the containor.
 */
   void Clear();
/*!
 * Object data getter. \n
 * Getter retrieves summary information of the object.
 * @param What [const string&]
 */
   Tint GetObjSummary(const string& What );
   
/*!
 * Fasta record getter. \n
 * Getter retrieves a particular substring segment from a given fasta string.
 * @param Cap [const string&]
 * @param Start [const Tint]
 * @param Stop [const Tint]
 */
   string GetSubStr(const string& Cap, const Tint Start, const Tint Stop );

};



template <typename Tint>
Fasta<Tint>::Fasta():sid(0),NumOfSeq(0),TotSize(0){}


template <typename Tint>
Fasta<Tint>::Fasta(const string& File, const string& TaxId):sid(0),NumOfSeq(0),TotSize(0){

   LoadFasta(File,TaxId);
}

template <typename Tint>
Fasta<Tint>::Fasta(const string& File):sid(0),NumOfSeq(0),TotSize(0){
   sid=0;
   LoadFasta(File,"INDEXED");
}

template <typename Tint>
void Fasta<Tint>::LoadFastaRec(const string& Cap, const string& Corp){
   this->LoadCap(Cap);
   this->LoadCorp(Cap, Corp);
}

template <typename Tint>
string Fasta<Tint>::LoadFastaRec(const string& Cap, const string& Corp,const string& TaxId){
   return LoadFastaRec(Cap,Corp, TaxId, "0");
}
template <typename Tint>
string Fasta<Tint>::LoadFastaRec(const string& Cap, const string& Corp,const string& TaxId, const string& Ss){
   string head = CapToIndex(Cap, TaxId, Ss);
   this->LoadCap(head);
   string SI(this->GetCapSiForCap(head));
   this->LoadCleanCorp(SI, Corp);
   return SI;
}

template <typename Tint>
void Fasta<Tint>::LoadFastaRec(unordered_map<string, string>& Records){
   for ( auto it = Records.begin(); it != Records.end(); ++it )
      LoadFastaRec(it->first,it->second);
   
}

template <typename Tint>
vector<string> Fasta<Tint>::LoadFastaRec(unordered_map<string, string>& Records, const string& TaxId){
   vector<string> ret;
   for ( auto it = Records.begin(); it != Records.end(); ++it )
      ret.push_back(LoadFastaRec(it->first,it->second,TaxId, "0"));
   return ret;
}

template <typename Tint>
void Fasta<Tint>::LoadFastaFile(const string& File, const string& TaxId){
   LoadFasta(File,TaxId);
}

template <typename Tint>
void Fasta<Tint>::LoadFastaFile(const string& File){
   LoadFasta(File,"INDEXED");
}

template <typename Tint>
unordered_map<string,string> Fasta<Tint>::GetCorp(const string& TaxId){
    vector <string> myIds =  this->GetCapSiForTi(TaxId);
    return this->GetCorpOnly(myIds);
}

template <typename Tint>
void Fasta<Tint>::DmpFastaAll(const string& File){
   Dmp(File,GetFastaAll());
}

template <typename Tint>
void Fasta<Tint>::DmpFastaAll(const string& File, const string& TaxId){
   unordered_map<string,string> ret = GetCorp(TaxId);
   Dmp(File,ret);
}

template <typename Tint>
void Fasta<Tint>::DmpFastaAllExcept(const string& File, const string& Cap){
   Dmp(File,GetFastaAllExcept(Cap));
}

template <typename Tint>
void Fasta<Tint>::DmpFastaAllExcept(const string& File, const vector<string>& Caps){
   Dmp(File,GetFastaAllExcept(Caps));
}

template <typename Tint>
void Fasta<Tint>::DmpFastaOnly(const string& File, const string& Cap){
   Dmp(File,GetFastaOnly(Cap));
}

template <typename Tint>
void Fasta<Tint>::DmpFastaOnly(const string& File, const vector<string>& Caps){
   Dmp(File,GetFastaOnly(Caps));
}



template <typename Tint>
unordered_map<string,string> Fasta<Tint>::GetFastaAll(){
   return this->GetCorpAll();
}

template <typename Tint>
unordered_map<string,string> Fasta<Tint>::GetFastaAll(const string& TaxId){
   return GetCorp(TaxId);
}

template <typename Tint>
unordered_map<string,string> Fasta<Tint>::GetFastaAllExcept( const string& Id){
   return this->GetCorpAllExcept(Id);
}

template <typename Tint>
unordered_map<string,string> Fasta<Tint>::GetFastaAllExcept(const vector<string>& Ids){
   return this->GetCorpAllExcept(Ids);
}

template <typename Tint>
unordered_map<string,string> Fasta<Tint>::GetFastaOnly( const string& Id){
   return this->GetCorpOnly(Id);
}

template <typename Tint>
unordered_map<string,string> Fasta<Tint>::GetFastaOnly(const vector<string>& Ids){
   return this->GetCorpOnly(Ids);
}




template <typename Tint>
void Fasta<Tint>::Dmp(const string& File, unordered_map<string,string> Seqs){
   sid = 0;
   fstream fs;
   string line, Cap, Corp;
   
   fs.open (File.c_str(), ios::out);
   if ( !fs.is_open())
      throw runtime_error ("Cannot open file: " + File );
   for ( auto it = Seqs.begin(); it != Seqs.end(); ++it ){
      fs << ">si|" << it->first << "|ti|" << this->GetCapTiForSi(it->first) << "|ss|" 
         << this->GetCapSsForSi(it->first) << "|\t"<<this->GetCapMetaForSi(it->first) <<endl;
      for (Tint i = 0; i < (it->second).size(); i += 80) {
         fs << (it->second).substr(i, 80) << endl;
      }
   }
   fs.close();
}

template <typename Tint>
void Fasta<Tint>::LoadFasta(const string& File, const string& TaxId){
   sid = 0;
   fstream fs;
   string line, Cap, Corp;
   
   fs.open (File.c_str(), ios::in);
   if ( !fs.is_open())
      throw runtime_error ("Cannot open file: " + File );
      
   while ( !getline (fs,line).eof() ){
      if (line[0] == '>'){
         if(Cap.size() > 0){
            string head;
            if(TaxId.compare("INDEXED") == 0)
               head = Cap;
            else
               head = CapToIndex(Cap, TaxId, "0");
            this->LoadCap(head);
            this->LoadCorp(this->GetCapSiForCap(head), Corp);
            NumOfSeq ++;
            TotSize += Corp.size();
            Corp = "";
         }
         Cap = line.substr(1);
      }else{
         Corp += line;
         
      }
   }
   if(Cap.size() > 0){
      string head;
      if(TaxId.compare("INDEXED") == 0)
         head = Cap;
      else
         head = CapToIndex(Cap, TaxId, "0");
      this->LoadCap(head);
      this->LoadCorp(this->GetCapSiForCap(head), Corp);
      NumOfSeq ++;
      TotSize += Corp.size();
   }
   
   fs.close();
}

template <typename Tint>
inline string Fasta<Tint>::CapToIndex(const string& Cap, const string& TaxId, const string& Si){

   sid++;
   stringstream ss;
   stringstream ssid;
   ssid << sid;
   string comp = TaxId + ssid.str()+"/"+Si;
   ss << setw(30) << setfill('0') << comp;
   
   return "si|"+ss.str()+"|ti|"+TaxId+"|ss|"+Si+"|\t"+Cap;
}


template <typename Tint>
Tint Fasta<Tint>::GetObjSummary(const string& What ){
   if(What.compare("TotSeq") == 0)
      return NumOfSeq;
   else if (What.compare("TotSeqSize") == 0)
      return TotSize;
}

template <typename Tint>
string Fasta<Tint>::GetSubStr(const string& Cap, const Tint Start, const Tint Stop ){
   unordered_map<string,string> ret = this->GetCorpOnly(Cap);
   return ret[Cap].substr((Start-1),(Stop-Start+1));
}

template <typename Tint>
void Fasta<Tint>::Clear(){
   sid = 0;
   NumOfSeq = 0;
   TotSize = 0;
   FastaCap::Clear();
   FastaCorp::Clear();
}

template <typename Tint>
Fasta<Tint>::~Fasta(){
   Clear();
}

}
