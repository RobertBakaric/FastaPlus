/*
 * FastaPlusDemo.cpp
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
#include <string>
#include <Fasta/Fasta.hpp>
#include <boost/program_options.hpp>




namespace po = boost::program_options;
using namespace std;
using namespace fastaplus;

void PrintLogo(){
  
   cout <<"\n\
\n\
  ______        _        _____  _           \n\
 |  ____|      | |      |  __ \\| |          \n\
 | |__ __ _ ___| |_ __ _| |__) | |_   _ ___ \n\
 |  __/ _` / __| __/ _` |  ___/| | | | / __|\n\
 | | | (_| \\__ \\ || (_| | |    | | |_| \\__ \\\n\
 |_|  \\__,_|___/\\__\\__,_|_|    |_|\\__,_|___/\n\
                                           \n\
                            by Robert Bakaric \n\
\n\
_________________________________________________________________v0.01\n\
**********************************************************************\n\
\n\
CONTACT:\n\
   Program Fasta has been written and is maintained by Robert Bakaric,\n\
   email: rbakaric@irb.hr , bakaric@evolbio.mpg.de                    \n\
                                                                      \n\
\n\
\n\
LICENSE:\n\
   The program is distributed under the GNU General Public License.   \n\
   You should have received a copy of the licence together  with this \n\
   software. If not, see http://www.gnu.org/licenses/                 \n\
______________________________________________________________________\n\
**********************************************************************\n\
" << endl;

                                                           
}


template <typename INT, typename CHARA>
po::variables_map SetOptions(INT& argc, CHARA& argv){

    try {
        int opt;
        string version = "0.01";
        po::options_description desc("Options");
        desc.add_options()
            ("help,h", "produce help message")
            ("version,v", "print version information")
            ("input-file,i", po::value< string >(), "input file")
            ("taxid,t", po::value< string >(), "taxonomy identifier")
            ("output-file,o", po::value< string >(), "output file");

        po::positional_options_description p;
        p.add("input-file,i", -1);
        
        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).
                  options(desc).positional(p).run(), vm);
        po::notify(vm);
    
        if (vm.count("help")) {
            cout << "Usage: ./program [options]\n\n";
            PrintLogo();
            cout << desc;
            exit(0);
        }else if (vm.count("version")) {
            cout << "Program version:  " << version << "\n";
            exit(0);
        }
        if (!vm.count("input-file")){
            cout << "Input file is not defined \n";
            exit(0);
        }
        if (!vm.count("taxid")){
            cout << "Taxonomy identifier not defined \n";
            exit(0);
        }
        return vm;    
    }
    catch(std::exception& e)
    {
        cout << e.what() << "\n";
        exit(0);
    }    
}

template<typename T>
ostream& operator<<(ostream& os, const vector<T>& v)
{
    copy(v.begin(), v.end(), ostream_iterator<T>(os, " ")); 
    return os;
}



int main(int argc, char** argv){
   
   po::variables_map arg;
   arg = SetOptions(argc,argv);

   string in =  arg["input-file"].as<string>();
   string taxid =  arg["taxid"].as<string>();
   string out =  arg.count("output-file") ? arg["output-file"].as<string>() : "out"; ;
   Fasta<int> NewFastaObj(in, taxid);
   
   cout << "This file contains:" << NewFastaObj.GetObjSummary("TotSeq") << " sequences\n";
   cout << "This file contains:" << NewFastaObj.GetObjSummary("TotSeqSize") << " sequence characters\n";
  
   /* Load from inside a program always assumes a clean record */
   string id = NewFastaObj.LoadFastaRec("MySeq", "HI MY NAME IS ROBERT and this Is my Sequence!!", "10001");
                                                                
   unordered_map<string,string> ret = NewFastaObj.GetFastaOnly(id);
   
   cout << "I added this sequence to my container: " << ret[id] << endl;
   cout << "My name is: "                            << NewFastaObj.GetSubStr(id,15,20) << endl;
   cout << "And you can locate my sequence in "      << out  << "directory under " << id << " ID\n"; 
   
   NewFastaObj.DmpFastaAll(out);
   
   return 0;
}



