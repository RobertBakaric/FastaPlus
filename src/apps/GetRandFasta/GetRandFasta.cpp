/*
 * GetRandFasta.cpp
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



 #include <string>
 #include <vector>
 #include <iostream>
 #include <fstream>
 #include <cstdlib>
 #include <unordered_map>
 #include <Utility/ConvertString.hpp>
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
   This program has been written and is maintained by Robert Bakaric,\n\
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
            ("taxid,t", po::value< string >(), "taxid")
            ("output-file,o", po::value< string >(), "output file")
            ("number,l", po::value< string >(), "The number of random sequences to be retrieved.")
        ;

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
        if (!vm.count("number")){
            cout << "The number of random sequences to be retrieved not specified \n";
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


 
 

 
 int main(int argc, char **argv){
   
   po::variables_map arg;
   arg = SetOptions(argc,argv);
  
   
try{
   string in      =  arg["input-file"].as<string>();
   string taxid   =  arg.count("taxid") ? arg["taxid"].as<string>() : "1";
   string num     =  arg["number"].as<string>();
   string output   = arg.count("output-file") ? arg["output-file"].as<string>() : "";  

   
   Fasta<int> NewFastaObj(in, taxid);
   unordered_map<string, string> fasta = NewFastaObj.GetFastaAll();
   unordered_map<string,string> keys; 
   
   if( StringToNumeric<int>(num) > fasta.size())
      throw runtime_error ("The number of fasta sequences is smaller than the number you have choosen:" +  num);
      
   ofstream fs;
   streambuf *backup;
   
   if(output.size()>0){
      backup = cout.rdbuf();
      fs.open (output.c_str());
      if ( !fs.is_open())
         throw runtime_error ("Cannot open file: " + output );
      cout.rdbuf(fs.rdbuf());
       
   }

   for(int x = 0; x< StringToNumeric<int>(num); ){
      srand(time(NULL)+x);
      auto random_it = next(begin(fasta), rand()%(fasta.size()-0 + 1) + 0);
      if(keys.find(random_it->first) == keys.end()){
         cout <<  NewFastaObj.GetCapMetaForSi(random_it->first) << endl;
         cout << random_it->second << endl;
         keys[random_it->first] = "check";
         x++;
      } 
   }
   
   if ( fs.is_open()){
      cout.rdbuf(backup);
      fs.close();
   }
}catch(runtime_error& e){
      cerr << e.what() << "\n";
}
   
   
    return 0;
 }
