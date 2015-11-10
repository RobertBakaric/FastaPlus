/*
 * CovnertString.hpp
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
#include <sstream>

using namespace std;

namespace fastaplus{

/*!
 * StringToNumeric function converts a string to a specified numeric value type
 * @param  str [const string&]
 */

template<typename Tnum>
inline Tnum StringToNumeric(const string& str){
   
   Tnum      num;
   stringstream ss;
   
   ss << str; ss >> num;
   return num;
}
/*!
 * StringToBool function converts a numeric value to a bool
 * @param  str [const string&]
 */
inline bool StringToBool(const string& str){

   return str.compare("t") == 0 || str.compare("T") == 0 ? true : false;
}

}
