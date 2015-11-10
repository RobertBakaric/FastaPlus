/*
 * SGE.hpp
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
 
 
 /* NOTE:
  * 
  * This library has been written according to seg(.c) program:
  * 
  *    Wootton, J.C., Federhen, S. (1993)  Statistics of local complexity
  *    in amino acid sequences and sequence  databases.  Computers &
  *    Chemistry 17: 149-163.
  * 
  * In order to preserve the logic, underlying  data structures and code, snippetss  
  * used here are based on the exact code fragments form seg.c (and associated) 
  * source files. Therefore  labeles are marked using same or similar variable names. 
  * 
  * */
 
 
#include <vector>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <math.h>
#include <unordered_map>
#include <Filters/LnFact.hpp>
#include <Utility/ConvertString.hpp>



using namespace std; 

/**
 * @brief The namespace of FastaPlus container.
 */

namespace fastaplus {

/*FIXME: Remove legacy code !!*/

/*FIXME: This needs to be fragmented. 
 * Modularity above all !!*/
/**
 * @brief SEG AA sequence filter.
 */
template <typename Tint>
class SEG {

/* Globals */
  Tint    SegWindow;
  double SegLocut;
  double SegHicut;
  Tint    MaxX;
  Tint    MaxTrim; 
  Tint    Period;
  Tint    MergeOverlaps;


/* Structures */
  struct Alphabet{  
    Tint   alphasize;           /* size */
    double lnalphasize;         /* ln(size) */
    Tint*  alphaindex;          /* ASCII code table*/ 
    unsigned char* alphaflag;   /* array of indicators if the characer is (not) AA */
  } ;

  struct CSeq{  
    struct CSeq* parent;        /* current one */
    char*  seq;                 /* AA sequence */
    Alphabet* palpha;           /* alphabet info */
    Tint   start;               /* starting for seg. */
    Tint   length;              /* sequence length */
    Tint   Xes;                 /* the number of X's */  
    Tint*  charfreq;            /* number of characters in a string */
    Tint*  state;                
    double entropy;
  } ;

  struct SeqSeg{ 
    Tint begin;          
    Tint end;            
    struct SeqSeg *next;  
  } ;
  

  Alphabet* alpha;
  
  

/* Functions */
   
  void     SegFree(SeqSeg* seg);
  void     MergeSegs(CSeq* seq, SeqSeg* segs);
  template<typename Targ>
  void     SetParamaters(Targ &arg);
  CSeq*    NewCSeq();
  double*  ComputeEntropy(CSeq* seq,Tint  first, Tint last, Tint downset, Tint upset);
  Tint     SegSeq(CSeq* seq, SeqSeg **segs, Tint offset);
  Tint     LocLow(Tint i, Tint limit, double* H);
  Tint     LocHigh(Tint i, Tint limit, double* H);
  CSeq*    OpenWin(CSeq* parent, Tint start, Tint length);
  void     StateOn(CSeq* win);
  void     CompOn(CSeq* win);
  Tint     Trim(CSeq* seq, Tint* leftend, Tint* rightend);
  void     CloseWin(CSeq* win);
  bool     ShiftWin1(CSeq* win);
  double   Entropy(Tint* sv);
  void     DecrementSV(Tint* sv, Tint clas);
  void     IncrementSV(Tint* sv, Tint clas);
  double   GetProb(const Tint* sv, Tint total, const Alphabet* palpha);
  double   LnPerm(const Tint* sv, Tint window_length);
  double   LnAss(const Tint* sv, Tint alphasize);
  double   lnFact(Tint n) ;
  void     CSeqFree(CSeq* seq);
  Alphabet* MakeAlpha();
  void     AlphaFree (Alphabet* alpha);
  template <typename T>
  void     SafeFree(T *x);
  template <typename T>
  void     SafeFree(T **x);
  void     EntropyOn(CSeq* win);
  
   
  public:

/*!
 * SEG class constructor. 
 * @param arg [unordered_map<string,string>&]
 */ 
  template <typename Targ>
    SEG(Targ& arg);
/*!
 * SEG class default constructor. 
 */
    SEG();

/*!
 * SEG class destructor. 
 */
  ~SEG();

/*!
 * Filter function identifies and masks (xXx) low complexity segments.
 * @param str [string] // AA sequence
 */
   string Filter(string str);

};
 


/* Constructors */
template <typename Tint>
SEG<Tint>::SEG(){
  unordered_map<string,string> para;
  SetParamaters(para);
  alpha= MakeAlpha();
}; 

template <typename Tint>
template <typename Targ>
SEG<Tint>::SEG(Targ& arg){
  SetParamaters(arg);
  alpha= MakeAlpha();
}
/* Explicite missing*/



/* Destructors */
template <typename Tint>
SEG<Tint>::~SEG(){
    AlphaFree(alpha);
}
/*Explicite missing*/


/* Functions  : Public */

template <typename Tint>
string SEG<Tint>::Filter(string str){

  
  CSeq* seq;
  SeqSeg* segs;
  bool params_allocated = false;
  Tint status = 0;
  SeqSeg  *seg;
  Tint begin, end, i;

/* old schoole - parse */

  seq = NewCSeq();
  seq->seq = &str[0];
  seq->length = str.size();
  seq->palpha = alpha;
  
  segs = (SeqSeg*) NULL;

/* compute lc segments */
   status = SegSeq (seq, &segs, 0);
   if (status < 0){
     seq->seq = NULL;
     CSeqFree (seq);
     throw runtime_error ("Low complexity segment computation could not be preformed!" ); 
   }
   
/* Create filtered sequence*/
    string filtstr(seq->seq);

/* merge segment if specified here you can completly omitt
 * this if raw positions are required by def is set to 1
 * - further testing required - 
 */
  if (MergeOverlaps == 1){
      MergeSegs(seq, segs);
      
/* Mask: X it*/
    for (seg=segs; seg!=NULL; seg=seg->next) {
      begin = seg->begin;
      end = seg->end;
      memset(&filtstr[0] + begin, 'X', end - begin +1);
    }
}

/* clean up */
   seq->seq = NULL;
   CSeqFree(seq);
   SegFree(segs);
   
   return filtstr;	
}


/* Functions  : Private */


template <typename Tint>
template <typename Targ>
void SEG<Tint>::SetParamaters(Targ &arg){
  
  SegWindow = (arg.find("window")   == arg.end() || StringToNumeric<Tint>(arg["window"])     <= 0)  ? 12   : StringToNumeric<Tint>(arg["window"]);
  SegHicut = (arg.find("hicut")     == arg.end() ||  StringToNumeric<double>(arg["hicut"])  < 0)  ?  2.5 : StringToNumeric<double>(arg["hicut"]);
  SegLocut = (arg.find("locut")     == arg.end() || StringToNumeric<double>(arg["locut"])   < 0)  ?  2.2 : StringToNumeric<double>(arg["locut"]);
  MaxX = (arg.find("maxXes")        == arg.end() || StringToNumeric<Tint>(arg["maxXes"])     < 0)  ?  0   : StringToNumeric<Tint>(arg["maxXes"]);
  MaxTrim = (arg.find("maxtrim")    == arg.end() || StringToNumeric<Tint>(arg["maxtrim"])    < 0)  ?  100 : StringToNumeric<Tint>(arg["maxtrim"]);
  Period = (arg.find("period")      == arg.end() || StringToNumeric<Tint>(arg["period"])     < 1)  ?  1   : StringToNumeric<Tint>(arg["period"]);
  MergeOverlaps = (arg.find("merge")== arg.end() || StringToNumeric<Tint>(arg["merge"])      < 1)  ?  1   : StringToNumeric<Tint>(arg["merge"]);
  
  
   if ( SegLocut > SegHicut)
       SegHicut = SegLocut;

   if (MaxX > SegWindow)
       MaxX = SegWindow;
}

template <typename Tint>
void SEG<Tint>::SegFree(SeqSeg* seg){
   SeqSeg* nextseg;
   while (seg){
      nextseg = seg->next;
      SafeFree(seg);
      seg = nextseg;
   }
}

template <typename Tint>
 void SEG<Tint>::MergeSegs(CSeq* seq, SeqSeg* segs){
	 
   SeqSeg* seg,* nextseg;          

   if (segs==NULL) return;

   if ((seq->length - 1 - segs->end) < 0) 
       segs->end = seq->length -1;

   seg = segs;
   nextseg = seg->next;

   while (nextseg!=NULL) {
      if (seg->begin - nextseg->end - 1 < 0) {
         if ((seg->end) < (nextseg->end)) seg->end = nextseg->end;
         if ((seg->begin) > (nextseg->begin)) seg->begin = nextseg->begin;
         seg->next = nextseg->next;
         SafeFree(nextseg);
      } else {
         seg = nextseg;
      }
      nextseg = seg->next;
   }
   if ((seg->begin) < 0) seg->begin = 0;
}

template <typename Tint>
typename SEG<Tint>::CSeq*  SEG<Tint>::NewCSeq(){
   CSeq* seq;

   seq = (CSeq*) calloc(1, sizeof(CSeq));
   if (seq==NULL)
      return(seq);

   seq->parent = (CSeq*) NULL;
   seq->seq = (char*) NULL;
   seq->palpha = (Alphabet*) NULL;
   seq->start = 0;
   seq->length = 0;
   seq->Xes =false;
   seq->charfreq = seq->state = (Tint*) NULL;
   seq->entropy = (double) 0.0;

   return(seq);
}

template <typename Tint>
double* SEG<Tint>::ComputeEntropy(CSeq* seq,Tint  first, Tint last, Tint downset, Tint upset){
	
   CSeq* win;
   double* H;
   Tint i;

   if (SegWindow>seq->length)
      return((double*) NULL);
     

   H = (double*) calloc(seq->length, sizeof(double));

   for (i=0; i<seq->length; i++)
      H[i] = -1.0;

   win = OpenWin(seq, 0, SegWindow);
   EntropyOn(win);

   for (i=first; i<=last; i++){
      if ((win->Xes )> MaxX){
         H[i] = -1.;
         ShiftWin1(win);
         continue;
        }
      H[i] = win->entropy;
      ShiftWin1(win);
     }

   CloseWin(win);
   return(H);
}

template <typename Tint>
Tint SEG<Tint>::SegSeq(CSeq* seq, SeqSeg **segs, Tint offset){
   SeqSeg* seg = (SeqSeg*) NULL;

   Tint downset, upset;
   Tint first, last, lowlim;
   Tint i;
   Tint leftend, rightend;
   double* H;
   Tint status = 0;

   downset = (SegWindow+1)/2 - 1;
   upset = SegWindow - downset;
   first = downset;
   last = seq->length - upset;
   lowlim = first;
   
   H = ComputeEntropy(seq,  first,  last,  downset,  upset);

   if (H == NULL) 
      return status;
   for (i=first; i<=last; i++){
      if (H[i] <= SegLocut && H[i] != -1.0){
         Tint loi = LocLow(i, lowlim, H); 
         Tint hii = LocHigh(i, last, H);
         CSeq* temp_seq = NULL;

         leftend = loi - downset;
         rightend = hii + upset - 1;

         temp_seq = OpenWin(seq, leftend, rightend-leftend+1);
         status = Trim(temp_seq, &leftend, &rightend); 

         if (status < 0) {
             CloseWin(temp_seq);
             break;
         }

         if (i+upset-1<leftend){
            Tint lend = loi - downset;
            Tint rend = leftend - 1;

            CSeq* leftseq = OpenWin(seq, lend, rend-lend+1);
            SeqSeg *leftsegs = (SeqSeg*) NULL;
            status = SegSeq(leftseq,  &leftsegs, offset+lend);
            if (status < 0)
              return status;

            if (leftsegs!=NULL){
               leftsegs->next = *segs;
               *segs = leftsegs;
            }
            CloseWin(leftseq);
         }

         seg = (SeqSeg*) calloc(1, sizeof(SeqSeg));
         seg->begin = leftend + offset;
         seg->end = rightend + offset;
         seg->next = *segs;
         *segs = seg;
         i = min(hii, rightend+downset);
         lowlim = i + 1;
        }
   }
   SafeFree(H);
   return status;
}

template <typename Tint>
Tint SEG<Tint>::LocLow(Tint i, Tint limit, double* H){
   Tint j;
   for (j=i; j>=limit; j--){
      if (H[j]==-1.0) break;
      if (H[j]>SegHicut) break;
   }
   return(j+1);
}

template <typename Tint>
Tint SEG<Tint>::LocHigh(Tint i, Tint limit, double* H){
   Tint j;
   for (j=i; j<=limit; j++){
      if (H[j]==-1.0) break;
      if (H[j]>SegHicut) break;
   }
   return(j-1);
}

template <typename Tint>
typename SEG<Tint>::CSeq* SEG<Tint>::OpenWin(CSeq* parent, Tint start, Tint length){
   CSeq* win;

   if (start<0 || length<0 || start+length>parent->length)
      return((CSeq*) NULL);

    win = (CSeq*) calloc(1, sizeof(CSeq));
    win->parent = parent;
    win->palpha = parent->palpha;
    win->start = start;
    win->length = length;
    win->seq = parent->seq + start;
    win->Xes = 0;
    win->entropy = -2.;
    win->state = (Tint*) NULL;
    win->charfreq = (Tint*) NULL;
	
    StateOn(win);

    return win;
}

template <typename Tint>
void SEG<Tint>::StateOn(CSeq* win){
	Tint letter, nel, c;
    Tint alphasize =  win->palpha->alphasize;

	if (win->charfreq == NULL)
		CompOn(win);
	win->state = (Tint*) calloc((alphasize+1), sizeof(win->state[0]));
	for (letter = nel = 0; letter < alphasize; ++letter) {
		if ((c = win->charfreq[letter]) == 0)
			continue;
		win->state[nel++] = c;
	}
	for (letter = nel; letter < alphasize+1; ++letter)
		win->state[letter] = 0;

    sort(win->state,  win->state+nel,greater<Tint>());
}

template <typename Tint>
void SEG<Tint>::CompOn(CSeq* win){
  Tint* comp;
  Tint letter;
  char* seq = win->seq; 
  char* seqmax = seq + win->length;
  Tint* alphaindex = win->palpha->alphaindex;
  unsigned char* alphaflag = win->palpha->alphaflag;
  Tint alphasize = win->palpha->alphasize;

  win->charfreq = comp = (Tint*) calloc(alphasize, sizeof(Tint));

  while (seq < seqmax) {
    letter = *seq++;
    if (!alphaflag[letter])
      comp[alphaindex[letter]]++;
    else 
      win->Xes++;
  }
}

template <typename Tint>
Tint SEG<Tint>::Trim(CSeq* seq, Tint* leftend, Tint* rightend){
   double prob, minprob = 1;
   Tint len;
   Tint lend=0 ;
   Tint rend =seq->length - 1;
   Tint minlen = 1;
   Tint status = 0;

   if ((seq->length-MaxTrim)>minlen) 
        minlen = seq->length-MaxTrim;
  
   for (len=seq->length; len>minlen; len--){
      bool shift = true;
      Tint i = 0;
      CSeq* win = OpenWin(seq, 0, len);

      while (shift){
         prob = GetProb(win->state, len, win->palpha);
         if (prob<minprob)
         {
            minprob = prob;
            lend = i;
            rend = len + i - 1;
         }
         shift = ShiftWin1(win);
         i++;
      }
      CloseWin(win);
   }

   *leftend = *leftend + lend;
   *rightend = *rightend - (seq->length - rend - 1);

   CloseWin(seq);
   return status;
}

template <typename Tint>
void SEG<Tint>::CloseWin(CSeq* win)
{
   if (win==NULL) return;

   if (win->state!=NULL)       SafeFree(win->state);
   if (win->charfreq!=NULL) SafeFree(win->charfreq);

   SafeFree(win);
   return;
}

template <typename Tint>
bool SEG<Tint>::ShiftWin1(CSeq* win){
	
  Tint j, length = win->length;
  Tint* comp = win->charfreq;
  Tint* alphaindex = win->palpha->alphaindex;
  unsigned char* alphaflag = win->palpha->alphaflag;

  if ((++win->start + length) > win->parent->length) {
    --win->start;
    return false;
  }

  if (!alphaflag[j = win->seq[0]])
    DecrementSV(win->state, comp[alphaindex[j]]--);
  else 
    win->Xes--;
  j = win->seq[length];
  ++win->seq;

  if (!alphaflag[j])
    IncrementSV(win->state, comp[alphaindex[j]]++);
  else 
    win->Xes++;

  if (win->entropy > -2.)
    win->entropy = Entropy(win->state);

	return true;
}

template <typename Tint>
double SEG<Tint>::Entropy(Tint* sv){
   double ent;
   Tint i, total = 0;

   for (i=0; sv[i]!=0; i++)
     {
      total += sv[i];
     }
   if (total==0) return(0.);

   ent = 0.0;
   const double ln2 = 0.69314718055994530941723212145818;
   for (i=0; sv[i]!=0; i++)
      ent += ((double)sv[i])*log(((double)sv[i])/(double)total)/ln2;   /// WHAT IS THIS
   double f = ent/(double)total; 
   ent = (f < 0.0) ? -(f) : (f);

   return(ent);
}

template <typename Tint>
void SEG<Tint>::DecrementSV(Tint* sv, Tint clas){
  Tint	svi;
  while ((svi = *sv++) != 0) {
    if (svi == clas && *sv < clas) {
      sv[-1] = svi - 1;
    break;
    }
  }
}

template <typename Tint>
void SEG<Tint>::IncrementSV(Tint* sv, Tint clas){
  for (;;) {
    if (*sv++ == clas) {
      sv[-1]++;
      break;
    }
  }
}

template <typename Tint>
double SEG<Tint>::GetProb(const Tint* sv, Tint total, const Alphabet* palpha){
   double  ans1, ans2 = 0, totseq;

   totseq = ((double) total) * (palpha->lnalphasize);

   ans1 = LnAss(sv, palpha->alphasize);
   if (ans1 > -100000.0 && sv[0] != (-2147483647-1))  // ncbi constant
      ans2 = LnPerm(sv, total);
   else
      cerr <<"Bad value!\n";
   
   return  ans1 + ans2 - totseq;
}
  
template <typename Tint>
double SEG<Tint>::LnPerm(const Tint* sv, Tint window_length){
   double ans;
   Tint i;
   
   ans = lnFact(window_length);
   for (i=0; sv[i]!=0; i++)
      ans -= lnFact(sv[i]);

   return(ans);
}

template <typename Tint>
double SEG<Tint>::LnAss(const Tint* sv, Tint alphasize){
  double	ans;
  Tint	svi, svim1;
  Tint	clas, total;
  Tint    i;

  ans = lnFactA[alphasize];
  if (sv[0] == 0)
    return ans;

  total = alphasize;
  clas = 1;
  svi = *sv;
  svim1 = sv[0];
  for (i=0;; svim1 = svi) {
    if (++i==alphasize) {
      ans -= lnFact(clas);
      break;
    }else if ((svi = *++sv) == svim1) {
      clas++;
      continue;
    }else {
      total -= clas;
      ans -= lnFact(clas);
      if (svi == 0) {
        ans -= lnFact(total);
        break;
      }else {
        clas = 1;
        continue;
      }
    }
  }
  return ans;
}
 
template <typename Tint>
double SEG<Tint>::lnFact(Tint n) {
  if (n < sizeof(lnFactA)/sizeof(*lnFactA))
     return lnFactA[n];
  else 
    return ((n+0.5)*log(n) - n + 0.9189385332);
}

template <typename Tint>
void SEG<Tint>::CSeqFree(CSeq* seq){
   if (seq==NULL) return;
   SafeFree(seq->seq);
   SafeFree(seq->charfreq);
   SafeFree(seq->state);
   SafeFree(seq);
}

template <typename Tint>
typename SEG<Tint>::Alphabet* SEG<Tint>::MakeAlpha (){
   Alphabet* palpha;
   Tint* alphaindex;
   unsigned char* alphaflag;
   Tint c,  i;
   const int kCharSet = 128;
   const double kLn20 = 2.9957322735539909;  // ncbi 

   palpha = (Alphabet*) calloc(1, sizeof(Alphabet));

   palpha->alphasize = 20;
   palpha->lnalphasize = kLn20;

   alphaindex = (Tint*) calloc(kCharSet , sizeof(Tint));
   alphaflag = (unsigned char*) calloc(kCharSet , sizeof(char));

   for (c=0, i=0; c<kCharSet; c++)
     {
        if (c == 65 || 
           (c >= 67 && c <= 73) || 
           (c >= 75 && c <= 78) || 
           (c >= 80 && c <= 84) || 
           (c >= 86 && c <= 87) ||  
           c == 89) {
			   
           alphaflag[c] = false; 
           alphaflag[c+32] = false;
           alphaindex[c] = i; 
           alphaindex[c+32] = i;
           ++i;
        } else {
           alphaflag[c] = true; 
           alphaindex[c] = 20;
        }
     }

   palpha->alphaindex = alphaindex;
   palpha->alphaflag = alphaflag;

   return (palpha);
}

template <typename Tint>
void SEG<Tint>::AlphaFree (Alphabet* alpha){
   SafeFree (alpha->alphaindex);
   SafeFree (alpha->alphaflag);
   SafeFree (alpha);
}

template <typename Tint>
void  SEG<Tint>::EntropyOn(CSeq* win){
   if (win->state==NULL) 
     StateOn(win);
   win->entropy = Entropy(win->state);
}

template <typename Tint>
template <typename T>
void SEG<Tint>::SafeFree(T **x){
    free(*x);
    *x = NULL;

}

template <typename Tint>
template <typename T>
void SEG<Tint>::SafeFree(T *x){
    free(x);
    x = NULL;

}


} // end fasta
