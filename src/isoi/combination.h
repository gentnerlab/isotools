// $Id: combination.h,v 1.1 2011/01/07 03:42:28 samn Exp $ 
// Combination algorithm implementation
//
// Copyright (C) 2004 - 2006, BenBear
//
// More information in http://www.bxmy.org

// This file is an algorithm of the combination. This library is free
// software; you can redistribute it and/or modify it under the terms
// of the GNU General Public License as published by the Free Software
// Foundation; either version 2, or (at your option) any later
// version.

// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this library; see the file COPYING.  If not, write to
// the Free Software Foundation, 59 Temple Place - Suite 330, Boston,
// MA 02111-1307, USA.

#ifndef __btb_combination_hpp_def
#define __btb_combination_hpp_def

#include <algorithm>

namespace btb
{
  ////////////////////////////////////////////////////////////////////
  // combination for STL permutation
  // 
  // combination_init    (first, middle, last)
  // combination_adjust  (first, middle, last)
  // combination         (first1, last1, first2, last2)
  // next_combination    (first, middle, last)
  // prev_combination    (first, middle, last)
  ////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////
  // twice_merge: merge [first1, last1) and [first2, last2), then fill
  //              the min element to [first1, last1), left to [first2,
  //              last2)
  ////////////////////////////////////////////////////////////////////
  template <typename Iter>
  void
  twice_merge (Iter first1, Iter last1, Iter first2, Iter last2)
  {
    typedef typename std::iterator_traits<Iter>::value_type value_type;
    typedef typename std::iterator_traits<Iter>::difference_type diff_t;

    if ((first1 == last1) || (first2 == last2))
      return;
    first1 = std::upper_bound(first1, last1, *first2);
    if (first1 == last1)
      return;
    last2 = std::lower_bound(first2, last2, *(last1-1));
    if (first2 == last2)
      return;

    diff_t len1 = std::distance(first1, last1);
    diff_t len2 = std::distance(first2, last2);
    bool min1 = len1 < len2;
    diff_t lent = min1 ? len1 : len2;
    value_type *tmp = new value_type[lent];

    if (min1)
      {
	std::copy(first1, last1, tmp);
	Iter p = first2;
	Iter q = tmp;
	Iter i;
	for (i = first1; i != last1; ++i)
	  if (*p < *q)
	    *i = *p++;
	  else
	    *i = *q++;
	for (i = first2; (i != last2) && (p != last2); ++i)
	  if (*p < *q)
	    *i = *p++;
	  else
	    *i = *q++;
	for (; i != last2; ++i, ++q)
	  *i = *q;
      }
    else
      {
	std::copy(first2, last2, tmp);
	Iter p = last1;
	Iter q = tmp+lent;
	Iter i;
	for (/*Iter*/ i = last2; i != first2;)
	  if (*(p-1) < *(q-1))
	    *--i = *--q;
	  else
	    *--i = *--p;
	for (i = last1; (i != first1) && (p != first1);)
	  if (*(p-1) < *(q-1))
	    *--i = *--q;
	  else
	    *--i = *--p;
	for (; i != first1;)
	  *--i = *--q;
      }

    delete[] tmp;
  }

  ///////////////////////////////////////////////////////////////////
  // __combination_sort: merge sort the [first, last)
  ///////////////////////////////////////////////////////////////////
  template <typename Iter>
  void
  __combination_sort (Iter first, Iter last)
  {
    typedef typename std::iterator_traits<Iter>::difference_type diff_t;
    diff_t len = std::distance(first, last);
    if (len <= 1)
      return;

    if (len == 2)
      {
	if (*first > *--last)
	  std::iter_swap (first, last);
      }
    else
      {
	Iter middle = first;
	std::advance(middle, len / 2);
	__combination_sort (first, middle);
	__combination_sort (middle, last);
	twice_merge (first, middle, middle, last);
      }
  }

  //////////////////////////////////////////////////////////////////////
  // combination_init: init the (first, midle, last) to the min or the
  //                   max combination
  //////////////////////////////////////////////////////////////////////
  template <typename Iter>
  void
  combination_init (Iter first, Iter middle, Iter last, bool min = true)
  {
    __combination_sort (first, middle);
    __combination_sort (middle, last);
    if (min)
      twice_merge (first, middle, middle, last);
    else
      twice_merge (middle, last, first, middle);
  }

  //////////////////////////////////////////////////////////////////////
  // combination_adjust: make the (first, middle, last) to a right
  //                     combination. [first, middle) are the elements
  //                     selected in, [middle, last) are selected out
  //////////////////////////////////////////////////////////////////////
  template <typename Iter>
  void
  combination_adjust (Iter first, Iter middle, Iter last)
  {
    __combination_sort (first, middle);
    __combination_sort (middle, last);
  }

  /////////////////////////////////////////////////////////////////////
  // combination: get next combination.
  //
  // [first1, last1): the elements selected in \\
  // [first2, last2): the elements selected out
  /////////////////////////////////////////////////////////////////////
  template <typename Iter>
  bool
  combination (Iter first1, Iter last1, Iter first2, Iter last2)
  {
    if ((first1 == last1) || (first2 == last2))
      return false;

    Iter qmax = last2;
    --qmax;
    Iter pout1 = std::lower_bound(first1, last1, *qmax);
    bool fin = pout1 == first1;
    Iter left1, left2;
    if (!fin)
      {
	Iter pout = pout1;
	--pout;
	Iter qin = std::upper_bound(first2, last2, *pout);
	std::iter_swap (pout, qin);
	left1 = pout;
	++left1;
	left2 = qin;
	++left2;
      }
    else
      {
	left1 = first1;
	left2 = first2;
      }
    twice_merge (left1, last1, left2, last2);
    return !fin;
  }

  /////////////////////////////////////////////////////////////////////
  // next_combination: get next combination.
  //
  // [first, middle): the elements selected in \\
  // [middle, last): the elements selected out
  /////////////////////////////////////////////////////////////////////
  template <typename Iter>
  inline bool
  next_combination (Iter first, Iter middle, Iter last)
  {
    return combination (first, middle, middle, last);
  }

  /////////////////////////////////////////////////////////////////////
  // prev_combination: get prev combination.
  //
  // [first, middle): the elements selected in \\
  // [middle, last): the elements selected out
  /////////////////////////////////////////////////////////////////////
  template <typename Iter>
  inline bool 
  prev_combination (Iter first, Iter middle, Iter last)
  {
    return combination (first, middle, middle, last);
  }
}

#endif
