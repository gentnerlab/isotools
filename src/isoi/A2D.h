// $Id: A2D.h,v 1.2 2011/01/08 01:12:12 samn Exp $ 
/*

isoi - This program calculates the Isolation Information (IsoI) cluster quality measures
described in the reference below. These measures were designed for clusters
obtained from neural extracellular recordings, but are applicable to an
arbitrary dataset.

Copyright (C) 2003-2011 Sam Neymotin & BioSignal Group

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

For more information contact Sam Neymotin ( samn at neurosim dot downstate dot edu )
 or Andre Fenton ( afenton at nyu dot edu ).

References:
 The methods used in this program are described in an article in press
  at The Journal of Neuroscience,
  Measuring the quality of neuronal identification in ensemble recordings
  by Neymotin SA, Lytton WW, Olypher AO, Fenton AA (2011).

*/

#ifndef A2D_H
#define A2D_H

template < typename T >
T **Allocate2DArray( int nRows, int nCols)
{
	if(nRows<=0 || nCols<=0) return 0x0;

    T **ppi;
    T *pool;
    T *curPtr;
    //(step 1) allocate memory for array of elements of column

    ppi = new T*[nRows];

    //(step 2) allocate memory for array of elements of each row
    pool = new T [nRows * nCols];

    // Now point the pointers in the right place
    curPtr = pool;
    for( int i = 0; i < nRows; i++)
    {
        *(ppi + i) = curPtr;
         curPtr += nCols;
    }
    return ppi;
}

template < typename T >
void Free2DArray(T** Array)
{
    delete [] *Array;
    delete [] Array;
}

template < class T >
class AutoFree2D
{
	T** m_p;
public:
	AutoFree2D(T** p)
		:m_p(p)
	{
	}
	~AutoFree2D()
	{
		Free2DArray(m_p);
	}
};

template < class T >
class A2D
{
	T** m_p;
	int m_iRows;
	int m_iCols;
	//A2D(const A2D&);

public:
	//get # rows
	int Rows() const { return m_iRows; }
	//get # cols
	int Cols() const { return m_iCols; }
	//construct empty
	A2D(T** p,int iRows,int iCols)
		:m_p(p),
		 m_iRows(iRows),
		 m_iCols(iCols)
	{
	}
	A2D()
		:m_p(0),
		m_iRows(0),
		m_iCols(0)
	{
	}
	A2D(const A2D& r)
	{	m_p=0;
		Init(r.Rows(),r.Cols());
		int i,j;
		for(i=0;i<m_iRows;i++)
			for(j=0;j<m_iCols;j++)
				m_p[i][j]=r[i][j];
	}
	A2D& operator=(const A2D& r)
	{	if(&r==this) return *this;
		Init(r.Rows(),r.Cols());
		int i,j;
		for(i=0;i<m_iRows;i++)
			for(j=0;j<m_iCols;j++)
				m_p[i][j]=r[i][j];
		return *this;
	}
	//init with specified size
	bool Init(int iRows,int iCols)
	{
		Clear();

		m_p = Allocate2DArray<T>(iRows,iCols);
		
		if(!m_p)
			return false;
		
		m_iRows = iRows;
		m_iCols = iCols;
		
		return true;
	}
	//construct with specified size
	A2D(int iRows,int iCols)
		:m_p(0)
	{
		Init(iRows,iCols);
	}
	//destructor
	virtual ~A2D()
	{
		Clear();
	}
	//free memory
	void Clear()
	{
		if(m_p)
			Free2DArray<T>(m_p);
		m_p=0;
		m_iRows=m_iCols=0;
	}
	//get 2D pointer
	T** GetP()
	{
		return m_p;
	}
	//get row pointer
	T* operator[](int i)
	{
		return m_p[i];
	}
	T* operator[](int i) const
	{
		return m_p[i];
	}
	//fill with val
	void Fill(const T& val)
	{
		int x,y;
		for(y=0;y<m_iRows;y++)
			for(x=0;x<m_iCols;x++)
				m_p[y][x]=val;
	}
};

#endif
