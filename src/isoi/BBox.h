/* $Id: BBox.h,v 1.1 2011/01/06 22:29:07 samn Exp $ */
#ifndef BBOX_H
#define BBOX_H

//struct to represent bounding box
//in a coord. system where top is 0
//and bottom of image is h

template <class T>
struct BBox
{

	//data
	T left,top,right,bottom;

	//functions
	BBox()
	{
		left=top=right=bottom=T(0);
	}
	BBox(const T& left_,const T& top_,const T& right_,const T& bottom_)
	{
		left = left_;
		top = top_;
		right = right_;
		bottom = bottom_;
	}
	template<class T2>
	BBox(const BBox<T2>& B)
	{
		left = B.left;
		top = B.top;
		right = B.right;
		bottom = B.bottom;
	}
	int area() const
	{
		return ( right - left + 1 ) * ( bottom - top + 1 );
	}
	int width() const
	{
		return right - left + 1;
	}
	//assumes bottom has higher value than top
	int height() const
	{
		return bottom - top + 1;
	}
	bool operator==(const BBox& oBox) const
	{
		return left==oBox.left &&
			   right==oBox.right &&
			   top==oBox.top &&
			   bottom==oBox.bottom;
	}
	bool operator!=(const BBox& oBox)
	{
		return !(*this == oBox);
	}
	BBox& operator/=(int i)
	{
		if(i)
		{
			left /= i;
			top /= i;
			right /= i;
			bottom /= i;
		}
		return *this;
	}
	BBox operator/(int i)
	{
		BBox tmp = *this;
		if(i)
		{
			tmp /= i;
		}
		return tmp;
	}
	//make sure coords are valid
	bool valid() const
	{
		return left >= 0 && right >= 0 &&
			   top >= 0 && bottom >= 0 &&
			   left <= right && bottom >= top;
	}

	template<class T2>
	BBox<T>& operator=(const BBox<T2>& B2)
	{
		left = B2.left;
		right = B2.right;
		top = B2.top;
		bottom = B2.bottom;
		return *this;
	}

	T MidX()
	{
		return left + width() / 2;
	}

	T MidY()
	{
		return top + height() / 2;
	}
};

//true iff box1 above box2
template <class T1, class T2>
inline bool Above(const BBox<T1>& box1,const BBox<T2>& box2)
{
	return box1.bottom < box2.top;
}

//true iff box1 below box2
template <class T1, class T2>
inline bool Below(const BBox<T1>& box1,const BBox<T2>& box2)
{
	return box1.top > box2.bottom;
}

template <class T1, class T2>
inline bool Left(const BBox<T1>& box1,const BBox<T2>& box2)
{
	return box1.right < box2.left;
}

template <class T1, class T2>
inline bool Right(const BBox<T1>& box1,const BBox<T2>& box2)
{
	return box1.left > box2.right;
}

template <class T>
inline int MidHeight(BBox<T>& box)
{
	return box.top + ( box.bottom - box.top ) / 2;
}


//do the two bounding boxes intersect?
template <class T1,class T2>
inline bool Intersect(const BBox<T1>& box1, const BBox<T2>& box2)
{
	return !( (box1.top		> box2.bottom) ||
			  (box1.bottom	< box2.top) ||
			  (box1.right	< box2.left) ||
			  (box1.left	> box2.right));
}

//get intersection of bounding boxes , assumes they intersect
template <class T>
inline BBox<T> GetIntersection(const BBox<T>& b1, const BBox<T>& b2)
{
	BBox<T> ret;
	ret.left = b1.left > b2.left ? b1.left : b2.left;
	ret.right = b1.right < b2.right ? b1.right : b2.right;
	ret.top = b1.top > b2.top ? b1.top : b2.top;
	ret.bottom = b1.bottom < b2.bottom ? b1.bottom : b2.bottom;
	return ret;
}

//get union of bounding boxes
template <class T>
inline BBox<T> GetUnion(const BBox<T>& b1, const BBox<T>& b2)
{
	BBox<T> ret;
	ret.left = b1.left < b2.left ? b1.left : b2.left;
	ret.right = b1.right > b2.right ? b1.right : b2.right;
	ret.top = b1.top < b2.top ? b1.top : b2.top;
	ret.bottom = b1.bottom > b2.bottom ? b1.bottom : b2.bottom;
	return ret;
}

template <class T>
inline bool Zero(BBox<T>& oB)
{
	return !(oB.left || oB.right ||
		     oB.top || oB.bottom);
}

typedef BBox<int> IBBox;
typedef BBox<unsigned int> UBBox;

struct BBoxD
{
	int left, top, right, bottom, rotate;

	BBoxD()
	{
		left = top = right = bottom = 1;
		rotate = 0;
	}

	BBoxD(int l, int t, int r, int b)
	{
		left = l;
		top = t;
		right = r;
		bottom = b;
		rotate = 0;
	}

	BBoxD(int l, int t, int r, int b, int rot)
	{
		left = l;
		top = t;
		right = r;
		bottom = b;
		rotate = rot;
	}
};

#endif
