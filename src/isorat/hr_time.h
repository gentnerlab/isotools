/* $Id: hr_time.h,v 1.1 2011/07/31 18:26:17 samn Exp $ */
// from http://allmybrain.com/2008/06/10/timing-cc-code-on-linux/ 
#ifndef __HR_TIME_H
#define __HR_TIME_H

#ifdef WIN32
#include <windows.h>

typedef struct {
    LARGE_INTEGER start;
    LARGE_INTEGER stop;
} stopWatch;

class CStopWatch {

private:
	stopWatch timer;
	LARGE_INTEGER frequency;
	double LIToSecs( LARGE_INTEGER & L);
public:
	CStopWatch();
	void startTimer( );
	void stopTimer( );
	double getElapsedTime();
};

#else
#include <sys/time.h>

typedef struct {
	timeval start;
	timeval stop;
} stopWatch;

class CStopWatch {

private:
	stopWatch timer;
public:
	CStopWatch() {};
	void startTimer( );
	void stopTimer( );
	double getElapsedTime();
};

#endif

#endif
