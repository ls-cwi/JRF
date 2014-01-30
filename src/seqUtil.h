#ifndef __SEQUTIL_H__
#define __SEQUTIL_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define APPEND_LEN	256

//#ifdef __SEQUTIL_C__
//	void seqMemInit();
//	void* seqMalloc(int size);
//	void seqFree();
//	void seqFreeAll();
// void inputString(char *input, char **ppcStr, int *iLen, int *iMaxLen);
//#else

  void seqMemInit();
	 void* seqMalloc(int size);
	 void seqFreeAll();
	 void seqFree(void* pos);
	 void inputString(char *input, char **ppcStr, int *iLen, int *iMaxLen);
//#endif



#endif

