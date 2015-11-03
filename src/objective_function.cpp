/*
 * objective_function.cpp
 *
 *  Created on: Jan 19, 2014
 *      Author: ctsotskas
 */



#include "objective_function.h"

//
//int FileExists(char const *strFilename) { //// ORIGINAL was BOOL!!!!!
//	struct stat stFileInfo;
//	////bool blnReturn;
//	int intStat;
//
//	// Attempt to get the file attributes
//	intStat = stat(strFilename,&stFileInfo);
//	if(intStat == 0) {
//		// We were able to get the file attributes
//		// so the file obviously exists.
//		////	blnReturn = true; //// ORIGINAL
//		return 1;
//	} else {
//		// We were not able to get the file attributes.
//		// This may mean that we don't have permission to
//		// access the folder which contains this file. If you
//		// need to do that level of checking, lookup the
//		// return values of stat which will give you
//		// more details on why stat failed.
//		////	blnReturn = false; //// ORIGINAL
//		return 0;
//	}
//
//	//// return(blnReturn); //// ORIGINAL
//}
