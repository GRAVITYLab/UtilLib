#ifndef _IOUTILS_H  
#define _IOUTILS_H  

#include <iostream>
#include <cstring>
#include <cstdio>
#include <list>

using namespace std;

namespace utillib
{
	class IOUtils
	{
	private:

	public:
		static long int _ulib_io_getFileSize( FILE *filePtr );
		static long int _ulib_io_getFileSize( string s );
		static bool _ulib_io_feof ( FILE* f );

		static void _ulib_io_reverseByteOrder( string infilename, string outfilename );  
		static void _ulib_io_getArgs( char *argsFileName, list<string> *argNames, list<string> *argValues );
		static char* _ulib_io_getArgWithName( char* name, list<string> *argNames, list<string> *argValues );

	};
}

#endif
