#include "IOUtils.h"

namespace utillib
{
	void IOUtils::_ulib_io_getArgs( char *argsFileName, list<string> *argNames, list<string> *argValues )
	{
		char* argName = "";
		char* argVal = "";
		char nextLine[1000];

		// Open ascii file containing list of arguments
		FILE* argsFile = fopen( argsFileName, "r" );

		// Scan the file to create two lists of strings
		// List of argument names and list of argument values
		while( 1 )
		{
			// Read next line
			fgets( nextLine, 1000, argsFile );

			// Check for EOF
			if( feof( argsFile ) )
				break;

			// Skip to next line if this is a comment or a blank line
			if( *nextLine == '#' || *nextLine == '\n' )
				//printf( "Comment or blank line\n" );
	        		continue;

			// Tokenize string in to argument name and value
			argName = strtok( nextLine, " \n" );
			argVal = strtok( NULL, " \n" );
			string name( argName );
			string val( argVal );
	
			// Push strings in to lists
			argNames->push_back( name );
			argValues->push_back( val );
		}

		// Close file
		fclose( argsFile );

	}// end function

	char* IOUtils::_ulib_io_getArgWithName( char* name, list<string> *argNames, list<string> *argValues )
	{	
		list<string>::iterator iterName;
		list<string>::iterator iterVal;

		for( iterName = argNames->begin(), iterVal = argValues->begin(); iterName != argNames->end(); ++iterName, ++iterVal )
		{
			if( strcmp( (*iterName).c_str(), name ) == 0 )
				return (char*)( (*iterVal).c_str() );

		}// end for

		return "";

	}// end function

	long
	IOUtils::_ulib_io_getFileSize( FILE *filePtr )
	{
		if( filePtr == NULL )
			cout << "File Pointer Invalid" << endl;

		// get length of file:
		fseek( filePtr, 0, SEEK_END );
		long int flen = ftell( filePtr );
		
		// restore pointer back to the front of the file
		fseek( filePtr, 0, SEEK_SET );

		return flen;
	
	}// end function
	
	long
	IOUtils::_ulib_io_getFileSize( string s )
	{
		FILE* filePtr = fopen( s.c_str(), "rb" );
		if( filePtr == NULL )
		{
			fprintf( stderr, "%s: File Not found!\n", s.c_str() );
			return 0;
		}


		long fLen = IOUtils::_ulib_io_getFileSize( filePtr );
		fclose( filePtr );

		return fLen;
	}

	bool
	IOUtils::_ulib_io_feof ( FILE* f )
	{
	   int c = fgetc ( f );
	   if( c == EOF )
		   return true;

	   ungetc ( c, f );
	   return false;
	}


	void IOUtils::_ulib_io_reverseByteOrder( string infilename, string outfilename )
	{
		// Open input file 
		FILE *infile = fopen( infilename.c_str(), "rb" );
		FILE *outfile = fopen( outfilename.c_str(), "wb" );

		// Get input file size in bytes
		long int fileSize = _ulib_io_getFileSize( infile );
		cout << "File Size: " << fileSize << " bytes" << endl;

		// read and convert file
		char *first, *second, *third, *fourth;
		char *reverseArray;
		first = new char();
		second = new char();
		third = new char();
		fourth = new char();
		reverseArray = new char[4];

		long int i=0;
		while( i<fileSize )
		{
			// Read 4 bytes 1-by-1
			fread( first, 1, 1, infile );
			fread( second, 1, 1, infile );
			fread( third, 1, 1, infile );
			fread( fourth, 1, 1, infile );

			// Store them in an array in reverse order
			reverseArray[0] = *fourth;
			reverseArray[1] = *third;
			reverseArray[2] = *second;
			reverseArray[3] = *first;

			// write back the reversed array
			//fwrite( reverseArray, 1, 4, outfile );
			fwrite( ( float* )reverseArray, 4, 1, outfile );

			// increment i
			i += 4;
			//if( i%1000 == 0 ) cout << i << endl;

		}// end while

		// close files
		fclose( infile );
		fclose( outfile );

		delete first, second, third, fourth;
		delete [] reverseArray;

	}// end method

}// end namespace
