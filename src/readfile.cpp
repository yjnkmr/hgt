// Input from FASTA format file


#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <map>
#include <algorithm>

#include "readfile.h"

using namespace std;


////////////////////////////////////////////////////////////////////////////
/*
ReadFASTA::ReadFASTA(string f){
		filename = f;					// Input file name
	}
ReadFASTA::ReadFASTA(){
		cout << "File name is needed." << endl;
		exit (0);
	}
*/


// Member function
int  ReadFile::readDataLine ( string fn, vector<string>& vsd ) {
	
	ifstream file_in ;
	file_in.open ( fn.c_str() ) ;
	if ( ! file_in ) {
		cerr << fn << " doesn't exist." << endl ;
		return 0 ;
	}

	string string_datum ;
	max_line_len = 0 ;
	char buffer[ BUFFER_SIZE ] ;
	while ( ! file_in.eof() ) {
		
		file_in.getline ( buffer , BUFFER_SIZE ) ;
		string_datum = buffer ;

		if ( string_datum.size() == 0 ) { continue ; }		// Skip the blank lines
		if ( string_datum.at(0) == '#' ) { continue ; }		// Skip the sentences having "#" at the first position
		if ( string_datum.find("//") == 0 ) { continue ; }		// Skip the sentences having "//" at the first position

		if ( string_datum.size() > max_line_len ) { max_line_len = string_datum.size() ; }
		vsd.push_back ( string_datum ) ;
	}
	file_in.close () ;
	if ( vsd.size() == 0 ) {
		cerr << fn << " has No data." << endl ;
		return 0 ;
	}
	return 1 ;
}


int  ReadFile::readDataLine ( string fn, vector<string>& vsd, vector<string>& com ) {
	
	ifstream file_in ;
	file_in.open ( fn.c_str() ) ;
	if ( ! file_in ) {
		cerr << fn << " doesn't exist." << endl ;
		return 0 ;
	}

	string string_datum ;
	max_line_len = 0 ;
	char buffer[ BUFFER_SIZE ] ;
	while ( ! file_in.eof() ) {
		
		file_in.getline ( buffer , BUFFER_SIZE ) ;
		string_datum = buffer ;

		if ( string_datum.size() == 0 ) { continue ; }		// Skip the blank lines
		if ( string_datum.at(0) == '#' ) { 
			com.push_back ( string_datum ) ;	// Comment lines
			continue ; 
		}		// Skip the sentences having "#" at the first position
		if ( string_datum.find("//") == 0 ) { continue ; }		// Skip the sentences having "//" at the first position

		if ( string_datum.size() > max_line_len ) { max_line_len = string_datum.size() ; }
		vsd.push_back ( string_datum ) ;
	}
	file_in.close () ;
	if ( vsd.size() == 0 ) {
		cerr << fn << " has No data." << endl ;
		return 0 ;
	}
	return 1 ;
}



int ReadFile::readFASTA( string fn ){
	filename = fn ;
	vector<string> vsd ;
	if ( ReadFile::readDataLine ( filename , vsd ) == 0 ) { 
		return 0 ; 
	}

	int i, pos ;
	map<string , string> temp_datum ;
	
	for ( i = 0 ; i < vsd.size() ; i++ ) {
		if ( vsd[i].at(0) == '>' ) {
			//cout << vsd[i] << endl ;
			if ( i > 0  &&  vsd[i-1].at(0) == '>' ) {
				cerr << filename << endl << vsd[i-1] << endl << "This has NO sequence." << endl ;
				return 0 ;
			}
			if ( vsd[i].size() < 2 ) { 
				cerr << filename << endl << vsd[i] << endl << "NO name" << endl ;   
				return 0 ; 
			}
			pos = vsd[i].find ( " " );			// find first 'space'
			if ( pos == 1 ) { 
				cerr << filename << endl << vsd[i] << endl << "Space between \">\" and SEQNAME" << endl ;   
				return 0 ; 
			}

			temp_datum[ "name" ] = vsd[i].substr ( 1, pos -1 ) ;
			while ( vsd[i].at ( vsd[i].size() - 1 ) == ' ' ) { vsd[i].erase( vsd[i].size() - 1 ) ; }
			//temp_datum[ "fullname" ] = vsd[i].substr ( 1, -1 );
			temp_datum[ "fullname" ] = vsd[i].substr ( 1 );

			temp_datum[ "seq" ] = "" ;
			seq_data.push_back ( temp_datum ) ;		// Increment of sequence data

		} else {
			if ( seq_data.size() == 0 ) { 
				cerr << filename << " : The first sequence has NO name."  << endl ;
				return 0 ; 
			}
			seq_data[ seq_data.size() - 1 ][ "seq" ] += vsd[i] ;
		}
	}
	if ( seq_data.size() == 0 ) { 
		cerr << filename << " : No FASTA Sequences." << endl ;
		return 0 ;
	}
	if ( seq_data[ seq_data.size() - 1 ][ "seq" ].size() == 0 ) { 
		cerr << filename << " : The last sequence has NO sequence." << endl ;
		return 0 ;
	}
	return 1 ;
}

void ReadFile::clearData() {
	seq_data.resize(0) ;
}

int  ReadFile::seq_count() {
	return ( seq_data.size() ) ;
}

int  ReadFile::max_line_length() {
	return ( max_line_len ) ;
}

void ReadFile::show () {
	int i ;
	for ( i = 0 ; i < seq_data.size() ; i++ ) {
		cout << "Name=" << seq_data[i][ "name" ] << "#" <<endl ;
		cout << "Fullname=" << seq_data[i][ "fullname" ] << "#" <<endl ;
		cout << "Sequence=" << seq_data[i][ "seq" ] << "#" << endl ;
	}
}


map<string, string> ReadFile::getSeqDatum ( int i ) {
	if ( i < 0  || i > seq_data.size() - 1 ) { 
		cerr << i << " is out of vector." << endl ;
		exit ( 0 ) ;
	}
	return ( seq_data[i] ) ;
}

map<string, string> ReadFile::getSeqDatumHash ( string& str ) {
	int i = hash_keys[ str ]  - 1 ;
	if ( i < 0 ) {
		cerr << str << " is NOT in the map." << endl ;
		exit ( 0 ) ;
	}
	return ( seq_data[i] ) ;
}

string  ReadFile::getSeq ( int i ) {
	if ( i < 0  || i > seq_data.size() - 1 ) { 
		cerr << i << " is out of vector." << endl ;
		exit ( 0 ) ;
	}
	return ( seq_data[i][ "seq" ] ) ;
}

string  ReadFile::getName( int i ) {
	if ( i < 0  || i > seq_data.size() - 1 ) { 
		cerr << i << " is out of vector." << endl ;
		exit ( 0 ) ;
	}
	return ( seq_data[i][ "name" ] ) ;
}


string  ReadFile::getFullname( int i ) {
	if ( i < 0  || i > seq_data.size() - 1 ) { 
		cerr << i << " is out of vector." << endl ;
		exit ( 0 ) ;
	}
	return ( seq_data[i][ "fullname" ] ) ;
}





void ReadFile::makeHash() {
	typedef  pair<string , int>  key_and_num ;
	int i ;
	for ( i = 0 ; i < seq_data.size() ; i++ ) {
		hash_keys.insert ( key_and_num( seq_data[i][ "name" ] , i+1 ) ) ;		// Make hash
	}
}

		
