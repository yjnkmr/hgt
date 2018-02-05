// Definition of Functions for 'SEQ.h'

#include <string>
#include <vector>
#include "SEQ.h"

using namespace std;


string SEQ::reverse(string &str){
	string seq = "";
	int i;
	for (i = 0; i < str.size(); i++){
		seq = str.at(i) + seq;
	}
	return (seq);
}


// Complement DNA Sequence
// Note: This just outputs complement sequences. Therefore, for example, "AGGT" outputs "TCCA"
// The direction of DNA starnd (5'-, 3'-) is NOT changed !!!
// If you want biologically meaningful complement sequence, use "SEQ::reverse_complement"
string SEQ::complement(string &str){
	string seq = "";
	int i;
	for (i = 0; i < str.size(); i++){
		if      (str.at(i) == 'T'){ seq = seq + "A" ; }
		else if (str.at(i) == 'C'){ seq = seq + "G" ; }
		else if (str.at(i) == 'A'){ seq = seq + "T" ; }
		else if (str.at(i) == 'G'){ seq = seq + "C" ; }
		// Irregular Bases
		else if (str.at(i) == 'N'){ seq = seq + "N" ; }
		else if (str.at(i) == 'R'){ seq = seq + "Y" ; }
		else if (str.at(i) == 'Y'){ seq = seq + "R" ; }
		else if (str.at(i) == 'M'){ seq = seq + "K" ; }
		else if (str.at(i) == 'K'){ seq = seq + "M" ; }
		else if (str.at(i) == 'S'){ seq = seq + "S" ; }
		else if (str.at(i) == 'W'){ seq = seq + "W" ; }
		else if (str.at(i) == 'H'){ seq = seq + "D" ; }
		else if (str.at(i) == 'D'){ seq = seq + "H" ; }
		else if (str.at(i) == 'B'){ seq = seq + "V" ; }
		else if (str.at(i) == 'V'){ seq = seq + "B" ; }
	}
	return (seq);
}

string SEQ::reverse_complement(string &str){
	string seq = "";
	int i;
	for (i = 0; i < str.size(); i++){
		if      (str.at(i) == 'T'){ seq = "A" + seq ; }
		else if (str.at(i) == 'C'){ seq = "G" + seq ; }
		else if (str.at(i) == 'A'){ seq = "T" + seq ; }
		else if (str.at(i) == 'G'){ seq = "C" + seq ; }
		// Irregular Bases
		else if (str.at(i) == 'N'){ seq = "N" + seq ; }
		else if (str.at(i) == 'R'){ seq = "Y" + seq ; }
		else if (str.at(i) == 'Y'){ seq = "R" + seq ; }
		else if (str.at(i) == 'M'){ seq = "K" + seq ; }
		else if (str.at(i) == 'K'){ seq = "M" + seq ; }
		else if (str.at(i) == 'S'){ seq = "S" + seq ; }
		else if (str.at(i) == 'W'){ seq = "W" + seq ; }
		else if (str.at(i) == 'H'){ seq = "D" + seq ; }
		else if (str.at(i) == 'D'){ seq = "H" + seq ; }
		else if (str.at(i) == 'B'){ seq = "V" + seq ; }
		else if (str.at(i) == 'V'){ seq = "B" + seq ; }
	}
	return (seq);
}


// Encode DNA sequence into 10 digit Number
int SEQ::encode(string str){
	int encode = 0;
	int digit, i;
	string n;
	for (i = 0; i < str.size(); i++){
		n = str.at(i);
		if      (n == "T"){digit = 0;}
		else if (n == "C"){digit = 1;}
		else if (n == "A"){digit = 2;}
		else if (n == "G"){digit = 3;}
		else    {return (-1);}
		encode = encode * 4 + digit;
	}
	return (encode);
}


// Decode 10 digit number into DNA sequence (sequence length == len)
string SEQ::decode(int code, int len){
	int digit, i;
	string seq = "";
	for (i = 0; i < len; i++){
		digit = code % 4;
		if      (digit == 0){seq = "T" + seq;}
		else if (digit == 1){seq = "C" + seq;}
		else if (digit == 2){seq = "A" + seq;}
		else if (digit == 3){seq = "G" + seq;}
		code = code >> 2;
	}
	return (seq);
}

// Split Sequence
// If i == 3 and slide == 3, for example, sequence is split by codon triplet
vector<string> SEQ::split ( string& str , int tuple , int  slide ) {
	vector<string> seq_vec ;
	//if ( str.size() % slide  !=  0 ) { return ( seq_vec ) ; }
	int k ; 
	for ( k = 0 ; k + tuple - 1 < str.size() ; k = k + slide ) {
		seq_vec.push_back ( str.substr ( k, tuple ) ) ;
	}
	return ( seq_vec ) ;
}


string SEQ::join ( vector<string>& seq_vec ) {
	string str = "" ;
	int i ;
	for ( i = 0 ; i < seq_vec.size() ; i++ ) {
		str += seq_vec[i] ;
	}
	return ( str ) ;
}


vector<string>  SEQ::getSNP ( string str ) {
	int i, k ;
	string nuc ;
	vector<string>  snp_vec ;
	snp_vec.resize(0) ;
	for ( i = 0 ; i < str.size() ; i++ ) {
		if (str.at(i) == 'T'){
			nuc = str.substr(0, i) + "C" + str.substr(i+1) ;
			snp_vec.push_back ( nuc ) ;
			nuc = str.substr(0, i) + "A" + str.substr(i+1) ;
			snp_vec.push_back ( nuc ) ;
			nuc = str.substr(0, i) + "G" + str.substr(i+1) ;
			snp_vec.push_back ( nuc ) ;
		} else if (str.at(i) == 'C'){
			nuc = str.substr(0, i) + "T" + str.substr(i+1) ;
			snp_vec.push_back ( nuc ) ;
			nuc = str.substr(0, i) + "A" + str.substr(i+1) ;
			snp_vec.push_back ( nuc ) ;
			nuc = str.substr(0, i) + "G" + str.substr(i+1) ;
			snp_vec.push_back ( nuc ) ;
		} else if (str.at(i) == 'A'){
			nuc = str.substr(0, i) + "T" + str.substr(i+1) ;
			snp_vec.push_back ( nuc ) ;
			nuc = str.substr(0, i) + "C" + str.substr(i+1) ;
			snp_vec.push_back ( nuc ) ;
			nuc = str.substr(0, i) + "G" + str.substr(i+1) ;
			snp_vec.push_back ( nuc ) ;
		} else if (str.at(i) == 'G'){
			nuc = str.substr(0, i) + "T" + str.substr(i+1) ;
			snp_vec.push_back ( nuc ) ;
			nuc = str.substr(0, i) + "C" + str.substr(i+1) ;
			snp_vec.push_back ( nuc ) ;
			nuc = str.substr(0, i) + "A" + str.substr(i+1) ;
			snp_vec.push_back ( nuc ) ;
		}
	}
	return ( snp_vec ) ;
}

