// Copyright: Yoji Nakamura, NRIFS, FRA
// Output is Counts or Log scores
// For usual pseudocount, use the option "-p 0" in execution
// For pseudocount ignoring termination codons, use the option "-p CODE" in execution
// For log scores, use the option "-l" in execution (after compilation)

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <map>

#include "readfile.h"
#include "SEQ.h"

using namespace std;

const  double  FLOOR_SCORE = -100.0 ;

void  seqCounterInitialization ( int, int, int ) ;
void  seqCount ( ReadFile&, int, int, int ) ;
int outTrans ( string&, string& , int& ) ;
int outEmit ( string&, string& , int& ) ;

// Global Variables
string  gCommand ;
vector<string>  gSeq_files ;
ReadFile  gRf ;
vector<int>  gCounter1D ;
vector<vector<int> >  gCounter2D ;
int  gTotal1D = 0 ;
vector<int>  gTotal2D ;
int gFrame ;
int gScore ; 
int gPseudo ; 

vector<int>  terFlag ( int , int ) ;
vector<int>  ter_flag ;
vector<string>  ter_list ; 


int main ( int argc , char *argv[] ) {

	if ( argc < 6 ) { 
		//cerr << "Arguments are not enough." << endl ;
		cerr << "Usage: ./PROGRAM STATE_ID CUT_LEN STATE_LEN SLIDE_LEN CDS_FILE (-tefl)" << endl ;
		cerr << "\t-t TRANS: Output file of TRANSITION matrix.  Default=None" << endl ;
		cerr << "\t-e EMIT : Output file of EMISSION matrix.  Default=None" << endl ;
		cerr << "\t-f FRAME: Frame to be read (1-6). Default=1" << endl ;
		cerr << "\t-l      : Output as log scores" << endl ;
		cerr << "\t-p CODE : Pseudocount (ignoring termination codons in a genetic CODE)" << endl ;
		exit ( 0 ) ;
	}

	gCommand         = argv[0] ;
	string  state_id = argv[1] ;
	int  cut_len     = atoi( argv[2] ) ;	// Summation of the current and next state sizes
	int  state_len   = atoi( argv[3] ) ;	// Length of state. If SS + SS > CL, these two states are overlapped.
	int  slide       = atoi( argv[4] ) ;

	string  trans_file = "no" ;		// Output of Transition matrix ( Default ) 
	string  emit_file  = "no" ;		// Output of Emission matrix ( Default ) 
	gFrame = 1 ;
	gScore = 0 ;		// Default: Count(Integer)
	gPseudo = -1 ; 

	int  i ;
	string  option ;
	for ( i = 5 ; i < argc ; i++ ) {
		option = argv[i] ;
		if ( option == "-t" ) {
			i++ ;
			if ( i >= argc ) { cerr  <<  "Option \"-t\" TRANSITION matrix file" << endl  ;  exit ( 0 ) ; }
			trans_file = argv[i] ;	
		} else if ( option == "-e" ) {
			i++ ;
			if ( i >= argc ) { cerr  <<  "Option \"-e\" EMISSION matrix file" << endl  ;  exit ( 0 ) ; }
			emit_file = argv[i] ;
		} else if ( option == "-f" ) {	
			// Start position (1-6): 1-> True reading frame, 4-> 1st frame of complementary strand
			i++ ;
			if ( i >= argc ) { cerr  <<  "Option \"-f\" Count start position" << endl  ;  exit ( 0 ) ; }
			gFrame = atoi ( argv[i] ) ;
			if ( gFrame < 1 || gFrame > 6 ) {
				cerr << gFrame << " Start position is incorrect." << endl ;  exit (0) ; 
			}
		} else if ( option == "-l" ) {	
			// Output Log Score
			gScore = 1 ;
		} else if ( option == "-p" ) {
			i++ ; 
			if ( i >= argc ) { 
				cerr  <<  "Option \"-p\" Codon table for pseudocount (0=for all tuples)" << endl  ;  
				exit ( 0 ) ; 
			}
			gPseudo = atoi(argv[i]) ;
		} else {
			gSeq_files.push_back ( option ) ;
		}
	}
	if ( gSeq_files.size() == 0 ) { cerr  <<  "Sequence files are not defined." << endl  ;  exit ( 0 ) ; }

	if ( state_len >= cut_len ) {
		cerr << "STATE_LEN must be smaller than CUT_LEN." << endl ;
		exit ( 0 ) ;
	}
	if ( slide > state_len ) {
		cerr << "SLIDE_LEN must be equal to or smaller than state_len." << endl ;
		exit ( 0 ) ;
	}
	if ( trans_file == "no" && emit_file == "no" ) {
		cerr << "No output files defined.  At least one file must be defined by \"-t\" or \"-e\"." << endl ;
		exit( 0 ) ; 
	}


	if ( gPseudo > -1 ) {
		cout << "Usual pseudocount" << endl ;
		cout << "Pseudocount ignoring termination codons" << endl ;
		if ( gPseudo > 0 ) {
			if ( gPseudo == 1 || gPseudo == 11 || gPseudo == 4 || gPseudo == 25 ) {
			} else {
				cerr << "Current version can treat genetic codes=1, 4, 11, and 25" << endl ;
				exit( 0 ) ;
			}
			ter_list.push_back ( "TAA" ) ;
			ter_list.push_back ( "TAG" ) ;
			if ( gPseudo == 1 || gPseudo == 11 ) {
				ter_list.push_back ( "TGA" ) ;
			}
		}
		//ter_flag = terFlag ( state_len , 1 ) ;
		ter_flag = terFlag ( state_len , gFrame ) ;
	}

	seqCounterInitialization ( cut_len , state_len , slide ) ;
	for ( i = 0 ; i < gSeq_files.size() ; i++ ) {
		gRf.readFASTA ( gSeq_files[i] ) ;
	}
	seqCount ( gRf , cut_len , state_len , slide ) ;
	if ( trans_file != "no" ) {
		if ( outTrans( trans_file, state_id , state_len ) == 0 ) {
			cerr << "Failed to make a TRANSITION matrix file." << endl ;
			return 0 ;
		}
		cout << "TRANSITION matrix : " << trans_file << endl ;
	}
	if ( emit_file != "no" ) { 
		if ( outEmit( emit_file, state_id , state_len ) == 0 ) {
			cerr << "Failed to make an EMISSION matrix file." << endl ;
			return 0 ;
		}
		cout << "EMISSION matrix : " << emit_file << endl ;
	}

	return 1 ;
}

void  seqCounterInitialization ( int  cut_len , int  state_len , int slide ) {
	int  i, k, init_code, next_code ;
	string  seq, init_seq, next_seq ;
	SEQ s ;

	int  total_init = (int) pow ( 4.0 , (double)state_len ) ;
	int  total_next = (int) pow ( 4.0 , (double)state_len ) ;
	int  slippage = cut_len - state_len ;

	gCounter1D.resize( total_init ) ;
	gCounter2D.resize( total_init ) ;
	gTotal2D.resize( total_init ) ;
	
	for ( i = 0 ; i < total_init ; i++ ) {

		gCounter2D[i].resize( total_next ) ;

		if ( gPseudo > -1 ) {

			// ============  Pseudocount Ignoring Termination Codons  ============
			// Check if termination codons are contained in tuples

			init_seq = s.decode ( i , state_len ) ;

			init_seq = init_seq.substr ( slippage , init_seq.size() - slippage ) ;
			init_code = (int) pow ( 4.0 , (double)slippage ) ;
			for ( k = 0 ; k < init_code ; k++ ) {
				next_seq = s.decode ( k , slippage ) ;
				next_seq = init_seq + next_seq ;
				next_code = s.encode ( next_seq ) ;

				if ( ter_flag[ i ] == 0  &&  ter_flag[ next_code ] == 0 ) { 	// No termination codons
					if ( slippage < 1 ) { cerr << "Splippage length must be > 0." << endl ;  exit(0) ; } 
					// Add pseudocounts
					gCounter1D[ i ]++ ;
					gCounter2D[ i ][ next_code ]++ ; 
					gTotal1D++ ;
					gTotal2D[ i ]++ ;
				}
			}	
		}
	}
}


void  seqCount ( ReadFile& rf , int  cut_len , int  state_len , int slide ) {
	int  i , k , init_code , next_code  ;
	string  seq , init_seq , next_seq;
	SEQ s ;
	
	int  slippage = cut_len - state_len ;
	for ( i = 0 ; i < rf.seq_count() ; i++ ) {
		seq = rf.getSeq( i ) ;
		if ( gFrame < 4 ) {
			seq = seq.substr ( gFrame - 1, seq.size() - gFrame + 1 ) ;
		} else {
			seq = s.reverse_complement ( seq ) ;
			seq = seq.substr ( gFrame - 4, seq.size() - gFrame + 4 ) ;
		}
//cerr << "Sequence=" << seq << endl ; 
		for ( k = 0 ; k + cut_len -1 < seq.size() ; k = k + slide ) {
			init_seq = seq.substr ( k , state_len ) ;
			next_seq = seq.substr ( k + slippage , state_len ) ;
//cerr << "Cut sequence=" << init_seq << "\t" << next_seq << endl ; 
			init_code = s.encode ( init_seq ) ;
			next_code = s.encode ( next_seq ) ;
			if ( init_code == -1  || next_code  == -1  ) { continue ; }	// Irregular Bases
			gCounter1D[ init_code ]++ ;
			gCounter2D[ init_code ][ next_code ]++ ;
			gTotal1D++ ;
			gTotal2D[ init_code ]++ ;
		}
	}
	if ( gTotal1D == 0 ) { 
		cerr << "Sequence cannot be counted." << endl ;
		exit ( 0 ) ;
	}
}

int outTrans( string&  tr , string& id , int& state_len ) {

	int i , k ; 
	SEQ s ;

	ofstream  tr_file_out ;

	tr_file_out.open ( tr.c_str() ) ;
	if ( ! tr_file_out ) {
		cerr << tr << " doesn't exist." << endl ;
		return 0 ;
	}

	tr_file_out << "### Executed program : " << gCommand << endl << endl ;
	for ( i = 0 ; i < gSeq_files.size() ; i++ ) { 
		tr_file_out << "# Transition probability from : " << gSeq_files[i]  << endl ; 
	}
	tr_file_out << "# Data type : " ;
	if ( gScore == 0 ) { 
		tr_file_out << "Observed counts" << endl ;
	} else {
		tr_file_out << "Log scores converted from observed counts" << endl ;
		tr_file_out << "# Floor score : " << FLOOR_SCORE << endl ;
	}
	if ( gPseudo > -1 ) {
		if ( gPseudo > 0 ) {
			tr_file_out << "# Pseudocounts added (termination codons are ignored)" << endl ;
		} else {
			tr_file_out << "# Pseudocounts added" << endl ;
		}
	}
	tr_file_out << "# Total sequences : " << gRf.seq_count() << endl ;
	tr_file_out << "# Total counts  : " << gTotal1D << endl ;
	tr_file_out << "# Reading frame : " << gFrame << endl ;
	tr_file_out << endl ;
	tr_file_out << "#StateNo. ID (StateNo.:Trans_prob)*" << endl ;

	tr_file_out << 0 << " BEGIN" ;

	for ( i = 0 ; i < gCounter1D.size() ; i++ ) { 
		tr_file_out << " " << i + 1 ;
		if ( gScore == 0 ) {
			tr_file_out << ":" << gCounter1D[ i ] ;
		} else {
			if ( gCounter1D[ i ] == 0 ) {
				tr_file_out << ":" << FLOOR_SCORE ;
			} else {
				tr_file_out << ":" << log ( (double) gCounter1D[ i ]  / (double) gTotal1D ) ;
			}
		}
	}
	tr_file_out << endl ;


	for ( i = 0 ; i < gCounter1D.size() ; i++ ) { 
		tr_file_out << i + 1 << " " << id ;

		if ( gTotal2D[i] == 0 ) {

		} else {
			for ( k = 0 ; k < gCounter2D[ i ].size() ; k++ ) {
				if ( gCounter2D[i][k] == 0 ) {
					//log_score = FLOOR_SCORE ;
				} else {
					if ( gScore == 0 ) { 
						tr_file_out << " " << k + 1 << ":" << gCounter2D[ i ][k]  ;
					} else {
						tr_file_out << " " << k + 1 << ":" << log( (double) gCounter2D[ i ][k]  / (double) gTotal2D[i] )  ;	
					}

				}
				
			}
		}
		tr_file_out << endl ;
	}

	tr_file_out.close() ;

	return 1 ;
}


int outEmit ( string& em , string& id , int& state_len ) {

	int i , k ; 
	SEQ s ;

	ofstream  em_file_out ;
	em_file_out.open ( em.c_str() ) ;
	if ( ! em_file_out ) {
		cerr << em << " doesn't exist." << endl ;
		return 0 ;
	}

	em_file_out << "### Executed program : " << gCommand << endl << endl ;
	for ( i = 0 ; i < gSeq_files.size() ; i++ ) { 
		em_file_out << "# Transition probability from : " << gSeq_files[i]  << endl ; 
	}
	em_file_out << "# Data type : " ;
	if ( gScore == 0 ) { 
		em_file_out << "Observed counts" << endl ;
	} else {
		em_file_out << "Log scores converted from observed counts" << endl ;
		em_file_out << "# Floor score : " << FLOOR_SCORE << endl ;
	}
	em_file_out << "# Reading frame : " << gFrame << endl ;
	em_file_out << endl ;
	em_file_out << "#StateNo. Emitted_symbol Emit_prob" << endl ;

	em_file_out << 0 << " - " ;
	if ( gScore == 0 ) { 
		em_file_out << 1 << endl ;
	} else {
		em_file_out << 0 << endl ;
	}

	for ( i = 0 ; i < gCounter1D.size() ; i++ ) { 

		em_file_out << i + 1 ;

		if ( gScore == 0 ) { 
			em_file_out << " " << s.decode ( i , state_len ) << " " << 1 << endl ;
		} else {
			em_file_out << " " << s.decode ( i , state_len ) << " " << 0 << endl ;
		}	
	}

	em_file_out.close() ;

	return 1 ;
}


vector<int>  terFlag ( int  tuple , int frame ) {
	int  i , k , j , flag , start ;
	string  seq , sub_seq ;
	SEQ s ;
	vector<int>  ter_flag ;

	int  total = (int) pow ( 4.0 , (double)tuple ) ;
	ter_flag.resize( total ) ;

	if  ( frame >= 1 && frame <= 3 ) {
		start = frame - 1 ;
	} else if ( frame >= 4 && frame <= 6 ) {
		start = frame - 3 - 1 ;
	} else {
		cerr << "Frame number is incorrect." << endl  ;
		exit ( 0 ) ;
	}
	for ( i = 0 ; i < total ; i++ ) {
		seq = s.decode ( i , tuple ) ;
		if ( frame >= 4 && frame <= 6 ) {
			seq = s.reverse ( seq ) ;
			seq = s.complement ( seq ) ;
		}
		flag = 0 ;
		for ( k = start ; k + tuple -1 < seq.size() ; k = k + 3 ) {

			sub_seq = seq.substr ( k , 3 ) ;		// Cut out a single codon one by one

			for ( j = 0 ; j < ter_list.size() ; j++ ) {
				if ( ter_list[j]  ==  sub_seq ) { flag = 1 ; }
			}
		}
		if ( flag == 1 ) { 
			ter_flag[i] = 1 ;		// Containing termination codons !!!
		} else {
			ter_flag[i] = 0 ;
		}
	}
	return ( ter_flag ) ;
}

