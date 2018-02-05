// Contents of HMM class
// Data type is "Counts"(Usually "Integer", but "Real" also OK)
// Counts are converted into Log scores in this class

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <map>
#include <algorithm>
#include <string.h>

#include "hmm.h"
#include "readfile.h"
#include "SEQ.h"

using namespace std;


HMM::HMM() {
	iFrame = 1 ;	//Default Frame
}

HMM::~HMM() {
}

// ==================================================================================
void  HMM::setHMMparameters ( string trans_file , string emit_file ) {

	vector<string>  vsd, com ;
	ReadFile  seq ;
	if ( seq.readDataLine( trans_file, vsd, com ) == 0 ) {
		cerr  <<  "Failed to obtain data from " << trans_file << endl ;
		exit ( 0 ) ;
	}

	resize_hmm ( vsd.size() ) ;

	char buffer[ BUFFER_SIZE ] ;
	char delim[] = " ,:\t\n" ;
	char  *tok ;

	int i , k , n ;

	for ( i = 0 ; i < com.size() ; i++ ) { 
		if ( com[i].substr( 0, 15 ) == "# Reading Frame" ) {
			iFrame = atoi( &com[i].at(18) ) ;
			//cerr << trans_file << " : Reading Frame=" << iFrame << endl ;
		}
	}
	for ( i = 0 ; i < iState_n ; i++ ) {

		strcpy ( buffer , vsd[i].c_str() ) ;
		tok = strtok ( buffer , delim ) ;
		if ( i != atoi ( tok ) ) {
			cerr << i << " and " << tok << " are differnet state numbers." << endl ;
			exit ( 0 ) ;
		}
		viTrans_total[i] = 0 ;

		tok = strtok ( NULL , delim ) ;
		k = 1 ; 

		while ( tok != NULL ) {
			if ( k == 1 ) { 
				vsState[i] = tok ; 
			} else if ( k % 2  ==  0 ) {
				n = atoi(tok) ;				// Possible next state numbers
				if ( n >= 0  && n <= iState_n - 1 ) {
					vviTrans_state[i].push_back( n ) ; 
				} else {
					cerr << "State=" << i << " : the next state " << n << " is incorrect." << endl ; 
					exit ( 0 ) ;
				}
			} else {
				if ( atoi(tok) < 0 ) {
					cerr << trans_file << "Seq count is incorrect." << endl ; 
					exit ( 0 ) ;
				} 
				vviTrans_count[i].push_back( atoi(tok) ) ;	
				viTrans_total[i] += atoi(tok) ;
				
			}
			tok = strtok ( NULL , delim ) ;
			k++ ; 
		}
		if ( vviTrans_state[i].size()  >  iState_n ) {
			cerr << "State=" << i << " : Number of the next states are larger than original state number.\n" ;
			exit ( 0 ) ;
		}
		if ( vviTrans_count[i].size()  !=  vviTrans_state[i].size() ) {
			cerr << "State=" << i << " : One-to-one correspondence of the next states and scores are Not satisfied.\n" ;
			exit ( 0 ) ;
		}
	}



	vsd.clear() ;
	if ( seq.readDataLine( emit_file, vsd ) == 0 ) {
		cerr  <<  "Failed to obtain data from " << emit_file << endl ;
		exit ( 0 ) ;
	}
	if ( vsd.size()  !=  iState_n ) { 
		cerr  <<  "Numbers of states are different between " << trans_file << " " << emit_file << endl ;
		exit ( 0 ) ;
	}
	
	for ( i = 0 ; i < iState_n ; i++ ) {

		strcpy ( buffer , vsd[i].c_str() ) ;

		tok = strtok ( buffer , delim ) ;
		if ( i != atoi ( tok ) ) {
			cerr << i << " and " << tok << " are differnet state numbers." << endl ;
			exit ( 0 ) ;	
		}
		viEmit_total[i] = 0 ;

		tok = strtok ( NULL , delim ) ;
		k = 1 ; 
		while ( tok != NULL ) {
			if ( k % 2 == 1 ) { 
				vvsSymbol[i].push_back( tok ) ; 
			} else {
				if ( atoi(tok) < 0 ) {
					cerr << emit_file << "Seq count is incorrect." << endl ; 
					exit ( 0 ) ;
				} 
				vviEmit_count[i].push_back( atoi(tok) ) ; 
				viEmit_total[i] += atoi(tok) ;
			}
			tok = strtok ( NULL , delim ) ;
			k++ ; 
		}
		if ( vvsSymbol[i].size ()  !=  vviEmit_count[i].size() ) {
			cerr << "Error in the line : " << vsd[i] << endl ;
			cerr << "Numbers of symbols and scores are different." << endl ;
			exit ( 0 ) ;
		}
	}

	convert_to_logscore () ;
}


void  HMM::printHMMparameters () {
	int i , k ;
	cout << "-------------------------" <<  endl ;
	for ( i = 0 ; i < vvdTrans_prob.size() ; i++ ) {
		cout << "# State=" << i << "\tState_ID=" << vsState[i] << endl ;
		cout << "# " ; 
		for ( k = 0 ; k < vvdTrans_prob[i].size() ; k++ ) {
			cout << "->" << vviTrans_state[i][k] << "(" << vvdTrans_prob[i][k] << ")\t" ;
		}
		cout << endl ;
	}
	cout << endl ; 
	for ( i = 0 ; i < vvdEmit_prob.size() ; i++ ) {
		cout << "# State=" << i << " : " ; 
		for ( k = 0 ; k < vvdEmit_prob[i].size() ; k++ ) {
			cout << vvsSymbol[i][k] << "=" << vvdEmit_prob[i][k] << "\t" ;
		}
		cout << endl ;
	}
	cout << endl ;
}


int  HMM::getNextState ( int n , int m ) {
	if ( n < 0 ) { cerr << n << " is a negative value." << endl ;   exit ( 0 ) ; }
	if ( m < 0 ) { cerr << m << " is a negative value." << endl ;   exit ( 0 ) ; }
	if ( n >= vviTrans_state.size()    ) { cerr << n << " is out of state range." << endl ;   exit ( 0 ) ; }
	if ( m >= vviTrans_state[n].size() ) { cerr << m << " is out of next state range." << endl ;   exit ( 0 ) ; }
	return ( vviTrans_state[n][m] ) ;
}


string  HMM::getSymbol ( int n , int m ) {
	if ( n < 0 ) { cerr << n << " is a negative value." << endl ;   exit ( 0 ) ; }
	if ( m < 0 ) { cerr << m << " is a negative value." << endl ;   exit ( 0 ) ; }
	if ( n >= vvsSymbol.size()    ) { cerr << n << " is out of state range." << endl ;   exit ( 0 ) ; }
	if ( m >= vvsSymbol[n].size() ) { cerr << m << " is out of emit range." << endl ;   exit ( 0 ) ; }
	return ( vvsSymbol[n][m] ) ;
}

int  HMM::getState () { return ( iState_n ) ; }
vector<vector<int> > HMM::getTransState() { return ( vviTrans_state ) ; }
vector<vector<double> > HMM::getTransProbVector() { return ( vvdTrans_prob ) ; }
vector<vector<double> > HMM::getEmitProbVector() { return ( vvdEmit_prob ) ; }

vector<vector<int> > HMM::getTransCountVector() { return ( vviTrans_count ) ; }
vector<int> HMM::getTransTotal() { return ( viTrans_total ) ; }
vector<vector<int> > HMM::getEmitCountVector() { return ( vviEmit_count ) ; }
vector<int> HMM::getEmitTotal() { return ( viEmit_total ) ; }

vector<string> HMM::getStateVector() { return ( vsState ) ; }
vector<vector<string> > HMM::getSymbolVector() { return ( vvsSymbol ) ; }
int HMM::getFrame() { return ( iFrame ) ; }


void  HMM::resize_hmm ( int  n ) {
	iState_n = n ;

	if ( vviTrans_count.size() > 0 ) { 

		viTrans_total.clear() ;
		vviTrans_count.clear() ;
		vvdTrans_prob.clear() ;
		vviTrans_state.clear() ;
		vsState.clear() ;

		viEmit_total.clear() ;
		vviEmit_count.clear() ;
		vvdEmit_prob.clear() ;
		vvsSymbol.clear() ;

	}

	viTrans_total.resize( n ) ;
	vviTrans_count.resize ( n ) ;
	vvdTrans_prob.resize ( n ) ;

	vviTrans_state.resize( n ) ;
	vsState.resize ( n ) ;

	viEmit_total.resize( n ) ; 
	vviEmit_count.resize ( n ) ;
	vvdEmit_prob.resize ( n ) ;

	vvsSymbol.resize ( n ) ;
}


void HMM::reserve_hmm ( int n ) {

	if ( vviTrans_count.size() == 0 ) { 
		resize_hmm(n) ;
	}

	viTrans_total.reserve( n ) ;
	viEmit_total.reserve( n ) ;
	vsState.reserve( n ) ;

	vviTrans_state.reserve( n ) ;

	vviTrans_count.reserve( n ) ;
	vvdTrans_prob.reserve( n ) ;
	vviEmit_count.reserve( n ) ;
	vvdEmit_prob.reserve( n ) ;

	vvsSymbol.reserve( n ) ;

	int i ; 
	for ( i = 0 ; i < n ; i++ ) { 
		vviTrans_state[i].reserve(n) ; 
		vviTrans_count[i].reserve(n) ; 
		vvdTrans_prob[i].reserve(n) ; 
		vviEmit_count[i].reserve(n) ; 
		vvdEmit_prob[i].reserve(n) ; 
		vvsSymbol[i].reserve(n) ; 
	}
}



void  HMM::weightNormalization ( vector<vector<double> >& trans_weight ) {
	int  i , k ;
	vector<double>  trans_total ;
	trans_total.resize( trans_weight.size() ) ;

	for ( i = 0 ; i < trans_weight.size() ; i++ ) {
		if ( trans_weight[i].size() != trans_weight.size()  ) {
			cerr  <<  "Weight matrix ( TRANSITION ) is incorrect. ( weightNormalization )" << endl ;
			exit ( 0 ) ;
		}
		trans_total[i] = 0.0 ;
		for ( k = 0 ; k < trans_weight[i].size() ; k++ ) {
			trans_total[i] += trans_weight[i][k] ;
		}
	}
	
	for ( i = 0 ; i < trans_weight.size() ; i++ ) {
		for ( k = 0 ; k < trans_weight[i].size() ; k++ ) {
			if ( trans_weight[i][k] == 0 ) {
			 	trans_weight[i][k] = FLOOR_SCORE ;
			} else {
				trans_weight[i][k] = log ( trans_weight[i][k] / trans_total[i] ) ;
			}
		}	
	}
}

void  HMM::weightNormalization ( vector<double>& init_weight, vector<vector<double> >& trans_weight ) {
	int  i ;
	double init_total ;
	if ( init_weight.size() == 0 ) { 
		cerr  <<  "Weighted vector ( INITIATION ) is incorrect. ( weightNormalization )" << endl ;
		exit ( 0 ) ;
	}
	for ( i = 0 ; i < init_weight.size() ; i++ ) {
		init_total += init_weight[i] ;
	}
	if ( init_total < 1e-8 ) { 
		cerr  <<  "Weighted vector ( INITIATION ) is incorrect. init_total is much small. ( weightNormalization )" << endl ;
		exit ( 0 ) ;
	}
	
	for ( i = 0 ; i < init_weight.size() ; i++ ) {
		if ( init_weight[i] == 0 ) {
			init_weight[i] = FLOOR_SCORE ;
		} else {
			init_weight[i] = log ( init_weight[i] / init_total ) ;
		}	
	}
	weightNormalization ( trans_weight ) ;
}


void  HMM::printWeight ( vector<vector<double> >& trans_weight ) {
	int i, k ;
	i = trans_weight.size() - 1 ;
	cout << "# HMM No.\tTrans_weights (0-" << i << ")" << endl ;
	cout << "# Probability Weight : " << endl ;
	for ( i = 0 ; i < trans_weight.size() ; i++ ) {
		cout << "# " << i ;
		for ( k = 0 ; k < trans_weight[i].size() ; k++ ) {
			cout << "\t" << exp(trans_weight[i][k]) ;
		}	
		cout << endl ; 
	}
	cout << "# Converted to Log Scores : " << endl ;
	for ( i = 0 ; i < trans_weight.size() ; i++ ) {
		cout << "# " << i ;
		for ( k = 0 ; k < trans_weight[i].size() ; k++ ) {
			cout << "\t" << trans_weight[i][k] ;
		}	
		cout << endl ; 
	}
}
void  HMM::printWeight ( vector<double>& init_weight, vector<vector<double> >& trans_weight ) {
	int i ;
	i = init_weight.size() - 1 ;
	cout << "# HMM No.\tInit_weights (0-" << i << ")" << endl ;
	cout << "# Probability Weight : " << endl ;
	for ( i = 0 ; i < init_weight.size() ; i++ ) {
		cout << "# " << i << "\t" << exp(init_weight[i]) << endl ; 
	}
	cout << "# Converted to Log Scores : " << endl ;
	for ( i = 0 ; i < init_weight.size() ; i++ ) {
		cout << "# " << i << "\t" << init_weight[i] << endl ; 
	}
	cout << endl ;
	printWeight ( trans_weight ) ; 
}

void HMM::connectStates ( vector<HMM>& hmm_vec, vector<vector<double> >& trans_weight ) {
	if ( hmm_vec.size() == 0 ) { 
		cerr << "HMM objects are not defined. ( connectStates )" << endl ;
		exit ( 0 ) ;
	}

	int i , k , j , l , m, start ;
	int  n = 0 ;
	vector<int> corStateNums ; 

	// ==========  Check of Matirx Size  and  Set State Size  ==========
	if ( trans_weight.size() != hmm_vec.size() ) { 
		cerr << "Weight matirx is incorrect. ( connectStates )" << endl ;
		exit ( 0 ) ;
	}

	corStateNums.resize( hmm_vec.size() ) ; 
	corStateNums[0] = 0 ; 
	for ( i = 0 ; i < hmm_vec.size() ; i++ ) {

		if ( trans_weight[i].size() != hmm_vec.size() ) { 
			cerr << "Weight matirx is incorrect. ( connectStates )" << endl ;
			exit ( 0 ) ;
		}

		// Change of STATE Numbers
		if ( i > 0 ) {
			corStateNums[i] = n ; 
		}
		n += hmm_vec[i].getState() - 1 ;	// Except for "BEGIN"
	}

	n++ ;	// Add "BEGIN"
	resize_hmm( n ) ;		// "n" is total size of states

	n = 0 ;
	for ( i = 0 ; i < hmm_vec.size() ; i++ ) {
			
		if ( i == 0 ) { start = 0 ; } else { start = 1 ; }

		//for ( k = 0 ; k < hmm_vec[i].iState_n ; k++ ) {
		for ( k = start ; k < hmm_vec[i].iState_n ; k++ ) {

			vsState[n] = hmm_vec[i].vsState[k] ;	// Inheritance of STATE ID

			for ( j = 0 ; j < hmm_vec.size() ; j++ ) {
				for ( l = 0 ; l < hmm_vec[j].vvdTrans_prob[k].size() ; l++ ) {
					vvdTrans_prob[n].push_back ( hmm_vec[j].vvdTrans_prob[k][l] + trans_weight[i][j] ) ;
					//hmm.vviTrans_state[n].push_back ( hmm_vec[j].vviTrans_state[k][l] ) ;
					vviTrans_state[n].push_back ( hmm_vec[j].vviTrans_state[k][l] +  corStateNums[j] ) ;
				}
			}

			for ( m = 0 ; m < hmm_vec[i].vvsSymbol[k].size() ; m++ ) {
				vvdEmit_prob[n].push_back ( hmm_vec[i].vvdEmit_prob[k][m] ) ;
				vvsSymbol[n].push_back ( hmm_vec[i].vvsSymbol[k][m] ) ;
			}
			n++ ;
		}
	}
}

void HMM::connectStates ( vector<HMM>& hmm_vec, vector<double>& init_weight, vector<vector<double> >& trans_weight ) {
	if ( hmm_vec.size() == 0 ) { 
		cerr << "HMM objects are not defined. ( connectStates )" << endl ;
		exit ( 0 ) ;
	}

	int i , k , j , l , m, start ;
	int  n = 0 ;
	vector<int> corStateNums ; 

	// ==========  Check of Matirx Size  and  Set State Size  ==========
	if ( init_weight.size() != hmm_vec.size() ) { 
		cerr << "Initiation Weight vector is incorrect. ( connectStates )" << endl ;
		exit ( 0 ) ;
	}

	if ( trans_weight.size() != hmm_vec.size() ) { 
		cerr << "Transition Weight matirx is incorrect. ( connectStates )" << endl ;
		exit ( 0 ) ;
	}

	corStateNums.resize( hmm_vec.size() ) ; 
	corStateNums[0] = 0 ;
	for ( i = 0 ; i < hmm_vec.size() ; i++ ) {

		if ( trans_weight[i].size() != hmm_vec.size() ) { 
			cerr << "Weight matirx is incorrect. ( connectStates )" << endl ;
			exit ( 0 ) ;
		}

		// Change of STATE Numbers
		if ( i > 0 ) {
			corStateNums[i] = n ; 
		}
		n += hmm_vec[i].getState() - 1 ;	// Except for "BEGIN"
	}

	n++ ;	// Add "BEGIN"
	resize_hmm( n ) ;		// "n" is total size of states

	n = 0 ;
	for ( i = 0 ; i < hmm_vec.size() ; i++ ) {
			
		if ( i == 0 ) { start = 0 ; } else { start = 1 ; }

		for ( k = start ; k < hmm_vec[i].getState() ; k++ ) {

			vsState[n] = hmm_vec[i].vsState[k] ;	// Inheritance of STATE ID
			for ( j = 0 ; j < hmm_vec.size() ; j++ ) {
				for ( l = 0 ; l < hmm_vec[j].vvdTrans_prob[k].size() ; l++ ) {
					if ( n == 0 ) { 
						vvdTrans_prob[n].push_back ( hmm_vec[j].vvdTrans_prob[k][l] + init_weight[j] ) ;
					} else { 
						vvdTrans_prob[n].push_back ( hmm_vec[j].vvdTrans_prob[k][l] + trans_weight[i][j] ) ;
					}
					vviTrans_state[n].push_back ( hmm_vec[j].vviTrans_state[k][l] +  corStateNums[j] ) ;
				}
			}

			for ( m = 0 ; m < hmm_vec[i].vvsSymbol[k].size() ; m++ ) {
				vvdEmit_prob[n].push_back ( hmm_vec[i].vvdEmit_prob[k][m] ) ;
				vvsSymbol[n].push_back ( hmm_vec[i].vvsSymbol[k][m] ) ;
			}
			n++ ;
		}
	}
}


void HMM::outputTransMatrix ( string& tr ) {
	int i , k ;
	ofstream  tr_file_out ;
        
	tr_file_out.open ( tr.c_str() ) ;
	if ( ! tr_file_out ) {
		cerr << tr << " doesn't exist." << endl ;
		exit(0) ;
	}
        
	tr_file_out << "# This matrix is made from HMM::outputTransMatrix." << endl ; 
	tr_file_out << "# Therefore, the data are shown as log scores." << endl ; 
	tr_file_out << endl ; 
	tr_file_out << "#StateNo. ID (StateNo. Trans_Prob)*" << endl ;

	for ( i = 0 ; i < vvdTrans_prob.size() ; i++ ) {

		tr_file_out << i << " " << vsState[i] ;
		for ( k = 0 ; k < vvdTrans_prob[i].size() ; k++ ) {
			tr_file_out << " " << vviTrans_state[i][k] << " " << vvdTrans_prob[i][k] ;
		}
		tr_file_out << endl ;
	}
}

void HMM::outputEmitMatrix ( string& em ) {
	int i , k ;
	ofstream  em_file_out ;
	em_file_out.open ( em.c_str() ) ;
	if ( ! em_file_out ) {
		cerr << em << " doesn't exist." << endl ;
		exit(0) ;
	}
	em_file_out << "# This matrix is made from HMM::outputEmitMatrix." << endl ; 
	em_file_out << "# Therefore, the data are shown as log scores." << endl ; 
	em_file_out << endl ; 
	em_file_out << "#StateNo. Emitted_Symbol Emit_Prob" << endl ;
	for ( i = 0 ; i < vvdEmit_prob.size() ; i++ ) {
		em_file_out << i ;
		for ( k = 0 ; k < vvdEmit_prob[i].size() ; k++ ) {
			em_file_out << " " << vvsSymbol[i][k] << " " << vvdEmit_prob[i][k] ;
		}
		em_file_out << endl ;
	}
	em_file_out << endl ;
}



void HMM::convert_to_logscore () {
	int  i, k; 

	for ( i = 0 ; i < iState_n ; i++ ) {
		if ( viEmit_total[i] == 0 ) {
			cerr << "All of the emit probs. at the state=" << i << " are ZERO." << endl ;
			exit(0) ;
		}
		vvdEmit_prob[i].resize( vviEmit_count[i].size() ) ; 
		vvdTrans_prob[i].resize( vviTrans_count[i].size() ) ; 
		for ( k = 0 ; k < vviEmit_count[i].size() ; k++ ) {
			if ( vviEmit_count[i][k] < 1 ) { 
				vvdEmit_prob[i][k] = log ( 0.5 / viEmit_total[i] ) ;
			} else {
				vvdEmit_prob[i][k] = log ( (double) vviEmit_count[i][k] / viEmit_total[i] ) ;
			}
		}
		for ( k = 0 ; k < vviTrans_count[i].size() ; k++ ) {
			if ( vviTrans_count[i][k] < 1 ) { 
				vvdTrans_prob[i][k] = log ( 0.5 / viTrans_total[i] ) ;
			} else {
				vvdTrans_prob[i][k] = log ( (double) vviTrans_count[i][k] / viTrans_total[i] ) ;
			}
		}
	}
}

void HMM::convert_to_probability () {
	int  i, k; 

	for ( i = 0 ; i < iState_n ; i++ ) {
		if ( viEmit_total[i] == 0 ) {
			cerr << "All of the emit probs. at the state=" << i << " are ZERO." << endl ;
			exit(0) ;
		}
		vvdEmit_prob[i].resize( vviEmit_count[i].size() ) ; 
		vvdTrans_prob[i].resize( vviTrans_count[i].size() ) ; 
		for ( k = 0 ; k < vviEmit_count[i].size() ; k++ ) {
			vvdEmit_prob[i][k] = (double) vviEmit_count[i][k] / viEmit_total[i] ;
		}
		for ( k = 0 ; k < vviTrans_count[i].size() ; k++ ) {
			vvdTrans_prob[i][k] = (double) vviTrans_count[i][k] / viTrans_total[i] ;
		}
	}
}


// fr : reading frame (1-6)
void HMM::subtractCounts ( string& qseq, int fr, int cut_len, int state_len ) {
	int slide = cut_len - state_len ;
	int  i , k , init_code , next_code ;
	string  seq, init_seq , next_seq ;
	SEQ s ;
	int j ; 
	int  slippage = cut_len - state_len ;
	if ( fr < 4 ) {
		seq = qseq.substr ( fr - 1, qseq.size() - fr + 1 ) ;
	} else {
		seq = s.reverse_complement ( qseq ) ;
		seq = seq.substr ( fr - 4, seq.size() - fr + 4 ) ;
	}
	for ( k = 0 ; k + cut_len -1 < seq.size() ; k = k + slide ) {
		init_seq = seq.substr ( k , state_len ) ;
		next_seq = seq.substr ( k + slippage , state_len ) ;
		init_code = s.encode ( init_seq ) ;
		next_code = s.encode ( next_seq ) ;
		if ( init_code == -1  || next_code  == -1  ) { continue ; }	// Irregular Bases

		if ( viTrans_total[init_code + 1] == 1 ) { 
			cerr << "Trans total will be ZERO. Maybe Pseudocounts are necessary." << endl ;
			exit(0) ;
		}
//
// If init_code + 1 is the state num.
// Usual HMM training will meet this condition.
//
		for ( j = 0 ; j < vviTrans_state[init_code + 1].size() ; j++ ) { 
			if ( vviTrans_state[init_code + 1][j] == next_code + 1 ) { break ; } 
		}
		vviTrans_count[0][init_code + 1]-- ;
		vviTrans_count[init_code + 1][j]-- ;
		viTrans_total[init_code + 1]-- ;
	}
	convert_to_logscore () ; 
}

void HMM::copy(HMM& hmm) {
	//cerr << hmm.getState() << "=" << iState_n << endl ; 

	iFrame = hmm.getFrame() ; 

	vviTrans_count = hmm.getTransCountVector() ;
	vvdTrans_prob = hmm.getTransProbVector() ;
	viTrans_total = hmm.getTransTotal() ; 
	vviEmit_count = hmm.getEmitCountVector() ;
	vvdEmit_prob = hmm.getEmitProbVector() ;
	viEmit_total = hmm.getEmitTotal() ; 
	
	vviTrans_state = hmm.getTransState() ;
	vsState = hmm.getStateVector() ;
	vvsSymbol = hmm.getSymbolVector() ;

}

