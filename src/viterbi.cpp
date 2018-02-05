// Contents of HMM Class
// If you want skip impossible transition, compile by "gcc viterbi -D SKIP"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <map>
#include <algorithm>

#include "hmm.h"
#include "viterbi.h"

using namespace std;

Viterbi::Viterbi() {
	//iData_type = 0 ; 
}

Viterbi::~Viterbi() {
}

// ==================================================================================
void  Viterbi::calcViterbi ( vector<string>& sym_vec , int first , int last ) {

	vector<vector<int> > pre_ptr , ptr ;		// State k , Sequence of States 
	vector<double> pre_path_prob , path_prob ;		// State l
	vector<int>  pre_possible_state , possible_state ;
	int  i, k, j, l , n, nn, flag, pre_state ;

	// =========  Check of path range ( first - last )  ==========
	if ( first >= last ) { 	
		cerr << "viPath range is incorrect." << endl ;  
		exit ( 0 ) ; 
	}
	if ( first < 0 ) { 
		cerr << "The first position in the path : " << first << " < 0" << endl ;   
		exit ( 0 ) ;  
	}
	if ( last > sym_vec.size() - 1 ) { 
		cerr << "The last position in the path : " << last << " > " << sym_vec.size() - 1 << endl ;   
		exit ( 0 ) ;  
	}

	// =========  Initialization ==========
	viPath.resize( 0 ) ;
	iFirst_path = first ;
    iLast_path  = last ;
	iPath = 0 ; 
	vsStopped.resize(0) ; 
	dProb = 0.0 ;

	vdInit_prob.resize( iState_n ) ;
	pre_ptr.resize( iState_n ) ;
	ptr.resize( iState_n ) ;
	
	double max = 0.0 ;
	char pos[4] ;
	string spos ; 
	
	// ==========  Computation in Log Space  ==========
	// First Step is the Transition from "Begin" -> Initial State
	// The Last Step ( Any State -> "End" ) is Ignored Here ( because the transition probabilities are identical )

#ifdef  SKIP
i = first ;
while ( i <= last ) {
#endif

	for ( nn = 0 ; nn < vvdTrans_prob[0].size() ; nn++ ) {

		//if ( vvdTrans_prob[0][nn] <= FLOOR_SCORE ) { continue ; } 	// Transition is impossible !!! 

		l = vviTrans_state[0][nn] ;		// Transition from BEGIN
		for ( j = 0 ; j < vvsSymbol[l].size() ; j++ ) {
			if ( vvsSymbol[l][j]  ==  sym_vec[ first ] ) {		// Transition is possible !!! 
#ifdef  SKIP
				if ( pre_possible_state.size() > 0 ) {
					path_prob.push_back( max + vvdTrans_prob[0][nn] + vvdEmit_prob[l][j] ) ;
                    ptr[l] = ptr[ pre_state ] ;
				} else {
#endif
					path_prob.push_back( vvdTrans_prob[0][nn] + vvdEmit_prob[l][j] ) ;
					vdInit_prob[l] = vvdTrans_prob[0][nn] + vvdEmit_prob[l][j] ;
#ifdef  SKIP
				}
#endif
				ptr[l].push_back(l) ;
				possible_state.push_back( l ) ;
	
			}
		}
	}


	if ( path_prob.size() == 0 ) {
		sprintf (pos, "%d",  first + 1 ) ;
		spos = pos ; 
		vsStopped.push_back( spos + ":" + sym_vec[first] ) ; 

#ifdef  SKIP
		cerr << "Transition is impossible at the state " ;
		cerr << first << " : " << sym_vec[first] << " (Skipped)" << endl ;
		
		if ( pre_possible_state.size() > 0 ) {
			for ( n = 0 ; n < pre_possible_state.size() ; n++ ) {
				k = pre_possible_state[n] ;
				ptr[k].push_back(iState_n) ;
			}
		} else {
			for ( l = 1 ; l < iState_n ; l++ ) {
				ptr[l].push_back(iState_n) ;
			}
		}
		first++ ;
		i = first ;
		continue ;
#endif

		cerr << "Transition is impossible." << endl ;
		return  ;
	}
	iPath++ ; 

	for ( i = first + 1 ; i <= last ; i++ ) {
		pre_path_prob = path_prob ;
		pre_possible_state = possible_state ;
		pre_ptr = ptr ;

		path_prob.resize( 0 ) ;
		possible_state.resize( 0 ) ;

		//for ( l = 0 ; l < iState_n ; l++ ) {					// Current State
		for ( l = 1 ; l < iState_n ; l++ ) {
			for ( j = 0 ; j < vvsSymbol[l].size() ; j++ ) {			// Emission
				if ( vvsSymbol[l][j]  ==  sym_vec[i] ) {			// Emission is possible !!!
					flag = 0 ;
					for ( n = 0 ; n < pre_possible_state.size() ; n++ ) {		

						// k = Previous State Number
						k = pre_possible_state[n] ;
						
						// Transition check to Current state
						for ( nn = 0 ; nn < vvdTrans_prob[k].size() ; nn++ ) {

							// Transition is possible !!!  -> flag = 1
							if ( vviTrans_state[k][nn] == l ) {
							//if ( vviTrans_state[k][nn] == l && vvdTrans_prob[k][nn] > FLOOR_SCORE ) {

								if ( flag == 0 ) {
								    max = pre_path_prob[n] + vvdTrans_prob[k][nn] ;
									//cerr << max << "\t" << pre_path_prob[n] << "\t" << vvdTrans_prob[k][nn] << endl ;
								    pre_state= k ; 
								    flag = 1 ;
								} else {
								    if ( pre_path_prob[n] + vvdTrans_prob[k][nn] > max ) {
										max = pre_path_prob[n] + vvdTrans_prob[k][nn] ;
										pre_state = k ; 
								    }
								}
							}
						}
					}
					if ( flag == 1 ) {
						path_prob.push_back( max + vvdEmit_prob[l][j] ) ;
						possible_state.push_back( l ) ;
						ptr[l] = pre_ptr[ pre_state ] ;
						ptr[l].push_back(l) ;
					}
				}
			}
		}
		if ( path_prob.size() == 0 ) {
			sprintf (pos, "%d",  i + 1 ) ;
			spos = pos ; 
			vsStopped.push_back( spos + ":" + sym_vec[i] ) ; 
			if ( i == sym_vec.size() - 1 ) {
				// Perhaps termination codon 
				cerr << "Transition is stopped at the final symbol." << endl ;
				i++ ; 
				break ; 
			}
			cerr << "Transition is stopped at the state " << i << " : " << sym_vec[i]  ;

#ifdef  SKIP
			cerr << " (Skipped)" ;
			flag = 0 ;
			for ( n = 0 ; n < pre_possible_state.size() ; n++ ) {
				k = pre_possible_state[n] ;
				if ( flag == 0 ) {
					max = pre_path_prob[n] ;
					pre_state = k ;
					flag = 1 ;
				} else {
					if ( pre_path_prob[n] > max ) {
						max = pre_path_prob[n] ;
						pre_state = k ;
					}
				}
			}
			ptr[pre_state].push_back(iState_n) ;

			i++ ;
			first = i ;
#endif

			cerr << endl ;
			break ; 
		}
		iPath++ ; 
	}
#ifdef  SKIP
}
#endif

	max = path_prob[0] ;
	nn = 0 ;
	for ( n = 1 ; n < path_prob.size() ; n++ ) {
		if ( path_prob[n] > max ) { 
			max = path_prob[n] ;
			nn = n ;
		}
	}

	l = possible_state[nn] ;			// Actual End State
	viPath = ptr[l] ;
	dProb  = max ;

}


// Uncap option
void  Viterbi::uncap () {
	//dProb -= vvdTrans_prob[0][viPath[0]] ;
//int i ;
//cerr << viPath[0] ;
//for ( i = 1 ; i < viPath.size() ; i++ ) { 
//	cerr << "," << viPath[i] ;
//}
//cerr << endl ;
	dProb -= vdInit_prob[viPath[0]] ;
	iPath-- ;
}


/*
void  Viterbi::statistics () {
	vector<string>  nr_state_id ;
	vector<int>  nr_state_count ;
	int i , pos ;
	string str = "@" ;
	for ( i = 0 ; i < viPath.size() ; i++ ) { str += hmm.vsState[ viPath[i] ] ; }
	//nr_state_id.push_back( hmm.vsState[ viPath[0] ] ) ;
	//nr_state_count.push_back( 1 ) ;
cout << "#" << str << "#" << endl ;
	for ( i = 0 ; i < viPath.size() ; i++ ) {
		pos = str.find( hmm.vsState[ viPath[i] ] ) ;
cout << "@" << pos <<  ;
		if ( pos > 0 ) {
			continue ;
		} else {
			nr_state_id.push_back( hmm.vsState[ viPath[i] ] ) ;
			nr_state_count.push_back( str.count( hmm.vsState[ viPath[i] ] ) ; ) ;
		}
	}
	pos = 0 ;
cout << endl ;
cout << nr_state_id.size() << endl ;

	for ( i = 0 ; i < nr_state_id.size() ; i++ ) {
		cout << i << " " << nr_state_id[i] << " " << nr_state_count[i] << endl ;
		pos += nr_state_count[i] ;
	}
	cout << "Total States=" << pos <<endl ;


}
*/


std::string  Viterbi::getStateSequence () {
	int i ;
	std::stringstream ss ;
#ifdef  SKIP
	if ( viPath[0] == iState_n ) {
		ss << "-" ;
	} else {
#endif
		ss << vsState[ viPath[0] ] ;
#ifdef  SKIP
	}
#endif

	for ( i = 1 ; i < viPath.size() ; i++ ) {
#ifdef  SKIP
		if ( viPath[i] == iState_n ) {
			ss << ",-" ;
			continue ; 
		}
#endif
		ss << "," << vsState[ viPath[i] ] ;
	}
	return ( ss.str() ) ;
}

std::string  Viterbi::getAllProbs ( vector<string>& sym_vec ) {
	int i, k, l, m, nn ; 
	int flag ; 

	std::stringstream ss ;
	ss.precision(7) ; 
	k = viPath[0] ;
#ifdef  SKIP
	if ( k == iState_n ) {
		ss << "-" ;
	} else { 
#endif
		ss << vdInit_prob[k] ;
#ifdef  SKIP
	}
#endif

	for ( i = 1 ; i < viPath.size() ; i++ ) {
#ifdef  SKIP
		if ( viPath[i] == iState_n ) {
			ss << ",-" ;
			continue ; 
		}
#endif
		flag = 0 ; 
		m = viPath[i] ;
		for ( nn = 0 ; nn < vviTrans_state[k].size() ; nn++ ) { 	
			if ( vviTrans_state[k][nn] != m ) { continue ; } 
			for ( l = 0 ; l < vvsSymbol[m].size() ; l++ ) { 
				if ( vvsSymbol[m][l] != sym_vec[i] ) { continue ; } 
				flag = 1 ; 
				break ;
			}
			break ; 
		}
		if ( flag == 0 ) {
			cerr << "Pos=" << i+1 << " : " << sym_vec[i] << " : " << "Emission is impossible." << endl ;
			return ( "" ) ;
		}
		ss << "," << vvdTrans_prob[k][nn] + vvdEmit_prob[m][l] ;
		k = m ;
	}
	return ( ss.str() ) ;
}


// Fisrt_pos  Last_pos  Stopped_or_skipped  Transit_path_length  Prob_log_score  HT_index  (Viterbi path)
std::string  Viterbi::getSummary ( int flag ) {
	int i ;
	std::stringstream ss ;
	if ( iPath < 1 ) { 
		cerr << "Path length=ZERO : Maybe errors has occurred." << endl ;
		return ( ss.str() ) ;
	}
	ss << iFirst_path + 1 << "\t" << iLast_path + 1 << "\t" ;
	if ( vsStopped.size() > 0 ) { 
		ss << vsStopped[0] ;
		for ( i = 1 ; i < vsStopped.size() ; i++ ) { ss << "," << vsStopped[i] ; }
	} else { 
		ss << "None" ;
	}
	ss.precision(7) ; 
	ss << "\t" << iPath << "\t" << dProb << "\t" << dProb / iPath ;
	if ( flag == 1 ) {
		ss << "\t" << vsState[ viPath[0] ] ;
		for ( i = 1 ; i < viPath.size() ; i++ ) { ss << "," << vsState[ viPath[i] ] ; }
	}		
	return ( ss.str() ) ;
}

int  Viterbi::getPathLen () { return ( iPath ) ; } 


