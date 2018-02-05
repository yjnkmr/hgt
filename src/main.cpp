/* Copyright: Yoji Nakamura, NRIFS, FRA */
/* HMM HGT Detection Program */

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <map>
#include <dirent.h>
#include <algorithm>

#include "hmm.h"
#include "viterbi.h"
#include "readfile.h"
#include "SEQ.h"


using namespace std;

/* Argument : DNA_SEQ_FILE TUPLE_in_STATE SLIDE (TRANSITION_FILE EMISSION_FILE)* -ITtemMq WEIGHT_PARAMETERS */
int main ( int argc , char *argv[] ) {

	if ( argc < 6 ) { 
		cerr << "Usage: ./PROGRAM CDS_FILE CUT_LEN STATE_LEN (-teITmM) (TRANS EMIT)" << endl ;
		cerr << "\t-t TRANS EMIT_FILES: One TRANSITION file and multiple EMISSION files" << endl ; 
		cerr << "\t-e EMIT TRANS_FILES: One EMISSION file and multiple TRANSITION files" << endl ;
		cerr << "\tTRANS EMIT         : Default setting of TRANSITION and EMISSION files" << endl ;
//		cerr << "\t-q : Quiet output (=No header)" << endl ;
		cerr << "\t-s : Subtract nucleotide frequencies of query sequence itself" << endl ;
		cerr << "\t-u : Uncap computation" << endl ; 
		cerr << "\t-S : Output state sequence" << endl ;
		cerr << "\t-P : Output all probabilities" << endl ;
		exit ( 0 ) ;
	}

	vector<string>  trans_files ;
	vector<string>  emit_files ;
	int i , k, frame ;

	string seq_file = argv[1] ;
	int cut_len    = atoi( argv[2] ) ;
	int state_len  = atoi( argv[3] ) ;
	if ( cut_len - state_len < 1 ) { 
		cerr << "CUT_LEN must be larget than STATE_LEN" << endl ; 
		exit(0) ;
	}

	string  arg1, arg2 ;
	int init_w_ratio  = 5 ;
    int trans_w_ratio = 5 ;
	int subtract_flag = 0 ; 
	int show_stateseq = 0 ;
	int show_allprobs = 0 ;
	int header_off = 1 ;
	int uncap = 0 ; 
	
	vector<int> pre_hmm_nums, hmm_nums, hmm_weights, weight_flags ;

	for ( i = 4 ; i < argc ; i++ ) {
		arg1 = argv[i] ;
		if ( arg1 == "-s" ) { 
			subtract_flag = 1 ; 
		} else if ( arg1 == "-q" ) {
			header_off = 1 ; 
		} else if ( arg1 == "-u" ) {
			uncap = 1 ;
		} else if ( arg1 == "-S" ) { 
			show_stateseq = 1 ;
		} else if ( arg1 == "-P" ) {
			show_allprobs = 1 ;
		} else if ( arg1 == "-I" ) {
			if ( ++i >= argc ) { cerr  <<  "arg1 \"-i\" INITIALTION WEIGHT" << endl  ;  exit ( 0 ) ; }
			init_w_ratio = atoi( argv[i] ) ;	
		} else if ( arg1 == "-T" ) {
			if ( ++i >= argc ) { cerr  <<  "arg1 \"-t\" TRANSITION WEIGHT" << endl  ;  exit ( 0 ) ; }
			trans_w_ratio = atoi( argv[i] ) ;
		} else if ( arg1 == "-m" ) {
			if ( i + 3 >= argc ) { cerr  <<  "arg1 \"-m\" X Y TRANSITION WEIGHT" << endl  ;  exit ( 0 ) ; }
			pre_hmm_nums.push_back ( atoi( argv[++i] ) ) ;
			hmm_nums.push_back ( atoi( argv[++i] ) ) ;	
			hmm_weights.push_back ( atoi( argv[++i] ) ) ;
			weight_flags.push_back ( 0 ) ;
		} else if ( arg1 == "-M" ) {			// Matually weighed
			if ( i + 3 >= argc ) { cerr  <<  "arg1 \"-M\" X Y TRANSITION WEIGHT" << endl  ;  exit ( 0 ) ; }
			pre_hmm_nums.push_back ( atoi( argv[++i] ) ) ;
			hmm_nums.push_back ( atoi( argv[++i] ) ) ;	
			hmm_weights.push_back ( atoi( argv[++i] ) ) ;
			weight_flags.push_back ( 1 ) ;
		} else if ( arg1 == "-t" ) {
			// One TRANS file is SHARED among models, each of which has a different EMIT file.
			// -t TRANS  EMIT_1 EMIT_2 ... EMIT_N
			if ( i + 1 >= argc ) { cerr << "TRANSITION File of " << arg1 << " is not defined." << endl ;  exit ( 0 ) ; }
			arg1 = argv[++i] ;
			while ( i + 1 < argc ) {
				arg2 = argv[++i] ; 
				if ( arg2.at(0) == '-' ) {
					i-- ;
					break ; 
				} 
				trans_files.push_back ( arg1 ) ;
				emit_files.push_back ( arg2 ) ;
			}
		} else if ( arg1 == "-e" ) {
			// One EMIT file is SHARED among models, each of which has a different TRANS file.
			// -e EMIT  TRANS_1 TRANS_2 ... TRANS_N
			if ( i + 1 >= argc ) { cerr << "EMISSION File of " << arg1 << " is not defined." << endl ;  exit ( 0 ) ; }
			arg1 = argv[++i] ;
			while ( i + 1 < argc ) {
				arg2 = argv[++i] ; 
				if ( arg2.at(0) == '-' ) {
					i-- ;
					break ; 
				} 
				trans_files.push_back ( arg2 ) ;
				emit_files.push_back ( arg1 ) ;
			}
		} else if ( arg1 == "-d" && i + 1 < argc ) {
			arg1 = argv[++i] ;
			DIR *dir = opendir ( arg1.c_str() ) ;
			struct dirent *dent ;
			if( dir ){
				while( dent = readdir( dir ) ) {
					arg2 = dent->d_name ;
					//cout << dent->d_name << endl;
					if ( arg2.size() > 6  &&  arg2.substr( arg2.size() - 6, 6 ) == ".trans" ) { 
						trans_files.push_back ( arg1 + "/" + arg2 ) ;
					} else if ( arg2.size() == 6  &&  arg2.substr( 0, 5 ) == "frame" ) { 
						trans_files.push_back ( arg1 + "/" + arg2 ) ;
					} else if ( arg2.size() > 5  &&  arg2.substr( arg2.size() - 5, 5 ) == ".emit" ) { 
						emit_files.push_back ( arg1 + "/" + arg2 ) ;
					} else if ( arg2.size() == 4 &&  arg2.substr( 0, 4 ) == "emit" ) { 
						emit_files.push_back ( arg1 + "/" + arg2 ) ;
					}
				}
				sort( trans_files.begin(), trans_files.end() ) ;
				if ( emit_files.size() != 1 ) {
					cerr << arg1 << " : Multiple Emit files" << endl ; 
					exit(0) ;
				}
				for ( k = 1 ; k < trans_files.size() ; k++ ) { emit_files.push_back ( emit_files[0] ) ; } 
			} else {
				cerr << "Filed to open directory '" << arg1 << "'" << endl ;
				exit (0) ;
			}
			closedir( dir );

		} else {
			// TRANS and EMIT files are set for models, separately.
			trans_files.push_back ( arg1 ) ;
			if ( i + 1 >= argc ) { cerr << "EMISSION File of " << arg1 << " is not defined." << endl ;  exit ( 0 ) ; }
			emit_files.push_back ( argv[++i] ) ;

		}
	}
	if ( trans_files.size() == 0 ) {
		cerr << "No trans files are defined !!" << endl ;  
		exit ( 0 ) ; 
	}
	if ( emit_files.size() == 0 ) {
		cerr << "No emit files are defined !!" << endl ;  
		exit ( 0 ) ; 
	}
	if ( init_w_ratio < 1 || trans_w_ratio < 1 ) { 
		cerr << "Weight must be larger than 0 !!" << endl ;  
		exit ( 0 ) ; 
	} 
	
	if ( header_off == 0 ) { 
		for ( k = 0 ; k < trans_files.size() ; k++ ) {
			cout << "#HMM No.=" << k << "\tTransition Matrix File : " << trans_files[k] << endl ;
			cout << "#         \t  Emission Matrix File : " << emit_files[k] << endl ;
		}
		cout << endl ;
	}


	vector<HMM> hmm_vec ;
	vector<HMM> hmm_vec_tmp ;
	hmm_vec.resize( trans_files.size() ) ;
	
	for ( i = 0 ; i < trans_files.size() ; i++ ) {
		hmm_vec[i].setHMMparameters ( trans_files[i] , emit_files[i] ) ;
		//hmm_vec[i].printHMMparameters () ; 
	}
	
	Viterbi  vd ;
	ReadFile rf ;
	SEQ seq ;

	// ==========  Weight Matrix and Normalization  ==========
	vector<double> init_weight ;
	vector<vector<double> > trans_weight ;
	init_weight.resize( hmm_vec.size() ) ;
	trans_weight.resize( hmm_vec.size() ) ;
	for ( i = 0 ; i < hmm_vec.size() ; i++ ) {
		init_weight[i] = 1 ;
		trans_weight[i].resize( hmm_vec.size() ) ;
		for ( k = 0 ; k < hmm_vec.size() ; k++ ) {
			trans_weight[i][k] = 1 ;
		}
	}
	init_weight[0] = init_w_ratio ;
	for ( i = 0 ; i < hmm_vec.size() ; i++ ) {
		trans_weight[i][i] = trans_w_ratio ;
	}

	if ( header_off == 0 ) { 
		cout << "# Initiation Weight : " << init_w_ratio << " (BEGIN -> Others)" << endl ; 
	}
	if ( hmm_weights.size() > 0 ) {
		if ( header_off == 0 ) { 
			cout << "# Transition Weight : Manually-defined" << endl ;
			cout << endl ;
		}
		for ( i = 0 ; i < hmm_weights.size() ; i++ ) {
			if ( pre_hmm_nums[i] < 0 ||  pre_hmm_nums[i] >= hmm_vec.size() ) {
				cerr << "HMM No." << pre_hmm_nums[i] << "is NOT defined." << endl ;     exit ( 0 ) ;
			}
			if ( hmm_nums[i] < 0 ||  hmm_nums[i] >= hmm_vec.size() ) {
				cerr << "HMM No." << hmm_nums[i] << "is NOT defined." << endl ;    exit ( 0 ) ;
			}
			if ( hmm_weights[i] < 0 ) { cerr << "HMM weight " << hmm_weights[i] << "is incorrect." << endl ;   exit ( 0 ) ; }
			trans_weight[ pre_hmm_nums[i] ][ hmm_nums[i] ] = hmm_weights[i]  ;
			if ( weight_flags[i] > 0 ) {
				trans_weight[ hmm_nums[i] ][ pre_hmm_nums[i] ] = hmm_weights[i]  ;
			}
		}
	} else if ( header_off == 0 ) { 
		cout << "# Transition Weight : " << trans_w_ratio << endl ;
		cout << endl ;
	}

	vd.weightNormalization ( init_weight, trans_weight ) ;
	if ( subtract_flag == 0 ) { 
		vd.connectStates ( hmm_vec, init_weight, trans_weight ) ;
	} else if ( subtract_flag == 1 ) { 
		hmm_vec_tmp.resize( trans_files.size() ) ;
		k = 0 ;
		for ( i = 0 ; i < hmm_vec_tmp.size() ; i++ ) {
			hmm_vec_tmp[i].reserve_hmm(hmm_vec[i].getState()) ;		// Allocation of possible states
			k += hmm_vec[i].getState() - 1 ;
		}
		vd.reserve_hmm(k+1) ;	// Allocation of all possible states (including BEGIN)
	}

	//vd.printHMMparameters () ;

	rf.readFASTA ( seq_file ) ;
	if ( header_off == 0 ) { 
		vd.printWeight ( init_weight, trans_weight ) ;
	}
		cout << "# Sequence file : " << seq_file << endl ;
		if ( uncap == 1 ) { cout << "# Uncap option is used." << endl ; }
		cout << endl ; 
		//cout << endl << "# Sequence File : " << seq_file << endl << endl ; 

		cout << "# Sequence_name\tStart_pos\tEnd_pos\tStopped_or_skipped\tTransit_path_length\tProb_log_score\tHT_index" ;
		if ( show_stateseq == 1 ) { cout << "\tState_sequence" ; }
		if ( show_allprobs == 1 ) { cout << "\tAll_log_probs" ; }
		cout << endl ;
//	}
	
	string  each_seq ;
	vector<string> seq_vec ;
	//cout.precision(4) ;
	for ( i = 0 ; i < rf.seq_count() ; i++ ) {
		each_seq = rf.getSeq ( i ) ;
		if ( each_seq.size() % 3 != 0 ) { 
			cerr << rf.getName( i ) << " : Sequence length cannot be devided by 3.  So skipped." << endl ; 
			continue ; 
		}

		//cout << ">" << rf.getName( i ) << endl ;
		seq_vec.clear() ;
		seq_vec = seq.split ( each_seq , state_len , cut_len - state_len ) ;

		if ( subtract_flag == 1 ) { 
			for ( k = 0 ; k < hmm_vec.size() ; k++ ) {
				frame = hmm_vec[k].getFrame() ;
				hmm_vec_tmp[k] = hmm_vec[k] ;
				hmm_vec_tmp[k].subtractCounts ( each_seq, frame, cut_len, state_len ) ;
			}
			vd.connectStates ( hmm_vec_tmp, init_weight, trans_weight ) ;
			//vd.printHMMparameters () ;
		}
 
		vd.calcViterbi ( seq_vec , 0 , seq_vec.size() - 1 ) ;
		if ( uncap == 1 ) { vd.uncap() ; }
		cout << rf.getName( i ) << "\t" ; 
		cout << vd.getSummary (show_stateseq) ;
    	if ( show_allprobs == 1 ) { cout << "\t" << vd.getAllProbs ( seq_vec ) ; } 
		cout << endl ;
	}

	return 1 ;
}


