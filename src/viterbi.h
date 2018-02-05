// Class declarations for Viterbi computation

#ifndef  VITERBI_H
#define  VITERBI_H

class Viterbi : public HMM {

	std::vector<int> viPath ;	// Path p -> State n
	double  dProb ;		// Probability(Score) of Viterbi Path
	int  iFirst_path ;		// Number of the first positiion in the path computed
	int  iLast_path ;		// Number of the last position in the path computed
	int  iPath ;
	std::vector<double> vdInit_prob ;
	std::vector<std::string> vsStopped ;	// Information about transition-stopped or skipped positions

    public :
	Viterbi() ; 
	~Viterbi() ; 
	virtual  void  calcViterbi ( std::vector<std::string>& ,  int , int ) ;
	virtual  std::string  getSummary ( int ) ;
	virtual  std::string  getStateSequence () ;
	virtual  std::string  getAllProbs ( std::vector<std::string>& ) ;
	virtual  void  uncap () ;
	virtual  int  getPathLen () ; 
	//virtual  void  statistics () ; 

} ;

#endif

