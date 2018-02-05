// Class declarations for HMM

#ifndef HIDDEN_MARKOV_H
#define HIDDEN_MARKOV_H

const  double  FLOOR_SCORE = -100.0 ;		// Log score

class HMM {

    protected :
	int iState_n ;								// Number of states
	int iFrame ;	// Reading frame 

	std::vector<int> viEmit_total ;	
	std::vector<int> viTrans_total ;					

	std::vector<std::vector<int> >  vviTrans_state ;				// State n , Possible next state -> State n+1
	std::vector<std::vector<double> >  vvdTrans_prob ;				// State n , Possible next state -> Probability(as Log Score)
	std::vector<std::vector<double> >  vvdEmit_prob ;				// State n , Emission k -> Probability(as Log Score)
	std::vector<std::vector<int> >  vviTrans_count ;				// State n , Possible next state -> Count
	std::vector<std::vector<int> >  vviEmit_count ;				// State n , Emission k -> Count
	
	std::vector<std::string>  vsState ;							// State n -> State_Identifier
	std::vector<std::vector<std::string> >  vvsSymbol ;				// State n , Emission k  -> Emitted_Symbol


    public :
	HMM() ;
	~HMM() ;

	virtual  int  getState () ;
	virtual  int  getNextState ( int , int ) ;
	virtual  std::string  getSymbol ( int , int ) ;
	virtual  std::vector<std::vector<int> > getTransState() ;
	virtual  std::vector<std::vector<double> > getTransProbVector() ;
	virtual  std::vector<std::vector<double> > getEmitProbVector() ;

	virtual  std::vector<std::vector<int> > getTransCountVector() ;
	virtual  std::vector<int> getTransTotal() ;
	virtual  std::vector<std::vector<int> > getEmitCountVector() ;
	virtual  std::vector<int> getEmitTotal() ;
	virtual  std::vector<std::string> getStateVector() ;
	virtual  std::vector<std::vector<std::string> > getSymbolVector() ;


	virtual  void  resize_hmm ( int ) ;	
	virtual  void  reserve_hmm ( int ) ; 
					
	virtual  void  setHMMparameters ( std::string , std::string ) ;
	virtual  void  printHMMparameters () ;
	virtual  void  weightNormalization ( std::vector<std::vector<double> >& ) ;
	virtual  void  weightNormalization ( std::vector<double>& , std::vector<std::vector<double> >& ) ;
	virtual  void  printWeight ( std::vector<std::vector<double> >& ) ;
	virtual  void  printWeight ( std::vector<double>&, std::vector<std::vector<double> >& ) ;

	virtual  void  convert_to_logscore () ;
	virtual  void  convert_to_probability () ;
	virtual  void  subtractCounts ( std::string&, int, int, int ) ;
	virtual  int  getFrame () ;
	virtual  void  copy(HMM&) ;

	virtual  void  outputTransMatrix ( std::string& ) ;
	virtual  void  outputEmitMatrix ( std::string& ) ;

	virtual void  connectStates ( std::vector<HMM>&, std::vector<std::vector<double> >& ) ;
	virtual void  connectStates ( std::vector<HMM>&, std::vector<double>&, std::vector<std::vector<double> >& ) ;
} ;	

#endif

