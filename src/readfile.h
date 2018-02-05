// Header File for Read Data (FASTA format etc.)

#ifndef  READFILE_H
#define READFILE_H

//const  int  BUFFER_SIZE = 1024 ;
//const  int  BUFFER_SIZE = 4096 ;
const  int  BUFFER_SIZE = 16384 ;
//const  int  BUFFER_SIZE = 262144 ;

class ReadFile {

    private:
	std::string filename ;
	int  max_line_len ;
	std::vector<std::map<std::string, std::string> > seq_data ;
	std::map<std::string, int>  hash_keys ;

	//std::vector<std::string> sequence;
	//std::vector<std::string> seq_name;

    public:

	// Constructor
	//SEQ (std::string f);
	//SEQ ();

	// ==========  Member functions  ==========
	int  readDataLine ( std::string, std::vector<std::string>& ) ;
	int  readDataLine ( std::string, std::vector<std::string>&, std::vector<std::string>& ) ;	// obtain comments,too

	int  readFASTA( std::string ) ;		// Cumulatively add the data into "seq_data"
	void clearData() ;			// Empty "seq_data"
	int  seq_count () ;
	int  max_line_length() ;	
	void show () ;	
	std::map<std::string, std::string> getSeqDatum ( int ) ;
	std::map<std::string, std::string>  getSeqDatumHash ( std::string& ) ;
	//std::vector<std::string>  splitSeq ( int , int ) ;
	void  makeHash() ;
	
	std::string getSeq ( int );
	std::string getName( int );
	std::string getFullname( int );

} ;					


#endif

