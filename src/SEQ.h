#ifndef  SEQ_H
#define SEQ_H

class SEQ
{	
	
	public:

	// Reverse of Sequence
	virtual std::string reverse( std::string& );

	// Complement DNA Sequence (Non-reverse!)
	virtual std::string complement( std::string& );

	// Complement DNA Sequence (Non-reverse!)
	virtual std::string reverse_complement( std::string& );

	// Encode DNA sequence into 10 digit Number
	virtual int encode( std::string str );

	// Decode 10 digit number into DNA sequence (sequence length == len)
	virtual std::string decode( int code, int len );

	// Split Sequence
	virtual std::vector<std::string>  split ( std::string& , int , int ) ;

	// Join Sequence
	virtual std::string  join ( std::vector<std::string>& ) ;

	// Get SNP Sequence ( Only 1 Nucleotide Different )
	virtual std::vector<std::string>  getSNP ( std::string ) ;


};


#endif

