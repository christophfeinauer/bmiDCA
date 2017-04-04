#define BOOST_DISABLE_ASSERTS
#define DNDEBUG
#include "./bioio/bioio.hpp"
#include "boost/multi_array.hpp"	
#include <functional>
#include <map>
#include "msastruct.hpp"

struct ReadError : public std::runtime_error{
	ReadError(std::string const & message)
		: std::runtime_error(message)
	{}
};

MSAStruct::MSAStruct(std::string& fn, 
		double theta, 
		std::string charmap_string,
		char na_char) 
		: 
		msa_ptr(nullptr),
		w_ptr(nullptr), 
		na_char(na_char){

		if (na_char != '\0') na_exist=true;
		if (charmap_string != "") set_charmap(charmap_string); 
		if (charmap_string != "") letter2num = [this](char c){return letter2num_charmap(c);};
		else letter2num = [this](char c){return this->letter2num_aa_fasta(c);}; 
		if (charmap_string!=""){q = charmap_string.size(); default_index=q-1; };
		if (charmap_string!="") print_charmap();  
		
		read_msa(fn); 
		reweight(theta); 
};

void MSAStruct::get_frequencies(tens3_ptr& f3tens_ptr, tens2_ptr& f2tens_ptr, double theta, std::string storage_order){
	if (msa_ptr==nullptr)
		throw ReadError("Trying to calculate frequencies but not MSA read");
	if (w_ptr==nullptr)
		reweight(theta);
	Dvec& w = *w_ptr;
	MSA& msa = *msa_ptr;
	llu lenbn2 = (len*(len-1))/2;
	if (f3tens_ptr==nullptr){
		if (storage_order=="fortran") f3tens_ptr.reset(new tens3(boost::extents[q][q][lenbn2],boost::fortran_storage_order()));
		else f3tens_ptr.reset(new tens3(boost::extents[q][q][lenbn2]));
		std::fill_n((*f3tens_ptr).data(), (*f3tens_ptr).num_elements(), 0.0);
	}
	if (f2tens_ptr==nullptr){
		if (storage_order=="fortran") f2tens_ptr.reset(new tens2(boost::extents[q][len],boost::fortran_storage_order()));
		else f2tens_ptr.reset(new tens2(boost::extents[q][len]));
		std::fill_n((*f2tens_ptr).data(), (*f2tens_ptr).num_elements(), 0.0);
	}
	tens3& f3tens = *f3tens_ptr;
	tens2& f2tens = *f2tens_ptr;
	
	llu l=0;
	double x = 0.0;
	for (int m=0; m<M; ++m,l=0)
		for (int i=0; i<len; ++i){
			f2tens[ msa[m][i] ][ i ] += w[m];
			for (int j=i+1; j<len; ++j)
				f3tens[ msa[m][i] ][ msa[m][j] ][l++] += w[m];
	}

	// Normalize
	double Meff = get_Meff();
	double* f3tens_data = f3tens_ptr->data();
	double* f2tens_data = f2tens_ptr->data();
	for (llu l=0; l<f3tens.num_elements(); ++l)
		f3tens_data[l]/=Meff;
	for (llu l=0; l<f2tens.num_elements(); ++l)
		f2tens_data[l]/=Meff;
		
} 

void MSAStruct::set_charmap(std::string charmap_string){

	for (int i=0; i<charmap_string.size(); ++i)
		charmap[charmap_string[i]] = i+1;

}

void MSAStruct::print_charmap(){

	printf("Charmap: ");
	std::map<char,int>::iterator it = charmap.begin();
	printf("%c -> %i",it->first,it->second-1);
	++it;
	while (it!=charmap.end()){
		printf(", %c -> %i",it->first,it->second-1);
		++it;
	}
	printf(", default->%d\n",default_index);

}

int MSAStruct::letter2num_aa_fasta(char c){ 

	if (!isupper(c) && c!='-'){
		std::stringstream ss;
		ss << "Cannot parse character: " << c;
		throw ReadError(ss.str());
	}

	switch(c){
		case 'A': return  0;
		case 'B': return 20;
		case 'C': return  1;
		case 'D': return  2;
		case 'E': return  3;
		case 'F': return  4;
		case 'G': return  5;
		case 'H': return  6;
		case 'I': return  7;
		case 'J': return 20;
		case 'K': return  8;
		case 'L': return 9;
		case 'M': return 10;
		case 'N': return 11;
		case 'O': return 20;
		case 'P': return 12;
		case 'Q': return 13;
		case 'R': return 14;
		case 'S': return 15;
		case 'T': return 16;
		case 'U': return 20;
		case 'V': return 17;
		case 'W': return 18;
		case 'X': return 20;
		case 'Y': return 19;	
		default: return 20;
	}
}

void MSAStruct::reweight(double theta){
	printf("Reweighting: M=%lu, Cutoff=%lf",M,theta);
	fflush(stdout);
	if (msa_ptr==nullptr)
		throw ReadError("Cannot do reweighting before reading msa");
	MSA& msa = *msa_ptr;
	w_ptr.reset(new std::vector<double>);
	Dvec& w = *w_ptr;
	w.resize(M);
	std::fill(w.begin(),w.end(),1.0);

	int same;	
	int threshold = (int) ceil((1.0-theta) * len);
	int m1,m2,i;
	
	for (m1=0; m1<M; ++m1)
		for (m2=m1+1; m2<M; ++m2){
			same=0;
			for (i=0; i<len; ++i)
				if (msa[m1][i]==msa[m2][i]) same++;
			if (same > threshold) {w[m1]++; w[m2]++;};
		}

	for (int m=0; m<M; ++m)
		w[m]=1.0/w[m];

	printf(", Meff=%lf\n",get_Meff());
}

double MSAStruct::get_Meff(){
	if (w_ptr==nullptr)
		throw ReadError("Cannot compute Meff with weights not set");
	double Meff = 0.0;		
	Dvec& w = *w_ptr;
	for (int m=0; m<M; ++m)
		Meff+=w[m];
	return Meff;
}

int MSAStruct::letter2num_charmap(char c){ 
	int r = charmap[c];
	return r==0 ? default_index : r-1;
}


void MSAStruct::add_sequence(int m, int len, std::string& seq, MSA& msa){
	char c;
	for (int i=0; i<len; ++i){
		c = seq[i];	
		if (na_exist && c==na_char){
			(*na_ptr)[m][i]=true;
			msa[m][i] = default_index;	
			continue;
		}
		msa[m][i] = letter2num(seq[i]);
	}

}

void MSAStruct::make_na(){

	printf("Using NA character %c\n",na_char);
	na_ptr.reset(new BinaryMSA(boost::extents[M][len]));
	std::fill(na_ptr->data(), na_ptr->data()+na_ptr->num_elements(), false);
}

llu MSAStruct::count_na(){

	bool* pdata = na_ptr->data();
	llu c = 0;
	for (llu i=0; i<na_ptr->num_elements(); ++i)
		if (pdata[i]) ++c;
	return c;
}

void MSAStruct::read_msa(std::string& fn){

	std::ifstream ifs;
	M = bioio::count_fasta_records(fn);	

	ifs.open(fn);

	std::string seq = bioio::read_fasta_seqs(ifs,1)[0];
	len = seq.size();
	if (na_char != '\0') make_na();
	msa_ptr.reset(new MSA(boost::extents[M][len]));
	MSA& msa = *msa_ptr;
	msa[0][0]=0;
	add_sequence(0,len,seq,msa);
	for (size_t m =	1; m<M; ++m){
		seq = bioio::read_fasta_seqs(ifs,1)[0];
		if ( len != seq.size() ) throw ReadError("Sequence lengths inconsistent");	
		add_sequence(m,len,seq,msa);
	}

	ifs.close();

	if (na_char != '\0') printf("Fraction NA: %lf\n",(double) count_na() / (double) (na_ptr->num_elements()));
}


