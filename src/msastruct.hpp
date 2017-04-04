typedef boost::multi_array<int,2> MSA;
typedef std::shared_ptr<MSA> MSA_ptr;
typedef boost::multi_array<bool,2> BinaryMSA; // binary MSA to contain na character information
typedef std::shared_ptr<BinaryMSA> BinaryMSA_ptr;
typedef std::vector<double> Dvec;
typedef std::shared_ptr<Dvec> Dvec_ptr;
typedef long long unsigned llu;
typedef boost::multi_array<double,3> tens3;
typedef boost::multi_array<double,2> tens2;
typedef std::shared_ptr<tens3> tens3_ptr;
typedef std::shared_ptr<tens2> tens2_ptr;
#include <map>

struct MSAStruct{
	MSA_ptr msa_ptr = nullptr;
	BinaryMSA_ptr na_ptr = nullptr;
	size_t M,len;
	size_t q = 21;
	int default_index = 20;
	Dvec_ptr w_ptr = nullptr;
	MSAStruct() : msa_ptr(nullptr), w_ptr(nullptr) {};
	char na_char = '\0';
	bool na_exist = false;
	llu count_na();
	//std::map<char,int> charmap;
	//MSAStruct(std::string& fn, double theta) : msa_ptr(nullptr),w_ptr(nullptr) {read_msa(fn); reweight(theta); letter2num = std::bind(&MSAStruct::letter2num_aa_fasta,this,std::placeholders::_1);};
	MSAStruct(std::string& fn, double theta, std::string charmap_string = "", char na_char = '\0'); 
	void read_msa(std::string&);
	void set_charmap(std::string);
	void print_charmap();
	std::map<char,int> charmap;
	void reweight(double);
	int letter2num_aa_fasta(char);
	int letter2num_charmap(char);
	std::function<int(char)> letter2num;
	void add_sequence(int,int,std::string&,MSA& msa);
	double get_Meff();
	void get_frequencies(tens3_ptr&, tens2_ptr&, double, std::string="fortran");
	void make_na(); // makes the array that contains information about na characters
	
};

