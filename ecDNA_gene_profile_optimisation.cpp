
# include <iostream>
# include <fstream>
# include <sstream>
# include <stdlib.h>
# include <stdio.h>
# include <cmath>
# include <math.h>
# include <random>
# include <cstdlib>
# include <dirent.h>
# include <string>
# include <getopt.h>
# include <algorithm>
# include <initializer_list>
# include <vector>
# include <iterator>




using namespace std;



/*******************************************************************************/



// Define global variables
int seed, initial_copyNumber, Ntot, N_ecDNA_hot, iter, doubled_ecDNA_copyNumber, dying_cell_index, daughter_ecDNA_copyNumber1, daughter_ecDNA_copyNumber2, daughter_1_current_index, daughter_2_current_index, Nmax, occupancy, max_ecDNA_size, number_of_oncogenes, num_fusions, num_fissions, current_loop_index, partition_start, partition_end, ecDNA_count, cellStartingCopyNumber;
double t, dt, time_at_confluency, r_birth_normalised, r_death_normalised, total_unnormalised_division_rate, total_unnormalised_death_rate, total_unnormalised_division_rate_logisticRescaled, rand_double, cumulative_division_rate , cumulative_death_rate, selection_coeff, sigmoid_a, sigmoid_b, ec_length_sum, x1, x2, x3, x4, ecDNA_size_multiplier, ecDNA_size_multiplier_factor, fusion_probability, fission_probability, partition_prob, base_rate;
bool verbose_flag, BIRTH, DEATH, added_to_occupancy_vector, allocated_to_daughter_1, counted;
string new_ecDNA_postFusion, new_ecDNA_postFission;
vector<string> daughter_1_ec, daughter_2_ec, mother_ec;
vector<int> daughter_1_ec_indices, mother_cell_indices, all_cell_indices, resampling_indices;

const char gene_A_char = 'A';
const char gene_B_char = 'B';


/*******************************************************************************/



// Define a cell. Genes are denoted A and B, with A being the oncogene and B being the passenger
class Cell
{
	public:
		vector<string> ecDNA;
		double division_rate;
		double death_rate;

		// Constructors for Cell object
		Cell()
		{
			vector<string> ec;
			set_ecDNA(ec);
			set_division_rate(ec);
			set_death_rate(ec);
		}

		Cell(vector<string> ec)
		{
			set_ecDNA(ec);
			set_division_rate(ec);
			set_death_rate(ec);
		}

		// Set() and get() methods
		
		void set_ecDNA(vector<string> ec)
		{
			this->ecDNA = ec;
		}

		void set_division_rate(vector<string> ec)
		{
			// Sigmoid selection function
			// Compute number of oncogenes (gene A) in cell
			base_rate = log(2.0);		// Corresponds to a doubling time of 1 day
			number_of_oncogenes = 0;
			for (int i = 0; i < ec.size(); ++i)	// Loop over each ecDNA in cell
			{
				for (int j = 0; j < ec[i].size(); ++j)	// Loop through each gene on ecDNA
				{
					if (ec[i][j] == gene_A_char) number_of_oncogenes += 1;
				}
			}

			if (number_of_oncogenes == 0) 
			{
				this->division_rate = base_rate;
			}
			else if (number_of_oncogenes >= sigmoid_b)
			{
				this->division_rate = (1.0 + selection_coeff) * base_rate;
			}
			else
			{
				x1 = selection_coeff;
				x2 = 1.0 + ((sigmoid_b - (double)number_of_oncogenes)/(sigmoid_b - sigmoid_a));
				x3 = ((double)number_of_oncogenes / sigmoid_b);
				x4 = (sigmoid_b / (sigmoid_b - sigmoid_a));

				this->division_rate = (1.0 + (x1 * x2 * (pow(x3 , x4)))) * base_rate;
			}

			// Constant selection function
			// if (ec.size() == 0) this->division_rate = 1.0;
			// else this->division_rate = 1.0 + selection_coeff;
		}


		void set_death_rate(vector<string> ec)
		{
			base_rate = log(2.0);		// Corresponds to a doubling time of 1 day

			// Base death rate varies between 0.5->1.0 times the cell's birth rate, depending on total ecDNA burden
			ec_length_sum = 0.0;
			for (int i = 0; i < ec.size(); ++i) ec_length_sum += (double)ec[i].size();

			this->death_rate = (0.5 * base_rate) + (ecDNA_size_multiplier_factor/100.0*ec_length_sum);

		}
};








//-----------------------







// Define an element of the occupancy vector
class Occupancy
{
	public:
		string ecDNA;
		int multiplicity;

		// Constructors for occupancy object
		Occupancy(string i , int j)
		{
			set_ecDNA(i);
			set_multiplicity(j);
		}

		// Set() and get() methods
		
		void set_ecDNA(string n)
		{
			this->ecDNA = n;
		}

		void set_multiplicity(int n)
		{
			this->multiplicity = n;
		}


};





//-----------------------






// Define a tuple for iterting over A and B pairs 
class Pair
{
	public:
		double a;
		double b;

		
		Pair()
		{
			set_a(0.0);
			set_b(0.0);
		}

		Pair(double x , double y)
		{
			set_a(x);
			set_b(y);
		}

		void set_a(double x)
		{
			this->a = x;
		}

		void set_b(double y)
		{
			this->b = y;
		}


};






//-----------------------







// Select cell to divide, based on individual division rates
void cell_division(vector<Cell> &tissue , double *total_unnormalised_division_rate , double *total_unnormalised_death_rate , int *Ntot , int *N_ecDNA_hot , double fusion_probability , double fission_probability , mt19937_64 *generator)
{

	cumulative_division_rate = 0.0;

	rand_double = drand48();

	for (int i = 0; i < *Ntot; ++i)
	{

		cumulative_division_rate += tissue[i].division_rate;

		if (rand_double <= cumulative_division_rate/(*total_unnormalised_division_rate))
		{

			// cout << "\n\n" << endl;

			// cout << "Mother cell ecDNA (pre-ecDNA replication): ";
			// for (int j = 0; j < tissue[i].ecDNA.size(); ++j)
			// {
			// 	cout << tissue[i].ecDNA[j] << ", ";
			// }
			// cout << endl;


			// mother_ec is a full, explicit list of ecDNA in the mother cell after ecDNA replication
			mother_ec.clear();
			mother_ec.insert(mother_ec.begin(), tissue[i].ecDNA.begin(), tissue[i].ecDNA.end() );
			mother_ec.insert(mother_ec.begin(), tissue[i].ecDNA.begin(), tissue[i].ecDNA.end() );


			// cout << "Mother cell ecDNA (post-ecDNA replication): ";
			// for (int j = 0; j < mother_ec.size(); ++j)
			// {
			// 	cout << mother_ec[j] << ", ";
			// }
			// cout << endl;




			// Evolve ecDNA. First simulate fusions of ecDNAs, then fissions
			// Fusion
			vector<string> mother_ec_postFusion;
			while (1)
			{
				// End when all ecDNAs from mother_ec have been allocated to mother_ec_postFissionAndFusion
				if (mother_ec.size() == 0) break;

				new_ecDNA_postFusion = mother_ec[0];
				mother_ec.erase(mother_ec.begin());


				// draw poisson random number for number of ecDNA this ecDNA groups (fuses) with
				poisson_distribution<int> draw_ecDNA_fusion_number(fusion_probability);
				num_fusions = draw_ecDNA_fusion_number(*generator);


				for (int j = 0; j < min((int)mother_ec.size() , num_fusions); ++j)	// Perhaps should revise how I deal with edge cases here
				{
					new_ecDNA_postFusion += mother_ec[0];
					mother_ec.erase(mother_ec.begin());
				}

				// Add fused ecDNA to new vector
				mother_ec_postFusion.push_back(new_ecDNA_postFusion);

			}


			// cout << "Mother cell ecDNA (after fusion): ";
			// for (int j = 0; j < mother_ec_postFusion.size(); ++j)
			// {
			// 	cout << mother_ec_postFusion[j] << ", ";
			// }
			// cout << endl;


			// Fission
			vector<string> mother_ec_postFissionAndFusion;
			for (int j = 0; j < mother_ec_postFusion.size(); ++j)
			{

				// number of breaks in ecDNA is proportional to its length
				binomial_distribution<int> draw_ecDNA_fission_number(mother_ec_postFusion[j].size()-1 , fission_probability);
				num_fissions = draw_ecDNA_fission_number(*generator);

				//num_fissions = min((int)mother_ec_postFusion[j].size()-1 , num_fissions);


				// Determine where to partition the ecDNA 
				partition_prob = 1.0/(double)(mother_ec_postFusion[j].size()-1);
				vector<int> partition_indices;
				current_loop_index = 0;
				while(1)
				{
					if (partition_indices.size() == num_fissions) break;
					//cout << "Finding " << num_fissions << " partitions for ecDNA " << mother_ec_postFusion[j] << ". current_loop_index=" << current_loop_index << ". partition_indices.size()=" << partition_indices.size() << endl;

					if ((drand48() < partition_prob) && (count(partition_indices.begin(), partition_indices.end(), current_loop_index) == 0)) partition_indices.push_back(current_loop_index);

					current_loop_index = (current_loop_index + 1)%(int)(mother_ec_postFusion[j].size()-1);

				}


				// Sort partition indices in ascending order
				sort(partition_indices.begin(), partition_indices.end());

				// if (partition_indices.size() > 0)
				// {
				// 	cout << "Partition indices: ";
				// 	for (int j = 0; j < partition_indices.size(); ++j)
				// 	{
				// 		cout << partition_indices[j] << endl;
				// 	}
				// }



				partition_start = 0;
				for (int k = 0; k < partition_indices.size(); ++k)
				{
					// Separate ecDNA
					partition_end = partition_indices[k] + 1;
					new_ecDNA_postFission = "";
					for (int l = partition_start; l < partition_end; ++l) new_ecDNA_postFission += mother_ec_postFusion[j][l];

					// Add separted ecDNA to new vector
					mother_ec_postFissionAndFusion.push_back(new_ecDNA_postFission);


					partition_start = partition_indices[k] + 1;
				}


				partition_end = mother_ec_postFusion[j].size();
				new_ecDNA_postFission = "";
				for (int k = partition_start; k < partition_end; ++k) new_ecDNA_postFission += mother_ec_postFusion[j][k];

				// Add separted ecDNA to new vector
				mother_ec_postFissionAndFusion.push_back(new_ecDNA_postFission);



			}








			// // Fission (no prior fusion)
			// vector<string> mother_ec_postFissionAndFusion;
			// for (int j = 0; j < mother_ec.size(); ++j)
			// {

			// 	// number of breaks in ecDNA is proportional to its length
			// 	binomial_distribution<int> draw_ecDNA_fission_number(mother_ec[j].size()-1 , fission_probability);
			// 	num_fissions = draw_ecDNA_fission_number(*generator);

			// 	//num_fissions = min((int)mother_ec_postFusion[j].size()-1 , num_fissions);


			// 	// Determine where to partition the ecDNA 
			// 	partition_prob = 1.0/(double)(mother_ec[j].size()-1);
			// 	vector<int> partition_indices;
			// 	current_loop_index = 0;
			// 	while(1)
			// 	{
			// 		if (partition_indices.size() == num_fissions) break;
			// 		//cout << "Finding " << num_fissions << " partitions for ecDNA " << mother_ec_postFusion[j] << ". current_loop_index=" << current_loop_index << ". partition_indices.size()=" << partition_indices.size() << endl;

			// 		if ((drand48() < partition_prob) && (count(partition_indices.begin(), partition_indices.end(), current_loop_index) == 0)) partition_indices.push_back(current_loop_index);

			// 		current_loop_index = (current_loop_index + 1)%(int)(mother_ec[j].size()-1);

			// 	}


			// 	// Sort partition indices in ascending order
			// 	sort(partition_indices.begin(), partition_indices.end());

			// 	// if (partition_indices.size() > 0)
			// 	// {
			// 	// 	cout << "Partition indices: ";
			// 	// 	for (int j = 0; j < partition_indices.size(); ++j)
			// 	// 	{
			// 	// 		cout << partition_indices[j] << endl;
			// 	// 	}
			// 	// }



			// 	partition_start = 0;
			// 	for (int k = 0; k < partition_indices.size(); ++k)
			// 	{
			// 		// Separate ecDNA
			// 		partition_end = partition_indices[k] + 1;
			// 		new_ecDNA_postFission = "";
			// 		for (int l = partition_start; l < partition_end; ++l) new_ecDNA_postFission += mother_ec[j][l];

			// 		// Add separted ecDNA to new vector
			// 		mother_ec_postFissionAndFusion.push_back(new_ecDNA_postFission);


			// 		partition_start = partition_indices[k] + 1;
			// 	}


			// 	partition_end = mother_ec[j].size();
			// 	new_ecDNA_postFission = "";
			// 	for (int k = partition_start; k < partition_end; ++k) new_ecDNA_postFission += mother_ec[j][k];

			// 	// Add separted ecDNA to new vector
			// 	mother_ec_postFissionAndFusion.push_back(new_ecDNA_postFission);



			// }





			



			// cout << "Mother cell ecDNA (after fission): ";
			// for (int j = 0; j < mother_ec_postFissionAndFusion.size(); ++j)
			// {
			// 	cout << mother_ec_postFissionAndFusion[j] << ", ";
			// }
			// cout << endl;



			// Every resulting ecDNA is copied once
			doubled_ecDNA_copyNumber = mother_ec_postFissionAndFusion.size();



			// Distribute ecDNA between two daughter cells according to binomial
			binomial_distribution<int> draw_new_ecDNA_copyNumber(doubled_ecDNA_copyNumber , 0.5);
			daughter_ecDNA_copyNumber1 = draw_new_ecDNA_copyNumber(*generator);
			daughter_ecDNA_copyNumber2 = doubled_ecDNA_copyNumber - daughter_ecDNA_copyNumber1;



			// Book-keeping
			*total_unnormalised_division_rate -= tissue[i].division_rate;
			*total_unnormalised_death_rate -= tissue[i].death_rate;



			// Allocate ecDNA to new daughter cells 
			daughter_1_ec.assign(daughter_ecDNA_copyNumber1 , "");
			daughter_2_ec.assign(daughter_ecDNA_copyNumber2 , "");




			// Define sequence of indices up to size of mother cell ecDNA vector 
			mother_cell_indices.assign(doubled_ecDNA_copyNumber , 0);
			for (int j = 0; j < doubled_ecDNA_copyNumber; ++j)
			{
				mother_cell_indices[j] = j;
			}



			// Randomly select n of these (n = daughter_ecDNA_copyNumber1)
			daughter_1_ec_indices.clear();
			sample(mother_cell_indices.begin(), mother_cell_indices.end(), back_inserter(daughter_1_ec_indices), daughter_ecDNA_copyNumber1, *generator);



			// Allocate ecDNA depending on if their index is contained within random selection vector 
			daughter_1_current_index = 0;
			daughter_2_current_index = 0;
			for (int k = 0; k < mother_ec_postFissionAndFusion.size(); ++k)
			{
				allocated_to_daughter_1 = false;
				for (int j = 0; j < daughter_1_ec_indices.size(); ++j)
				{
					if (daughter_1_ec_indices[j] == k)
					{
						daughter_1_ec[daughter_1_current_index] = mother_ec_postFissionAndFusion[k];
						daughter_1_current_index += 1;
						allocated_to_daughter_1 = true;
						break;
					}
				}

				// If index not in list of indices for daughter 1, ecDNA goes into daughter 2
				if (allocated_to_daughter_1 == false)
				{
					daughter_2_ec[daughter_2_current_index] = mother_ec_postFissionAndFusion[k];
					daughter_2_current_index += 1;
				}
			}


			// cout << "Daughter 1 cell ecDNA: ";
			// for (int j = 0; j < daughter_1_ec.size(); ++j)
			// {
			// 	cout << daughter_1_ec[j] << ", ";
			// }
			// cout << endl;



			// cout << "Daughter 2 cell ecDNA: ";
			// for (int j = 0; j < daughter_2_ec.size(); ++j)
			// {
			// 	cout << daughter_2_ec[j] << ", ";
			// }
			// cout << endl;


		


			// Create daughter cells
			tissue[i] = Cell(daughter_1_ec);
			tissue[*Ntot] = Cell(daughter_2_ec);



			// Book-keeping
			*total_unnormalised_division_rate += tissue[i].division_rate + tissue[*Ntot].division_rate;
			*total_unnormalised_death_rate += tissue[i].death_rate + tissue[*Ntot].death_rate;

			if ((daughter_ecDNA_copyNumber1 > 0) && (daughter_ecDNA_copyNumber2 > 0)) *N_ecDNA_hot += 1;
			*Ntot += 1;

			return;
		}
	}


	// If we get this far, there is a problem with the Gillespie rates. 
	cout << "Problem with Gillespie rates encountered when choosing cell to divide. Exiting..." << endl;
	exit(0);
	
}






//-----------------------






// Select cell for death, based on individual death rates
void cell_death(vector<Cell> &tissue , double *total_unnormalised_division_rate , double *total_unnormalised_death_rate , int *Ntot , int *N_ecDNA_hot , mt19937_64 *generator)
{

	cumulative_death_rate = 0.0;

	rand_double = drand48();

	for (int i = 0; i < *Ntot; ++i)
	{

		cumulative_death_rate += tissue[i].death_rate;

		if (rand_double <= cumulative_death_rate/(*total_unnormalised_death_rate))
		{

			// Book-keeping
			if (tissue[i].ecDNA.size() > 0) *N_ecDNA_hot -= 1;
			*total_unnormalised_division_rate -= tissue[i].division_rate;
			*total_unnormalised_death_rate -= tissue[i].death_rate;

			// Cell dies. Replace with last cell in tissue vector so that all vector indices from 0 to Ntot-1 are occupied by cells
			tissue[i] = Cell(tissue[*Ntot-1].ecDNA);
			*Ntot -= 1;
			
			return;
		}
	}


	// If we get this far, there is a problem with the Gillespie rates. 
	cout << "Problem with Gillespie rates encountered when choosing cell to die. Exiting..." << endl;
	exit(0);
	
}







//-----------------------








// Parse command line arguments (Flags and numerical arguments)
void parse_command_line_arguments(int argc, char** argv , bool *verbose_flag , int *seed , int *Nmax , int *initial_copyNumber , double *selection_coeff , double *fusion_probability , double *fission_probability , double *ecDNA_size_multiplier_factor)
{
	int c;
	int option_index;
	char* arg_long = nullptr;
	int verbose = 0;

	static struct option long_options[] =
	{
		{"verbose", no_argument, &verbose, 1},
	}; 

	while ((c = getopt_long(argc, argv, "x:n:k:s:a:b:p:q:c:", long_options, &option_index)) != -1)
	switch (c)
	{
		case 0:
		{
			arg_long = optarg;
			break;
		}

		// Random seed
		case 'x':
			*seed = atoi(optarg);	
			break;


		// Selection coefficient
		case 'n':
			*Nmax = atoi(optarg);
			break;


		// ecDNA copy number in initial cell
		case 'k':
			*initial_copyNumber = atoi(optarg);
			break;


		// ecDNA copy number in initial cell
		case 's':
			*selection_coeff = atof(optarg);
			break;


		// // a parameter for sigmoid fitness function
		// case 'a':
		// 	*sigmoid_a = atof(optarg);
		// 	break;


		// // b parameter for sigmoid fitness function
		// case 'b':
		// 	*sigmoid_b = atof(optarg);
		// 	break;



		// parameter defining probability of ecDNA fusing 
		case 'p':
			*fusion_probability = atof(optarg);
			break;



		// parameter defining probability of ecDNA separating
		case 'q':
			*fission_probability = atof(optarg);
			break;


		// ecDNA burden penalty
		case 'c':
			*ecDNA_size_multiplier_factor = atof(optarg);	
			break;


		case '?':
			if (optopt == 'c')
				fprintf (stderr, "Option -%c requires an argument.\n", optopt);
			else if (isprint (optopt))
				fprintf (stderr, "Unknown option `-%c'.\n", optopt);
			else
				fprintf (stderr,
				"Unknown option character `\\x%x'.\n",
				optopt);
		//return 1;
		default:
		abort ();
	}



	// Set boolean values for verbose flag
	if (verbose == 1) *verbose_flag = true;



	if (*initial_copyNumber < 0)
	{
		cout << "Initial copy number must be 0 or greater. Exiting." << endl;
		exit(0);
	}


	if (*fusion_probability < 0)
	{
		cout << "Invalid ecDNA fusion probability. Exiting..." << endl;
		exit(0);
	}

	if (*fission_probability < 0)
	{
		cout << "Invalid ecDNA fission probability. Exiting..." << endl;
		exit(0);
	}

	// if (*sigmoid_a >= *sigmoid_b)
	// {
	// 	cout << "Must have a < b. Exiting..." << endl;
	// 	exit(0);
	// }
}






//-----------------------






// Set up tissue (i.e. array of cells)
vector<Cell> initialise_tissue(int Nmax , int initial_cell_number , double *total_unnormalised_division_rate , double *total_unnormalised_death_rate , int *Ntot , int *N_ecDNA_hot , int initial_copyNumber , mt19937_64 *generator)
{

	if (verbose_flag) cout << " " << endl;


	// Set up the vector of cells, called tissue
	vector<Cell> tissue(Nmax*1.5); 
	if (verbose_flag) printf(" Initialising tissue... Done.\r");
	if (verbose_flag) cout << " " << endl;


	// Sample from distribution of initial ecDNA copy number to create list of initial cells
	normal_distribution draw_ecDNA_copy_number{(double)(initial_copyNumber), (double)(initial_copyNumber)*0.5}; 	// Choose wide distribution for initial copy numbers
		

	for (int i = 0; i < initial_cell_number; ++i)
	{
		// Draw cell's ecDNA copy number from Gaussian distribution
		cellStartingCopyNumber = round(max(0.0 , draw_ecDNA_copy_number(*generator)));
		vector<string> initial_ec(cellStartingCopyNumber , "AB");

		// Seed first tissue cell
		tissue[i] = Cell(initial_ec);

		// Book-keeping
		*Ntot += 1;
		if (initial_ec.size() > 0) *N_ecDNA_hot += 1;
		*total_unnormalised_division_rate += tissue[i].division_rate;
		*total_unnormalised_death_rate += tissue[i].death_rate;
	}

	return tissue;
}






//-----------------------






// Compute normalised birth and death rates 
void compute_normalised_birth_and_death_rates(int Ntot , int N_ecDNA_hot , double total_unnormalised_division_rate , double total_unnormalised_death_rate , double *r_birth_normalised , double *r_death_normalised)
{

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	// Compute un-normalised reaction rates, constant death rate for all cells
	//total_unnormalised_death_rate = 0.5*Ntot;
	//total_unnormalised_death_rate = 0.0;

	// If you want cell death to be proportional to cell birth rate, compute the line below and invoke cell_death() method
	//total_unnormalised_death_rate = 0.5*total_unnormalised_division_rate;

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



	// Compute normalised reaction rates
	*r_birth_normalised = total_unnormalised_division_rate/(total_unnormalised_division_rate + total_unnormalised_death_rate);
	*r_death_normalised = total_unnormalised_death_rate/(total_unnormalised_division_rate + total_unnormalised_death_rate);

}






//-----------------------






// Choose next event in Gillespie algorithm
void choose_next_event(bool *BIRTH , bool *DEATH , double r_birth_normalised , double r_death_normalised)
{
	*BIRTH = false;
	*DEATH = false;

	rand_double = drand48();
	if (rand_double <= r_birth_normalised)
	{
		*BIRTH = true;
	}
	else if (rand_double <= r_birth_normalised + r_death_normalised)
	{
		*DEATH = true;
	}
	else
	{
		cout << "Problem with Gillespie rates..." << endl;
		exit(0);
	}
}








//-----------------------







// // Kill cell and remove from lattice
// void kill_cell(vector<Cell> &tissue , int *Ntot , int *N_ecDNA_hot , double *total_unnormalised_division_rate)
// {

// 	// Randomly select one cell to die (uniform probability)
// 	dying_cell_index = round((drand48() * (*Ntot)) - 0.5);

	
// 	// Book-keeping
// 	if (tissue[dying_cell_index].ecDNA > 0) *N_ecDNA_hot -= 1;
// 	*total_unnormalised_division_rate -= tissue[dying_cell_index].division_rate;

// 	// Cell dies. Replace with last cell in tissue vector so that all vector indices from 0 to Ntot-1 are occupied by cells
// 	tissue[dying_cell_index] = Cell(tissue[*Ntot-1].ecDNA);
// 	*Ntot -= 1;

// }











// /*******************************************************************************/















int main(int argc, char** argv)
{



	// Reset time and tissue size variables
	t = 0.0;
	Ntot = 0;
	N_ecDNA_hot = 0;
	total_unnormalised_division_rate = 0.0;
	total_unnormalised_death_rate = 0.0;

	double num_replatings = 1;
	int resampling_number = 20000;



	// Define list of A and B combinations to iterate over
	vector<Pair> ab_pairs(6);
	ab_pairs[0] = Pair(5,10);
	ab_pairs[1] = Pair(10,20);
	ab_pairs[2] = Pair(15,30);
	ab_pairs[3] = Pair(25,50);
	ab_pairs[4] = Pair(40,80);
	ab_pairs[5] = Pair(60,120);



	// Set initial A and B parameters
	sigmoid_a = ab_pairs[0].a;
	sigmoid_b = ab_pairs[0].b;



	//================== Parse command line arguments ====================//
	parse_command_line_arguments(argc , argv , &verbose_flag , &seed , &Nmax , &initial_copyNumber , &selection_coeff , &fusion_probability , &fission_probability , &ecDNA_size_multiplier_factor);


	if (Nmax < resampling_number)
	{
		cout << "Nmax cannot be lower than resampling number. Exiting..." << endl;
		exit(0);
	}


	cout << "Outputs for " << argv[0] << " -n " << Nmax << " -n_replate " << resampling_number << " -k " << initial_copyNumber << " -s " << selection_coeff << " -A " << sigmoid_a << " -B " << sigmoid_b << " -C " << ecDNA_size_multiplier_factor << " -p " << fusion_probability << " -q " << fission_probability << " -x " << seed << endl;
	cerr << "Errors for " << argv[0] << " -n " << Nmax << " -n_replate " << resampling_number << " -k " << initial_copyNumber << " -s " << selection_coeff << " -A " << sigmoid_a << " -B " << sigmoid_b << " -C " << ecDNA_size_multiplier_factor << " -p " << fusion_probability << " -q " << fission_probability << " -x " << seed << endl;


	// Seed random number generator
	srand48(seed);
	mt19937_64 generator;
	generator.seed(seed);





	//================== Initialise tissue ====================//
	vector<Cell> tissue = initialise_tissue(Nmax , resampling_number , &total_unnormalised_division_rate , &total_unnormalised_death_rate , &Ntot , &N_ecDNA_hot , initial_copyNumber , &generator);





	//================== Simulate tissue growth ==================//
	iter = 0;

	for (int repeat = 0; repeat < ab_pairs.size(); ++repeat)
	{
		// Set new A and B parameters for fitness function
		sigmoid_a = ab_pairs[repeat].a;
		sigmoid_b = ab_pairs[repeat].b;


		// Update birth rates and re-compute total unnormalised birth rates under new fitness curve
		total_unnormalised_division_rate = 0.0;
		for (int i = 0; i < Ntot; ++i)
		{
			tissue[i].set_division_rate(tissue[i].ecDNA);
			total_unnormalised_division_rate += tissue[i].division_rate;
		}

		//cout << repeat << endl;
		
		do
		{



			// Timestamp for simulation optimisation purposes
			//clock_t iter_start_time, iter_end_time;

			//iter_start_time = clock();

			//if (iter == 20) exit(0);

			//cout << "---- iter = " << iter << " ---------- n_ecDNA_hot = " << N_ecDNA_hot << " ----------------------------" << endl;

			
			++iter;



			// For logistic growth, modulate all death rates according to distance from carrying capacity (which is Nmax)
			//total_unnormalised_death_rate_logisticRescaled = total_unnormalised_death_rate + ((total_unnormalised_division_rate - total_unnormalised_death_rate) * (sqrt((double)Ntot/(double)Nmax)));
			total_unnormalised_division_rate_logisticRescaled = total_unnormalised_division_rate - ((total_unnormalised_division_rate - total_unnormalised_death_rate) * (sqrt((double)Ntot/(double)Nmax)));



			//cout << "Ntot = " << Ntot << ", Nmax = " << Nmax << " | Death rescaled " << total_unnormalised_death_rate << " -> " << total_unnormalised_death_rate_logisticRescaled << " | Total birth rate minus death rate per cell = " << (total_unnormalised_division_rate-total_unnormalised_death_rate_logisticRescaled)/(double)Ntot << endl;






			// Re-evaluate birth and death rates
			compute_normalised_birth_and_death_rates(Ntot , N_ecDNA_hot , total_unnormalised_division_rate_logisticRescaled , total_unnormalised_death_rate , &r_birth_normalised , &r_death_normalised);





			// Choose birth or death based on normalised rates
			choose_next_event(&BIRTH , &DEATH , r_birth_normalised , r_death_normalised);




			// If division:
			if (BIRTH)
			{
				// Cell divides
				cell_division(tissue , &total_unnormalised_division_rate , &total_unnormalised_death_rate , &Ntot , &N_ecDNA_hot , fusion_probability , fission_probability , &generator);
			}


			// If death:
			if ((DEATH) && (Ntot > 1))
			{
				// Cell dies (if constant death rate used)
				//kill_cell(tissue , &Ntot , &N_ecDNA_hot , &total_unnormalised_division_rate);

				// Cell dies (if CN dependent death rate used)
				cell_death(tissue , &total_unnormalised_division_rate , &total_unnormalised_death_rate , &Ntot , &N_ecDNA_hot , &generator);
			}







			if (iter%100000 == 0)
			{			
				//cout << t << "," << Ntot << endl; 

				if (verbose_flag) cout << "Iteration #" << iter << " -- N = " << Ntot << " -- N_ecDNA_hot = " << N_ecDNA_hot << " -- total time = " << t << " days --  time at confluency = " << time_at_confluency << endl;
			
				// If all ecDNA lost from population, exit Gillespie loop and write output data
				if (N_ecDNA_hot == 0)
				{
					//if (repeat%10 == 0)
					//{
						// stringstream f;
						// f.str("");
						// f << "./RESULTS/Nmax=" << Nmax << "_k=" << initial_copyNumber << "_s=" << selection_coeff << "_A=" << sigmoid_a << "_B=" << sigmoid_b << "_C=" << ecDNA_size_multiplier_factor << "_p=" << fusion_probability << "_q=" << fission_probability << "/seed=" << seed;
						// DIR *dir = opendir(f.str().c_str());
						// if(!dir)
						// {
						// 	f.str("");
						// 	f << "mkdir -p ./RESULTS/Nmax=" << Nmax << "_k=" << initial_copyNumber << "_s=" << selection_coeff << "_A=" << sigmoid_a << "_B=" << sigmoid_b << "_C=" << ecDNA_size_multiplier_factor << "_p=" << fusion_probability << "_q=" << fission_probability << "/seed=" << seed;
						// 	system(f.str().c_str());
						// }

						// ofstream tissue_file;
						// f.str("");
						// f << "./RESULTS/Nmax=" << Nmax << "_k=" << initial_copyNumber << "_s=" << selection_coeff << "_A=" << sigmoid_a << "_B=" << sigmoid_b << "_C=" << ecDNA_size_multiplier_factor << "_p=" << fusion_probability << "_q=" << fission_probability << "/seed=" << seed << "/tissue_resample" << repeat <<".csv";
						// tissue_file.open(f.str().c_str());


						// if (verbose_flag) cout << " " << endl;
						// if (verbose_flag) cout << "Created output files..." << endl;


						// // Loop over all cells 
						// for (int i = 0; i < tissue.size(); ++i)
						// {
							
						// 	tissue_file << "," << endl;
							
						// }
					//}


					exit(0);
				}
			}




			// Progress time variable
			dt = (-1.0 / (total_unnormalised_division_rate_logisticRescaled + total_unnormalised_death_rate)) * log( drand48() );
			t += dt;


			// Define confluency reached as 95% of carrying capacity
			if ((double)Ntot/(double)Nmax >= 0.95) time_at_confluency += dt;
			else time_at_confluency = 0.0;




		} while (time_at_confluency < 14.0);		// Exit once confluency maintained for at least 14 days








		//================== Open data files & write final system data ==================//

		stringstream f;
		f.str("");
		f << "./RESULTS/Nmax=" << Nmax << "_resampleSize=" << resampling_number << "_k=" << initial_copyNumber << "_s=" << selection_coeff << "_C=" << ecDNA_size_multiplier_factor << "_p=" << fusion_probability << "_q=" << fission_probability << "_Asequence=";
		for (int i = 0; i < ab_pairs.size(); ++i) f << ab_pairs[i].a << "_";
		f << "Bsequence=";
		for (int i = 0; i < ab_pairs.size(); ++i) f << ab_pairs[i].b << "_";
		f << "seed=" << seed;



		DIR *dir = opendir(f.str().c_str());
		if(!dir)
		{
			stringstream g;
			g.str("");
			g << "mkdir -p " << f.str();
			system(g.str().c_str());
		}

		ofstream tissue_file;
		f << "/A=" << sigmoid_a  << "_B=" << sigmoid_b << "_tissue.csv";
		tissue_file.open(f.str().c_str());


		if (verbose_flag) cout << " " << endl;
		if (verbose_flag) cout << "Created output files..." << endl;






		// Sample small number of cells to repeat growth experiment
		all_cell_indices.assign(Ntot , 0);
		for (int i = 0; i < Ntot; ++i)
		{
			all_cell_indices[i] = i;
		}


		// Randomly select n of these (n = resampling_number)
		resampling_indices.clear();
		sample(all_cell_indices.begin(), all_cell_indices.end(), back_inserter(resampling_indices), resampling_number, generator);


		// Re-populate tissue vector and reset book keeping variables and regrow to Nmax cells
		vector<Cell> tissue_resample(resampling_number);
		Ntot = resampling_number;
		N_ecDNA_hot = 0;
		total_unnormalised_division_rate = 0.0;
		total_unnormalised_death_rate = 0.0;

		for (int i = 0; i < resampling_number; ++i)
		{
			tissue_resample[i] = tissue[resampling_indices[i]];
			
			if (tissue[resampling_indices[i]].ecDNA.size() > 0) N_ecDNA_hot += 1;
			total_unnormalised_division_rate += tissue[resampling_indices[i]].division_rate;
			total_unnormalised_death_rate += tissue[resampling_indices[i]].death_rate;
		}

		for (int i = 0; i < resampling_number; ++i)
		{

			// Add re-sampled cell to new tissue array
			tissue[i] = tissue_resample[i];

			if (tissue_resample[i].ecDNA.size() > 0) 
			{

				vector<Occupancy> allOccupancies;
				for (int j = 0; j < tissue_resample[i].ecDNA.size(); ++j)	// Loop over all ecDNA in cell
				{

					counted = false;

					// Check if ecDNA already counted in occupancy vector
					for (int k = 0; k < allOccupancies.size(); ++k)
					{
						if (allOccupancies[k].ecDNA == tissue_resample[i].ecDNA[j])
						{
							counted = true;
							break;
						}
					}


					// Skip if this ecDNA type has already been counted
					if (counted) continue;	

					ecDNA_count = 0;
					for (int jj = 0; jj < tissue_resample[i].ecDNA.size(); ++jj)
					{
						if (tissue_resample[i].ecDNA[j] == tissue_resample[i].ecDNA[jj]) ecDNA_count += 1;
					}

					allOccupancies.push_back(Occupancy(tissue_resample[i].ecDNA[j] , ecDNA_count));

				}


				// Loop back over populated occupancy vector and write data to file 
				for (int k = 0; k < allOccupancies.size()-1; ++k)
				{
					tissue_file << allOccupancies[k].ecDNA << "," << allOccupancies[k].multiplicity << ";";
				}
				tissue_file << allOccupancies[allOccupancies.size()-1].ecDNA << "," << allOccupancies[allOccupancies.size()-1].multiplicity << endl;

			}
		}
		tissue_file.close();







	}

	if (verbose_flag) cout << " " << endl;












	return 0;
}

















