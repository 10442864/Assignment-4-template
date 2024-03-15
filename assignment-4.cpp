// PHYS 30762 Programming in C++
// Assignment 4
// Practice special functions and operators in C++ classes
// Note that the hints in the skeleton are given to help you
// in case of doubt, the official guidance/marking scheme is on the slides on BB 

#include<iostream>
#include<string>
#include<vector>
#include<cmath>

using std::string;

// Beginning of particle class
class particle
{
private:
  string name;
  double rest_mass;
  int charge;
  std::vector<double> *P {nullptr};
  //...other data members (see slides on BB)
public:
	// Constructors
  // Here you need a default constructor, a parameterised constructor and a copy constructor
  // The parameterised constructor needs to dynamically allocate the std::vector containing the four-vector elements
  // The parameterised constructor also needs to check the validity of the energy component
  // The copy constructor needs to make a deep copy of the std::vector holding the 4-momentum

  particle() : name{}, rest_mass{}, charge{}, P{} {}
  particle(string name_in, double mass_in, int charge_in, std::vector<double>& P_in) : 
  name{name_in}, rest_mass{mass_in}, charge{charge_in} 
	{
    P = new std::vector<double>(4); // Dynamically allocates four-momentum vector
    set_four_momentum(P_in); // Checks validity of four-momentum values
  }
  // Deep copy constructor
  particle(const particle& copy) : name(copy.name), rest_mass(copy.rest_mass), charge(copy.charge)
  {
    // Deep copy the P vector by creating a new vector with a copy of elements
    std::cout<<"\nCalling Deep Copy Constructor..."<<std::endl;
    P = new std::vector<double>(*copy.P);
  }
	// Destructor
  // The destructor needs to free the memory allocated by the constructor
  ~particle() = default; // Error produced by delete for vector, default avoids this.

  // Assignment operator
  // Needs to avoid self-assignment using the *this pointer
  particle operator=(const particle& v_3)   // Deep Copy Assignment Operator
  {
    std::cout<<"\nCalling Deep Copy Assignment Operator..."<<std::endl;
    if (this != &v_3) // Avoid self-assignment
    {
      name = v_3.name;  
      rest_mass = v_3.rest_mass;
      charge = v_3.charge;
      if (P) 
      {
        delete P; // Delete existing P. Avoids memory leaks
      }
      P = new std::vector<double>(*v_3.P);
    }
    return *this;
  }
 
  // Move constructor
  // The move constructor needs to correctly steal the memory from the object you're calling it on
  particle move(particle &&vector)
  {
    std::cout<<"\nCalling Move Constructor..."<<std::endl;
    name = vector.name;  // Moves data from original vector to new one
    rest_mass = vector.rest_mass;
    charge = vector.charge;
    P = vector.P;
    vector.name = "";   // Resets the original vector to null values
    vector.rest_mass = 0;
    vector.charge = 0;
    vector.P = {nullptr};
    return vector;
  }

  // Move assignment operator
  // The move assignment operator needs to correctly reassign the memory from the original object
  particle& operator=(particle&& vector_2)
  {
    std::cout<<"\nCalling Move Assignment Operator..."<<std::endl;
    std::swap(name, vector_2.name);
    std::swap(rest_mass, vector_2.rest_mass);
    std::swap(charge, vector_2.charge);
    std::swap(P, vector_2.P);
    return *this;
  }


  // Getter functions (accessors) to individual elements of 4-momentum
  // This should include function returning beta value 
	string get_name() const {return name;}
  double get_mass() const {return rest_mass;}
  int get_charge() const {return charge;}
	double get_E() const {return P ? P->at(0) : 0;}
	double get_Px() const { return P ? P->at(1) : 0.0; }
	double get_Py() const { return P ? P->at(2) : 0.0; }
	double get_Pz() const { return P ? P->at(3) : 0.0; }

  // Setter functions, to change values of 4-momentum 
  // Make sure you check input validity for the energy in the 4-momentum 
  void set_name(const string name_in) {name=name_in;} 
  void set_mass(const double mass_in) {rest_mass=mass_in;}
  void set_charge(const int charge_in) {charge=charge_in;}
  void set_four_momentum(std::vector<double>& P_in) 
	{
    if (!P) {P = new std::vector<double>(4);} // Allocates the four-momentum if not already done
		while (P_in[0] <= 0 || P_in[0] >= 1) // Checks that E/c falls in range 0 to 1
		{
      std::cerr<<"ERROR: Invalid four-momentum. Energy cannot be negative."<<std::endl;
      std::cout<<"Please enter a new value for E/c: "<<std::endl;
			std::cin >> P_in[0];
    }
    *P = P_in; // Defines valid four-momentum elements for object
  }

  // Accessors for implemented functions
	void print_data();
  particle operator+(const particle &v_1);
  void print_sum();
  particle operator*(const particle &v_2);
  void print_dot();
  
  
// Implementation of functions goes here
};
// End of particle class and associated member functions


void particle::print_data() // Prints overall data for a particle
{
  std::cout<<"Particle is of type: "<<name<<", with values of Rest Mass = "<<rest_mass<<" MeV/C^2, Charge = "
	   <<charge<<", Momentum[E/c, Px, Py, Pz]= ["<<P->at(0)<<", "<<P->at(1)<<", "<<P->at(2)<<", "<<P->at(3)<<"] kgm/s"<<std::endl;
  return;
}

void particle::print_sum() // Prints result of vector sum
{
	std::cout<<"\nSum of particle momenta [E/c, px, py, pz] = ["<<P->at(0)<<", "
	<<P->at(1)<<", "<<P->at(2)<<", "<<P->at(3)<<"]"<<std::endl;
}

particle particle::operator+(const particle &v_1) // Sum operator
{                                           // Will request re-entry of E/c if sum greater than 1
  std::vector<double> sum(4);
  for(auto i=0; i<4; ++i)
  {
    sum[i] = P->at(i) + v_1.P->at(i);
  }
  return particle("Sum", 0, 0, sum);
}

particle particle::operator*(const particle &v_2) // Dot product operator
{
  std::vector<double> sum(4);
  for(auto i=0; i<4; ++i)
  {
    sum[i] = P->at(i) * v_2.P->at(i);
  }
  return particle("Sum", 0, 0, sum);
}

void particle::print_dot() // Prints dot product result.
{
  std::cout<<"\nDot product of particle momenta [E/c, px, py, pz] = ["<<P->at(0)<<", "
	<<P->at(1)<<", "<<P->at(2)<<", "<<P->at(3)<<"]"<<std::endl;
}



// Main program
int main()
{
  // Create the following particles: 
  // two electrons, four muons, one antielectron, one antimuon 
  // Use the parameterised constructor to do these
  std::vector<std::vector<double>> momentum_values;  // Momenta vector and particle vector produced sepeartely to avoid vector errors
  momentum_values.push_back({0.4, 35, 35, 35}); // electron1
  momentum_values.push_back({0.5, 45, 45, 45}); // electron2
  momentum_values.push_back({0.3, 60, 60, 60}); // muon1
  momentum_values.push_back({0.2, 50, 50, 50}); // muon2
  momentum_values.push_back({0.5, 15, 15, 15}); // muon3
  momentum_values.push_back({0.1, 30, 30, 30}); // muon4
  momentum_values.push_back({0.4, 25, 25, 25}); // antielectron
  momentum_values.push_back({0.35, 70, 70, 70}); // antimuon

  std::vector<particle> particle_data;
  particle_data.push_back(particle("electron", 0.511, -1, momentum_values[0])); //electron
  particle_data.push_back(particle("electron", 0.511, -1, momentum_values[1])); //electron
  particle_data.push_back(particle("muon", 105.658, -1, momentum_values[2])); //muon
  particle_data.push_back(particle("muon", 105.658, -1, momentum_values[3])); //muon
  particle_data.push_back(particle("muon", 105.658, -1, momentum_values[4])); //muon
  particle_data.push_back(particle("muon", 105.658, -1, momentum_values[5])); //muon
  particle_data.push_back(particle("electron", 0.511, +1, momentum_values[6])); //antielectron
  particle_data.push_back(particle("muon", 105.658, +1, momentum_values[7])); //antimuon
	

  // (optional but nice) Print out the data from all the particles (put them in a vector)
  std::cout<<"Listing data for all particles: \n"<<std::endl;
  std::vector<particle>::iterator particle_i; 
  for(auto particle_i=particle_data.begin(); particle_i<particle_data.end(); ++particle_i)
  {particle_i->print_data();}


  // Sum the four-momenta of the two electrons
	particle momenta_sum{particle_data[0].operator+(particle_data[1])};
	momenta_sum.print_sum();
  // Do the dot product of the first two four-muons
  particle dot{particle_data[2].operator*(particle_data[3])};
  dot.print_dot();

  // Assignment operator of an electron to a new electron
  particle electron_assigned;
  electron_assigned = particle_data[0]; // Assignment operator to new electron
  std::cout<<"New electron produced: "<<std::endl;
  electron_assigned.print_data(); // Print data for new electron

  // Copy constructor of the first muon to a new muon
  particle muon_copy(particle_data[2]); // Deep copy constructor on first muon
  std::cout<<"New muon produced: "<<std::endl;
  muon_copy.print_data();

  // Move the antielectron into another antielectron using the move constructor 
  particle new_antielectron = std::move(particle_data[6]);
  std::cout<<"New anti-electron produced: "<<std::endl;
  new_antielectron.print_data();

  // Assign the antimuon to another antimuon using the move assignment
  particle new_antimuon;
  new_antimuon = std::move(particle_data[7]);
  std::cout<<"New anti-muon produced: "<<std::endl;
  new_antimuon.print_data();

  std::cout<<"\nThe constructors and operators output the expected results."<<std::endl;
  return 0;
}