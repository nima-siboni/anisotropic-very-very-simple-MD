/*
 * Compile with:
 *  g++ -o mdsim mdsim.cxx -std=c++0x -O2 -Wall -pedantic
 */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <random>
#include <iomanip>
#include <map>
using namespace std;
extern const double PI2=3.14159265359;

//?? dont forget to exclude forces between walls

double logistic_map()
{
  static double x = 0.5; // start value
  double const mu = 3.99;
  x = mu*x*(1 - x);
  return x;
}

/**
 * this struct defines a coordinate in 2D
 * and will be used for the position, velocity
 * and force of each particle
 */

struct coordinate_type
{
  double x;
  double y;
  coordinate_type() : x(0), y(0) {}
  // copy-constructor
  coordinate_type(coordinate_type const& c)
    : x(c.x)
    , y(c.y) {}
};

/**
 * this struct contains all relevant information and configuration
 * about the MD simulation
 */
struct simulation_data
{
  int N;        // number of particles
  double dt;    // delta t
  double dteq;    // delta t for euilibration of OU process
  int    step;  // current simulation step
  double potential_energy; // per particle
  double kinetic_energy;   // per particle
  double temperature;
  double pressure;
  double virial_xy;
  double virial_yx;
  double virial_xx;
  double virial_yy;
  double virial_pair_xy;
  double virial_pair_xx;
  double virial_pair_yy;
  double virial_pair_yx;
  coordinate_type total_momentum;
  coordinate_type box_length;
  vector<coordinate_type> position;
  vector<coordinate_type> velocity;
  vector<coordinate_type> force;
  vector<double> angle; // the angle of each particle
  vector<double> omega; // the anglular velocity of each particle
  vector<double> torque; // the torque on each particle
  struct potential_data {
    double rr_c;        // square of cutoff of the potential cutoff
    double shift_epot;  // shift of the potential
  } potential;
  double thermostat_temp;

  /* wall related */
  vector<int> type; /* if fluid it is type 0, if upper wall or lower wall it is type 1 or 2 (which is which should be checked*/
  int wall_flag;
  double wall_thickness; /*thickness of the wall, it should be larger than the cutoff*/
  coordinate_type fwall_up;
  coordinate_type fwall_down;
  double shear_rate;
  double equilibration_before_wall_formation; // time for whic the system is equilibrated before the formation of walls
  int nr_liq; // number of particles in the liquid (all particles - wall particles)

};




/**
 * This function is used to calculate the interaction between two particles.
 *
 * Input parameter is the square of the distance between the interacting
 * particles: rr = dx*dx + dy*dy
 * Output is the normalized force f/r and the potential energy
 */
void pair_potential(double rr, double& fval, double& epot, simulation_data const& sim)
{
  // try to calculate the exponents of r in an effective manner:
  // reuse the expontentials of 1/r^2 for both 1/r^{12} and 1/r^{6}
  double const epsilon = 1;
  double rri = 1 / rr;          // 1/r^2
  double r6i = rri * rri * rri; // 1/r^6
  double eps_r6i = epsilon * r6i;
  fval = 48 * rri * eps_r6i * (r6i - 0.5);
  epot = 4 * eps_r6i * (r6i - 1) - sim.potential.shift_epot;
}

/**
 * Equates two simulation_data
 */
void equality(simulation_data& sim, simulation_data& sim0)
{

  for (int i = 0; i < sim.N; ++i) {
    coordinate_type& r1 = sim.position[i];
    coordinate_type& r0 = sim0.position[i];  
    r1.x=r0.x;
    r1.y=r0.y;
  }

}


/**
 * Iterate over all particle pairs and use calculate the force
 * and potential energy of the pair interaction
 */
void calculate_force_and_torque(simulation_data& sim)
{
  // reset the potential energy, as it will be recalcuated in this step
  double epot = 0;
  // reset the force and torque array
  sim.virial_pair_xx =  sim.virial_pair_yy =  sim.virial_pair_yx =  sim.virial_pair_xy = 0;
  sim.virial_xx =  sim.virial_yy =  sim.virial_yx =  sim.virial_xy = 0;

  for (int i = 0; i < sim.N; ++i) {
    sim.force[i].x = 0;
    sim.force[i].y = 0;
    sim.torque[i] = 0;
  }

  for (int i = 0; i < sim.N; ++i) {
    coordinate_type& x1 = sim.position[i];
    //double& theta_i=sim.angle[i]; //theta_i
    // calculate the interaction beween two particles NOT only once
    // (ie. NOT using newtons third law)
    for(int j = 0; j < sim.N; ++j) {
      double pot, fval;
      coordinate_type& x2 = sim.position[j];
      //double& theta_j=sim.angle[j]; //theta_j
      // distance between two particles
      coordinate_type dx;
      dx.x = x1.x - x2.x;
      dx.y = x1.y - x2.y;
      
      // This is gonna exclude the force from particles of one wall on the other one
      if ((sim.type[i]*sim.type[j])!=0) 	continue;
      
      // apply minimum image convention
      if(dx.x > sim.box_length.x/2) {
      	dx.x -= sim.box_length.x;
      }
      if(dx.x < -sim.box_length.x/2) {
       	dx.x += sim.box_length.x;
      }
      if(dx.y > sim.box_length.y/2) {
	dx.y -= sim.box_length.y;
      }
      if(dx.y < -sim.box_length.y/2) {
	dx.y += sim.box_length.y;
      }

      // calculate r^2
      double rr = dx.x*dx.x + dx.y*dx.y;

      // if we are over the cutoff, there is no force & torque contribution
      if(rr > sim.potential.rr_c || i==j)
	continue;

      // calculate pair interaction
      pair_potential(rr, fval, pot, sim);
      //pot is here the potential energy between the pair ij
      epot += pot;

      // convert F/r into vectorial value
      coordinate_type f;
      // f is here the force between pair ij
      f.x = fval*dx.x;
      f.y = fval*dx.y;

      // set the force value 
      sim.force[i].x += f.x;
      sim.force[i].y += f.y;
      sim.torque[i]+=0;

      // Here we calculate the different components of the virial expression for EACH PAIR;
      // i.e. we are calculating f_{ij}.r_{ij}
      // here I am including only those pair which are liquid
      if (j>i && sim.type[i]==0 && sim.type[j]==0){
	sim.virial_pair_xx += f.x*dx.x;
	sim.virial_pair_xy += f.x*dx.y;
	sim.virial_pair_yx += f.y*dx.x;
	sim.virial_pair_yy += f.y*dx.y;
      }

      
    }// end of j loop ; changing j to find the force and torque on i-th particle
    
  }// end of the loop over all particles

  /* Calculating virial by F_i.r_i; one should note that
     as soon as there are wall particles this expression 
     is not anymore equal to the pair_wise calculated 
     expression F_{ij}.r_{ij}: althought the sum is 
     over liq paritcles, in F_i there are forces
     of the wall particles as well. */
  
  for(int i = 0; i < sim.N; ++i) {
    coordinate_type const &f=sim.force[i];
    coordinate_type &r = sim.position[i];
    int &type=sim.type[i];
    
    if (type==0){
      sim.virial_xx += f.x*r.x;
      sim.virial_xy += f.x*r.y;
      sim.virial_yx += f.y*r.x;
      sim.virial_yy += f.y*r.y;
    }
  }
  
  // we want the potential energy per particle
  sim.potential_energy = epot/sim.N;//??

   double area = 0;
   if (sim.wall_flag==0 || sim.wall_flag ==1)
   area= sim.box_length.x * (sim.box_length.y-sim.wall_flag*2*sim.wall_thickness);
   else {cout<<" the wall_flag is supposed to be 0 pr 1"<<endl;}
   sim.virial_xx /= area;
   sim.virial_yy /= area;
   sim.virial_xy /= area;
   sim.virial_yx /= area;
   sim.virial_pair_xx /= area;
   sim.virial_pair_yy /= area;
   sim.virial_pair_xy /= area;
   sim.virial_pair_yx /= area;
  
}



void read_phase_space(simulation_data& sim, string fname)
{
  ifstream input(fname);
  for (int i = 0; i < sim.N; ++i) {
    input >> sim.position[i].x >> sim.position[i].y;
  }
}

double rand_normal(double mean, double stddev)
{//Box muller method
  static double n2 = 0.0;
  static int n2_cached = 0;
  if (!n2_cached)
    {
      double x, y, r;
      do
        {
	  x = 2.0*rand()/RAND_MAX - 1;
	  y = 2.0*rand()/RAND_MAX - 1;

	  r = x*x + y*y;
        }
      while (r == 0.0 || r > 1.0);
      {
	double d = sqrt(-2.0*log(r)/r);
	double n1 = x*d;
	n2 = y*d;
	double result = n1*stddev + mean;
	n2_cached = 1;
	return result;
      }
    }
  else
    {
      n2_cached = 0;
      return n2*stddev + mean;
    }
}


void md_step(simulation_data& sim)
{
  // this function assumes that the force has been calculated
  // the step before
  
  double const dt = sim.dt;
  double temp_old = 0;
  
  double inv_I = 1/1.0; //angular moment I as in torque = I * d2 theta/dt2
  double wall_v=sim.shear_rate*(sim.box_length.y/2.0-sim.wall_thickness);

  
  // first half-step of the verlet integration
  // update the positions and the velocities
  for (int i = 0; i < sim.N; ++i) {
    coordinate_type& pos = sim.position[i];
    coordinate_type& vel = sim.velocity[i];
    double& theta = sim.angle[i];
    double& angvel = sim.omega[i];
    coordinate_type const& force = sim.force[i];
    double const& torquei = sim.torque[i];

    int& particle_type = sim.type[i];

    if (particle_type==0){ //if particles are liquid particles
      pos.x += vel.x*dt + force.x*dt*dt*0.5;
      pos.y += vel.y*dt + force.y*dt*dt*0.5;

      vel.x += force.x*dt*0.5;
      vel.y += force.y*dt*0.5;

      theta += angvel*dt + torquei*inv_I*dt*dt*0.5;
      angvel += torquei*inv_I*dt*0.5;
    }

    if (particle_type==2){ //be caureful the order of v and x are changed with respect to velocity Verlet
      //if(force.x!=0) cout<<"force on particle i="<<i<<" in timestep "<< sim.step<< " equals "<<force.x<<endl;
      vel.x = wall_v;
      vel.y = 0;

      pos.x += vel.x*dt; //+ force.x*dt*dt*0.5;
      pos.y += 0.0; //vel.y*dt + force.y*dt*dt*0.5;
      
    }

    if (particle_type==1){ //be caureful the order of v and x are changed with respect to velocity Verlet
      vel.x = -1.0*wall_v;
      vel.y = 0;

      pos.x += vel.x*dt; //+ force.x*dt*dt*0.5;
      pos.y += 0.0; //vel.y*dt + force.y*dt*dt*0.5;

    }



    // apply periodic boundary conditions for positions
    if(pos.x > sim.box_length.x) {
        pos.x -= sim.box_length.x;
    }
    if(pos.x < 0) {
        pos.x += sim.box_length.x;
    }
    if(pos.y > sim.box_length.y) {
        pos.y -= sim.box_length.y;
    }
    if(pos.y < 0) {
        pos.y += sim.box_length.y;
    }

    //?? applying periodic boundary condition for angles
    
    temp_old += vel.x*vel.x + vel.y*vel.y;
  } // End of the first half step

  // update the forces as the positions and angles have changed
  // Be careful of bonds and such in periodic boundary conditions
  calculate_force_and_torque(sim);

  // second half-step of the verlet integration
  // updates the velocities and angular velocities
  for (int i = 0; i < sim.N; ++i) {
    coordinate_type& vel = sim.velocity[i];
    coordinate_type const& force = sim.force[i];
    double& angvel = sim.omega[i];
    double const& torquei = sim.torque[i];
    int& particle_type = sim.type[i];
    if (particle_type==0){
      vel.x += 0.5*dt*force.x;
      vel.y += 0.5*dt*force.y;
      angvel += torquei*inv_I*dt*0.5;
    }

  }
    
}


// write the positions and velocities into a single file
void write_phase_space(simulation_data const& sim)
{
  //ofstream output("phase_space" + to_string(sim.step) +".lammpstrj");
  if (sim.step==0)
    {
      ofstream output("phase_space.lammpstrj");
    }
  ofstream output("phase_space.lammpstrj",std::ofstream::app);
  // the following lines set the output precision very high
  output.precision(15);
  output.setf(ios::scientific, ios::floatfield);
  output<<"ITEM: TIMESTEP"<<endl;
  output<<sim.step<<endl;
  output<<"ITEM: NUMBER OF ATOMS"<<endl;
  output<<sim.N<<endl;
  output<<"ITEM: BOX BOUNDS xx yy zz"<<endl;
  output<<"0.0 "<<sim.box_length.x<<endl;
  output<<"0.0 "<<sim.box_length.y<<endl;
  output<<"0.0 1.0"<<endl;
  output<<"ITEM: ATOMS id x y z mux muy muz type"<<endl;
  // write the phase space data
  int liquid_or_wall;
  for (int i = 0; i < sim.N; ++i) {

    if (sim.type[i]==0) liquid_or_wall=1;
    if (sim.type[i]!=0) liquid_or_wall=2;
    output << i+1 << '\t' << sim.position[i].x << '\t' << sim.position[i].y <<'\t'<< " 0 " << '\t' << cos(sim.angle[i]) << '\t' << sin(sim.angle[i]) << '\t' << " 0 "<<'\t'<<liquid_or_wall<<endl;
  }
}

/**
 * Calculate the kinetic_energy and thermostating it with all positional degrees of freedom without removing CM; so be careful when using with shear
 */

void thermostat(simulation_data& sim)
{
  double ekin_yy = 0;
  double ekin_xx = 0;
  for (int i = 0; i < sim.N; ++i) {
    coordinate_type const& vel = sim.velocity[i];
    // assume unit mass
    if (sim.type[i]==0) {
      ekin_yy +=  vel.y*vel.y;
      ekin_xx +=  vel.x*vel.x;
    }
  }
  // missing factor of two, and divide by the number of particles
  ekin_xx /= 1.*sim.nr_liq;
  ekin_yy /= 1.*sim.nr_liq;
  sim.kinetic_energy = 0.5*(ekin_xx+ekin_yy);

  if (sim.shear_rate==0){
    sim.temperature = 0.5*(ekin_xx+ekin_yy);
  }
  else
    {
      sim.temperature = ekin_yy;
    }
  

  double scale=sqrt(sim.thermostat_temp/sim.temperature);
  for (int i = 0; i < sim.N; ++i) {
    coordinate_type & vel = sim.velocity[i];
    // assume unit mass
    if (sim.type[i]==0) {
      if (sim.shear_rate==0) {
	vel.x = vel.x*scale;
      }
      vel.y = vel.y*scale;
    }
  }
}
 




/**
 * the following things are being done in this function
 * - place all particles on a square lattice
 * - generate the velocities using the logistic map
 * - shift the centre of mass to to zero
 * - scale the particles velocities to temperature T=1
 */
void set_initial_phase_space(simulation_data& sim)
{
  // lattice constant
  double a = sqrt(sim.box_length.x * sim.box_length.y/sim.N);
  // place particles on lattice
  int Nx = (int) sqrt(sim.N);
  int Ny = (int) sqrt(sim.N);
  for(int i = 0; i < Nx; ++i) {
    for(int j = 0; j < Ny; ++j) {
      sim.position[i + Nx*j].x = a*(i+0.5);
      sim.position[i + Nx*j].y = a*(j+0.5);
    }
  }
  double pi = 3.141592;
  for (int i = 0; i < sim.N; ++i) {
    double& ang = sim.angle[i];
    ang = (logistic_map() - 0.5)*pi;
  }
  sim.total_momentum = coordinate_type(); // reset to zero
  double sum_angular_vels = 0;
  for (int i = 0; i < sim.N; ++i) {
    coordinate_type& v = sim.velocity[i];
    double& angvel = sim.omega[i];
    v.x = logistic_map() - 0.5;
    v.y = logistic_map() - 0.5;
    angvel = logistic_map() - 0.5;
    sum_angular_vels += angvel;
    sim.total_momentum.x += v.x;
    sim.total_momentum.y += v.y;
  }
  // shift total momentum to zero
  coordinate_type shift;
  shift.x = sim.total_momentum.x/sim.N;
  shift.y = sim.total_momentum.y/sim.N;
  double vv = 0;
  for (int i = 0; i < sim.N; ++i) {
    coordinate_type& v = sim.velocity[i];
    v.x -= shift.x;
    v.y -= shift.y;
    // calculate new v^2
    vv += v.x*v.x + v.y*v.y;
  }

  double shift_ang_vel = sum_angular_vels/sim.N;
  for (int i = 0; i < sim.N; ++i) {
    double& angvel = sim.omega[i];
    angvel -= shift_ang_vel;
  }
  
  // scale the velocities using the scaling factor
  double scale_factor = sqrt(2*sim.N*sim.temperature/vv);
  for (int i = 0; i < sim.N; ++i) {
    sim.velocity[i].x *= scale_factor;
    sim.velocity[i].y *= scale_factor;
  }

  //?? angular velocities should be also rescaled.


  for (int i = 0; i < sim.N; ++i) {
    sim.type[i]=0; /* all particles are by default liquid particles*/
  }
  sim.nr_liq = sim.N;

}

/**
 * calculate the msd of the sim
 */

double calculate_msd(simulation_data& sim, simulation_data& sim0)
{
  double dx, dy;
  double result = 0;
  for (int i = 0; i < sim.N; ++i) {
    coordinate_type& r1 = sim.position[i];
    coordinate_type& r0 = sim0.position[i];  
    dx = r1.x -r0.x;
    dy = r1.y -r0.y;
    //cout<<dx<< " " <<dy <<endl;
    result += (dx*dx + dy*dy)/sim.N;
  }
  return result;
}


/**
 * calculate the msd of the sim
 */
void print_displacements(simulation_data& sim, simulation_data& sim0,int time)
{
  double dx, dy;
  string str1("displacements.");
  str1.append(to_string(time));
  ofstream displacementfile(str1);
  
  for (int i = 0; i < sim.N; ++i) {
    coordinate_type& r1 = sim.position[i];
    coordinate_type& r0 = sim0.position[i];  
    dx = r1.x -r0.x;
    dy = r1.y -r0.y;
    displacementfile<<dx<<" "<<dy<<endl;
  }
}


/**
 * initializes the simulation (box size, number of
 * particles, potential) and calls set_initial_phase_space()
 * in order to place the particles on the lattice and set
 * the initial velocity according to sim.temperature
 */
void initialize(simulation_data& sim)
{
  sim.step = 0;
  //  sim.N = sim.N;
  sim.temperature = 1;
  // set up arrays
  sim.position.resize(sim.N);
  sim.velocity.resize(sim.N);
  sim.force.resize(sim.N);
  sim.angle.resize(sim.N); // the angle of each particle
  sim.omega.resize(sim.N); // the anglular velocity of each particle
  sim.torque.resize(sim.N); // the torque on each particle
  sim.type.resize(sim.N);
  sim.wall_flag = 0;
  if(sim.equilibration_before_wall_formation==0) sim.equilibration_before_wall_formation=1; //default value
 
  // sim.box_length.x = L;
  // sim.box_length.y = L;
  // intitialize the potential
  sim.potential.shift_epot = 0;
  // calculate the potential energy shift
  double r_c = pow(2, 1./6.);
  sim.potential.rr_c = r_c*r_c;
  double fval, epot_shift;
  pair_potential(r_c*r_c, fval, epot_shift, sim);
  sim.potential.shift_epot = epot_shift;

  // set the starting positions
  set_initial_phase_space(sim);
}


void define_walls(simulation_data& sim)
{
  // defining walls: the upper wall is 2, lower wall is 1, liquid between them is 0
  
  sim.wall_flag=1;
  int upperwall=0;
  int lowerwall=0;
  for (int i = 0; i < sim.N; ++i) {
    sim.type[i]=0; /*that is the default which is liquid*/
  }
  for (int i = 0; i < sim.N; ++i) {
    coordinate_type& pos = sim.position[i];
    if (pos.y>(sim.box_length.y-sim.wall_thickness)){
	sim.type[i]=2;
	upperwall++;
      }
    if (pos.y<sim.wall_thickness) {
      sim.type[i]=1;
      lowerwall++;
    }

  }
  cout<<" number of particles in the upper wall: "<<upperwall<<endl;
  cout<<" number of particles in the lower wall: "<<lowerwall<<endl;
  sim.nr_liq = sim.N - (upperwall+lowerwall);
}


void output_thermodynamic_variables(simulation_data& sim){
  ofstream thermodynamics;
  if (sim.step==0){
    thermodynamics.open("thermodynamics");
    thermodynamics<<"#1 time"
		  << "\n" << "#2- temp" 
		  << "\n" << "#3- kin energy"
		  << "\n" << "#4- potentialenergy??"
		  << "\n" << "#5- pressure"
		  << "\n" << "#6- virial_xx"
		  << "\n" << "#7- virial_xy"
		  << "\n" << "#8- virial_yx"
		  << "\n" << "#9- virial_yy"
		  << "\n" << "#10- virial_pair_xx"
		  << "\n" << "#11- virial_pair_xy"
		  << "\n" << "#12- virial_pair_yx"
		  << "\n" << "#13- virial_pair_yy"
		  << "\n" << "#14- total_momentum.x"
		  << "\n" << "#15- total_momentum.y"
		  << "\n" << "#16- fwall_up.x"
		  << "\n" << "#17- fwall_up.y"
		  << "\n" << "#18- fwall_down.x"
		  << "\n" << "#19- fwall_down.y"
		  << endl;
  }
  else{
    thermodynamics.open("thermodynamics",std::ofstream::app);
  }
  
  double ekin_xx = 0;
  double ekin_yy = 0;
  sim.total_momentum = coordinate_type();
  for (int i = 0; i < sim.N; ++i) {

    if (sim.type[i]==0) {
      coordinate_type const& vel = sim.velocity[i];
      sim.total_momentum.x += vel.x;
      sim.total_momentum.y += vel.y;
      ekin_yy +=  vel.y*vel.y;
      ekin_xx +=  vel.x*vel.x;
    }
  }
  // missing factor of two, and divide by the number of particles
  ekin_xx /= 1.*sim.nr_liq;
  ekin_yy /= 1.*sim.nr_liq;
  sim.kinetic_energy = 0.5*(ekin_xx+ekin_yy);

  if (sim.shear_rate==0){
    sim.temperature = 0.5*(ekin_xx+ekin_yy);
  }
  else
    {
      sim.temperature = ekin_yy;
    }
  
  double area=0;
  if (sim.wall_flag==0 || sim.wall_flag ==1)
    area= sim.box_length.x * (sim.box_length.y-sim.wall_flag*2*sim.wall_thickness);
  else {cout<<" the wall_flag is supposed to be 0 pr 1"<<endl;}
    
  sim.pressure = sim.temperature*sim.nr_liq/area + 0.5*(sim.virial_xx+sim.virial_yy);
  


  /*wall force calculation */
  sim.fwall_up.x = sim.fwall_up.y = 0;
  sim.fwall_down.x = sim.fwall_down.y = 0;

  for (int i = 0; i < sim.N; ++i) {
    int& type=sim.type[i];
    if (type==2) {
      sim.fwall_up.x += sim.force[i].x;
      sim.fwall_up.y += sim.force[i].y;
    }
    if (type==1) {
      sim.fwall_down.x += sim.force[i].x;
      sim.fwall_down.y += sim.force[i].y;
    }
  }

  sim.fwall_up.x /= sim.box_length.x;
  sim.fwall_up.y /= sim.box_length.x;
  sim.fwall_down.x /= sim.box_length.x;
  sim.fwall_down.y /= sim.box_length.x;

  
  thermodynamics << sim.step*sim.dt
		 << " " << sim.temperature 
		 << " " << sim.kinetic_energy 
		 << " " << sim.potential_energy
		 << " " << sim.pressure
		 << " " << sim.virial_xx
		 << " " << sim.virial_xy
		 << " " << sim.virial_yx
		 << " " << sim.virial_yy
    		 << " " << sim.virial_pair_xx
		 << " " << sim.virial_pair_xy
		 << " " << sim.virial_pair_yx
		 << " " << sim.virial_pair_yy
    		 << " " << sim.total_momentum.x
		 << " " << sim.total_momentum.y
		 << " " << sim.fwall_up.x
		 << " " << sim.fwall_up.y
    		 << " " << sim.fwall_down.x
    		 << " " << sim.fwall_down.y
		 << endl;
      
}


int main(int argc, char **argv)
{
  double Tmax;
  simulation_data sim;
  simulation_data sim0;
  if (argc!=10) {
    cout<<"Please enter Nparticle, system_size, dt,Tmax, temperature, wall thickness, shear rate, equilibration before wall formation, randomseed"<<endl;
    Tmax=0;
  }else{
    sim.N = atoi(argv[1]);
    cout<<"# N is "<<sim.N<<endl;
    sim.box_length.x = atof(argv[2]);
    sim.box_length.y = atof(argv[2]);
    cout<<"# Lx and Ly are "<<sim.box_length.x<<" "<<sim.box_length.y<<endl;
    sim.dt=atof(argv[3]); 
    cout<<"# dt is "<<sim.dt<<endl;
    Tmax=atof(argv[4]);
    cout<<"# Tmax is "<<Tmax<<endl;

    sim.thermostat_temp = atof(argv[5]);
    cout<<"# Tbath "<<sim.thermostat_temp<<endl;

    sim.wall_thickness=atof(argv[6]);
    cout<<"# wall thickness is "<<sim.wall_thickness<<endl;

    sim.shear_rate=atof(argv[7]);
    cout<<"# shear rate is "<<sim.shear_rate<<endl;

    sim.equilibration_before_wall_formation=atof(argv[8]);
    cout<<"# Equilibration before wall formation "<<sim.equilibration_before_wall_formation<<endl;

    srand(atoi(argv[9])+1);
    
    
    initialize(sim);
    
    sim0.N = atoi(argv[1]);
    sim0.box_length.x = atof(argv[2]);
    sim0.box_length.y = atof(argv[2]);
    
    initialize(sim0);

  }
  
  equality(sim0,sim);
  write_phase_space(sim);

  int steps = Tmax/sim.dt;
  int thermostating_interval=(int)(.2/sim.dt);
  int dump_interval=(int)(0.1/sim.dt);
  cout << "Simulating " << steps << " steps.." << endl;

  // calculate the inital force
  calculate_force_and_torque(sim);

  // print at max 1000 theromodynamic values
  int next_print_interval = 0; //max(steps/10, 1);

  for (int i = 0; i < steps; ++i) {
    //if (i%(int)(1/sim.dt)) cout<<" "<<i*sim.dt;
    if ((i%thermostating_interval)==0) {
      thermostat(sim);
    }
    if (i==(int) (sim.equilibration_before_wall_formation/sim.dt)){
	define_walls(sim); cout<<"# walls are formed at timestep "<<i<<endl;
    }
    if (i%dump_interval==0) write_phase_space(sim);

    if ((i == next_print_interval) || i == 0) {
      //print_displacements(sim,sim0,i);
      output_thermodynamic_variables(sim);
      next_print_interval+=(int)(0.01/sim.dt);
    }
    
    sim.step++;
    md_step(sim);
  }
  // write out the final configuration
  write_phase_space(sim);
  // the output can be plotted in gnuplot with:
  // set size square
  // 'phase_space_1000000' u 1:2:(0.1) with circle linecolor rgb '#000000' fillstyle solid title '10^6 steps'
  return 0;
}
