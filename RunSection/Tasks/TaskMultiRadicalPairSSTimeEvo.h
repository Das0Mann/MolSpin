/////////////////////////////////////////////////////////////////////////
// TaskMultiRadicalPairSSTimeEvo (RunSection module)
// ------------------
//
// Simple time-evolution calculation in Liouville space.
//
// -- Multi-system version: Allows transitions between SpinSystems --
//
// Molecular Spin Dynamics Software - developed by Claus Nielsen and Luca Gerhards.
// (c) 2019 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
#ifndef MOD_RunSection_TaskMultiRadicalPairSSTimeEvo
#define MOD_RunSection_TaskMultiRadicalPairSSTimeEvo

#include "BasicTask.h"
#include "SpinSpace.h"
#include "SpinSystem.h"
#include "SpinAPIDefines.h"

namespace RunSection
{
	struct SubSystemTransition
	{
		SpinAPI::transition_ptr transition;
		int type; //0 = Transition out of one subsystem, 1 = Transition between two subsytems
		std::string source;
		std::string target;
		SubSystemTransition()
			:transition(nullptr), type(-1), source(""), target("")
		{
		}
	};

	class TaskMultiRadicalPairSSTimeEvo : public BasicTask
	{
	private:
		double timestep;
		double totaltime;
		SpinAPI::ReactionOperatorType reactionOperators;

		void WriteHeader(std::ostream &); // Write header for the output file

		// Private method that gathers and outputs the results from a given time-integrated density operator
		void GatherResults(const arma::cx_mat &, const SpinAPI::SpinSystem &, const SpinAPI::SpinSpace &, std::vector<std::complex<double>>& traj);
		
		void StateYield(double, double&, const std::vector<std::complex<double>>&, std::vector<double>& ); //Method that calculates the yeild for a given state
		double simpson_integration(std::vector<double> x_list, std::vector<double> y_list);
		
		bool GenerateHamiltonian(const std::vector<SpinAPI::interaction_ptr> interactions, arma::sp_cx_mat& H, int dimension, std::shared_ptr<SpinAPI::SpinSpace> SpinSystem);
		bool GenerateReactionOperator(const std::vector<SubSystemTransition> transitions, arma::sp_cx_mat& K, int dimension, std::shared_ptr<SpinAPI::SpinSpace> SpinSystem, std::string source);

		bool ValidateSubSystems();

		//RungeKutta4 method 
		bool RungeKutta4(arma::sp_cx_mat& L, arma::cx_vec& RhoNaught, arma::cx_vec& drhodt, double timestep);
		arma::cx_vec ComputeRhoDot(arma::sp_cx_mat& L, arma::cx_vec& K, amra::cx_vec& RhoNaugt); //k in this case isn't refering to the reaction operator K but k_i, i E {1,2,3,4} for the RK4 method

	protected:
		bool RunLocal() override;
		bool Validate() override;

	public:
		// Constructors / Destructors
		TaskMultiRadicalPairSSTimeEvo(const MSDParser::ObjectParser &, const RunSection &); // Normal constructor
		~TaskMultiRadicalPairSSTimeEvo();												   // Destructor
	};
}

#endif
