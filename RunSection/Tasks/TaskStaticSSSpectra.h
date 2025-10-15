/////////////////////////////////////////////////////////////////////////
// TaskStaticSSSpectra (RunSection module)  developed by Irina Anisimova and Luca Gerhards.
// ------------------
//
// Simple quantum yield calculation in Liouville space, derived from the
// properties of the Laplace transformation.
//
// Molecular Spin Dynamics Software - developed by Irina Anisimova and Luca Gerhards.
// (c) 2022 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
#ifndef MOD_RunSection_TaskStaticSSSpectra
#define MOD_RunSection_TaskStaticSSSpectra

#include "BasicTask.h"
#include "SpinAPIDefines.h"
#include "SpinSpace.h"
#include "Utility.h"

namespace RunSection
{
	//Because the declaration of i is long define it as a new variable that is easier to use
	using SystemIterator = std::vector<SpinAPI::system_ptr>::const_iterator;

	class TaskStaticSSSpectra : public BasicTask
	{
	private:
		double timestep;
		double totaltime;
		SpinAPI::ReactionOperatorType reactionOperators;

		void WriteHeader(std::ostream &); // Write header for the output file
		static arma::cx_vec ComputeRhoDot(double t, arma::sp_cx_mat& L, arma::cx_vec& K, arma::cx_vec RhoNaught);
		bool ProjectAndPrintOutputLine(auto &_i, SpinAPI::SpinSpace &_space, arma::cx_vec &_rhovec, double &_printedtime, double _timestep, unsigned int &_n, bool &_cidsp, std::ostream &_data_stream, std::ostream &_log_stream);
		bool ProjectAndPrintOutputLineInf(auto &_i, SpinAPI::SpinSpace &_space, arma::cx_vec &_rhovec, double &_printedtime, double _timestep, bool &_cidsp, std::ostream &_datastream, std::ostream &_logstream);

	protected:
		bool RunLocal() override;
		bool Validate() override;

	public:
		// Constructors / Destructors
		TaskStaticSSSpectra(const MSDParser::ObjectParser &, const RunSection &); // Normal constructor
		~TaskStaticSSSpectra();	                                                  // Destructor
		
		bool GetEigenvectors_H0(SpinAPI::SpinSpace &_space, arma::vec &_eigen_val, arma::sp_cx_mat &_eigen_vec_sp) const;												
	};

}

#endif
