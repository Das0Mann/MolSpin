/////////////////////////////////////////////////////////////////////////
// Function class (SpinAPI Module)
// ------------------
// The Function class represents a mathmatical function.
//
// Molecular Spin Dynamics Software - developed by Claus Nielsen and Luca Gerhards.
// (c) 2019 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
#ifndef MOD_SpinAPI_Function
#define MOD_SpinAPI_Function

#include<string>
#include<iostream>

namespace SpinAPI
{
    class Function //handles prefactors like cos and sin but essentially any other mathmatical function
    	{
    	public:
    		typedef void*(*FuncPtr)(void*);

    		enum class ReturnType
    		{
    			d = 0, //double
    			cd, //complex double
    			undefined
    		};

			enum class InternalOperations
			{
				p = 0, //plus
				mi, //minus
				mu, //multiply
				d //divide
			};
    	private:
    		std::string m_FunctionName;
    		FuncPtr m_func; //function pointer
    		ReturnType m_funcType; //function type
    		std::vector<std::string> m_variables;
			std::vector<std::string> m_duplicates;
    		std::vector<double> m_factors;
			std::vector<InternalOperations> m_op;
			std::vector<int> m_VarDepth;
			std::string m_FunctionString;
		
		private:
			double EvaluateFuncValue(std::vector<void*>);
    	public:
    		arma::cx_double operator()(void* value); //function evaluation for one variable
			arma::cx_double operator()(std::vector<void*> value); //multiple variables
    		Function(FuncPtr, ReturnType, std::string, std::string var = "x", double factor = 1.0); //constructor for one variable
			Function(FuncPtr, ReturnType, std::string, std::vector<std::string> vars, std::vector<double> factors); //multiple variables
    		
			std::vector<std::string> GetVariable();
    		std::string GetName();

			void SetOp(std::vector<InternalOperations> op) { m_op = op;} //set operations list 
			void SetVarDepth(std::vector<int> vd) {m_VarDepth = vd; } //set bracket depth
			
			void SetFunctionString(std::string FuncString)
			{
				m_FunctionString = FuncString;
				std::cout << FuncString << std::endl;
			}
			std::string GetFunctionString()
			{
				std::cout << m_FunctionString << std::endl;
				return m_FunctionString;
			}
    		//using FunctionPtr = std::shared_ptr<Function>;
    	};

    std::shared_ptr<Function> FunctionParser(std::string&, std::string&);

    namespace MathematicalFunctions
	{
		void* sin(void*); //double
		void* cos(void*); //double
		void* scalar(void* = nullptr); //double
	}

}

#endif