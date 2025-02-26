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

			enum class VarType
			{
				d = 0, //double
				z = 1, //complex
				f = 2 //function
			};

			struct VariableDefinition
			{
				std::string name;
				VarType type; 
				
				//only useful for cd and f types
				std::string VariableString = "";
				std::vector<std::string> InternalVariables = {};

				//only useful for f types
				std::shared_ptr<Function> InternalFunction = nullptr;
			};
			
			//depricated 
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
			std::vector<VariableDefinition> m_variables;
			std::vector<std::string> m_duplicates;
			std::string m_FunctionString;
			std::string m_PostFixEquation;

			//depricated
    		std::vector<double> m_factors; 
			std::vector<InternalOperations> m_op;
			std::vector<int> m_VarDepth;
			//std::vector<std::string> m_variables;
		
		private:
			std::complex<double> EvaluateFuncValue(std::vector<void*>);
    	public:
    		arma::cx_double operator()(void* value); //function evaluation for one variable
			arma::cx_double operator()(std::vector<void*> value); //multiple variables
			Function(FuncPtr, ReturnType, std::string, std::vector<VariableDefinition>);
			Function(FuncPtr, ReturnType, std::string, VariableDefinition var = {"", VarType::d, "", {""}, nullptr});
    		
			std::vector<std::string> GetVariable();
    		std::string GetName();

			VariableDefinition FindVar(std::string);

			void SetFunctionString(std::string FuncString)
			{
				m_FunctionString = FuncString;
			}
			void SetPostFix(std::string PostFix)
			{
				m_PostFixEquation = PostFix;
			}
			std::string GetFunctionString()
			{
				return m_FunctionString;
			}
			std::string GetPostFixEq()
			{
				return m_PostFixEquation;
			}
    		//using FunctionPtr = std::shared_ptr<Function>;
    	};

    std::shared_ptr<Function> FunctionParser(std::string&, std::string&, int FuncStartNum = 0, bool ResetComplexNum = false); //function name i.e sin and function contents i.e 0.5x

    namespace MathematicalFunctions
	{
		void* sin(void*); //double
		void* cos(void*); //double
		void* scalar(void* = nullptr); //double
		void* scalarcx(void* = nullptr); //complex double
		void* exp(void*); //double
		void* expcd(void*); //complex double
	}

}

#endif