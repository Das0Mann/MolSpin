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

namespace SpinAPI
{
    class Function //handles prefactors like cos and sin but essentially any other mathmatical function
    	{
    	public:
    		typedef void*(*FuncPtr)(void*);

    		enum class ReturnType
    		{
    			i = 0, //int
    			d, //double
    			f, //float
    			cd, //complex double
    			undefined
    		};
    	private:
    		std::string m_FunctionName;
    		FuncPtr m_func; //function pointer
    		ReturnType m_funcType; //function type
    		std::string m_variable; //TODO: turn into vector
    		double m_factor;

    	public:
    		arma::cx_double operator()(void* value);
    		Function(FuncPtr, ReturnType, std::string, std::string var = "x", double factor = 1.0);
    		std::string GetVariable();
    		std::string GetName();
    		//using FunctionPtr = std::shared_ptr<Function>;
    	};

    std::shared_ptr<Function> FunctionParser(std::string&, std::string&);

    namespace MathematicalFunctions
	{
		void* sin(void*); //double
		void* cos(void*); //double
		void* scaler(void* = nullptr); //double
	}

}

#endif