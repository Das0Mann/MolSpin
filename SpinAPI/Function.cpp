/////////////////////////////////////////////////////////////////////////
// Function implementation (SpinAPI Module)
// ------------------
// The Function class represents a mathmatical function.
//
// Molecular Spin Dynamics Software - developed by Claus Nielsen and Luca Gerhards.
// (c) 2019 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////

#include <cctype>
#include <math.h> 
#include <memory>
#include <armadillo>

#include "Function.h"

namespace SpinAPI
{
    //Function class stuff
	arma::cx_double Function::operator()(void* value)
	{
		arma::cx_double ReturnValue;
		if(m_funcType == ReturnType::i)
		{
			int val = (int)m_factor * *(int*)value; 
			int _val = *(int*)m_func((void*)(int*)&val);
			ReturnValue = arma::cx_double((double)_val, 0);
		}
		else if(m_funcType == ReturnType::d)
		{
			double val = m_factor * *(double*)value; 
			double _val = *(double*)m_func((void*)(double*)&val);
			ReturnValue = arma::cx_double(_val, 0);
		}
		else if(m_funcType == ReturnType::f)
		{
			float val = (float)m_factor * *(float*)value; 
			float _val = *(float*)m_func((void*)(float*)&val);
			ReturnValue = arma::cx_double((float)_val, 0);
		}
		else if(m_funcType == ReturnType::cd)
		{
			arma::cx_double val = (arma::cx_double)m_factor * *(arma::cx_double*)value; 
			arma::cx_double _val = *(arma::cx_double*)m_func((void*)(arma::cx_double*)&val);
			ReturnValue = _val;
		}

		return ReturnValue;
	}

	Function::Function(FuncPtr FunctionPtr, ReturnType type, std::string name, std::string var, double factor)
	{
		m_func = FunctionPtr;
		m_funcType = type;
		m_FunctionName = name;
		m_variable = var;
		m_factor = factor;
	}

	std::string Function::GetVariable()
	{
		return m_variable;
	}

	std::string Function::GetName()
	{
		return m_FunctionName;
	}

	std::shared_ptr<Function> FunctionParser(std::string& func, std::string& var)
	{
		std::string buffer = "";
		bool decimal = false;
		std::pair<double,double> nums = {0.0, 0.0};
		bool value = true;
		for(auto c = var.cbegin(); c != var.cend(); c++)
		{
			if((*c) == '.')
			{
				decimal = true;
				nums.first = std::stod(buffer);
				buffer = "";
			}
			else if(std::isdigit((*c)) == 0 && value)
			{
				value = false;
				double val = std::stod(buffer);
				if(!decimal)
				{
					nums.first = val;
				}
				else
				{
					int DivideBy = buffer.length();
					val = val * std::pow(10.0, -1.0 * (double)DivideBy);
					nums.second = val;
				}
				buffer = (*c);
			}
			else
			{
				buffer += (*c);
			}
		}

		std::string variable = buffer;
		double factor = nums.first + nums.second;

		std::string FunctionName = "";
		std::vector<int> removed = {};
		for(auto c = func.cbegin(); c != func.cend(); c++)
		{
			auto a = func.cend() - func.cbegin();
			if((*c) == '+' || (*c) == '-')
			{
				continue;
			} 
			else
			{
				FunctionName += (*c);
				removed.push_back(c - func.begin());
			}
		}
		int r = 0;
		for(int i = 0; i < removed.size(); i++)
		{
			func.erase(func.begin() + removed[i] - r);
			r++;
		}

		std::shared_ptr<Function> _func = nullptr;
		if(FunctionName.compare("sin") == 0)
		{
			 _func = std::make_shared<Function>(MathematicalFunctions::sin, Function::ReturnType::d, FunctionName, variable, factor);
		}
		if(FunctionName.compare("cos") == 0)
		{
			_func = std::make_shared<Function>(MathematicalFunctions::cos, Function::ReturnType::d, FunctionName, variable, factor);
		}
		
		return _func;
	}

	namespace MathematicalFunctions
	{
		void* sin(void* value) //double
		{
			double val = *(double*)value;
			double _val = std::sin(val);
			double* v = &_val;
			return (void*)v;
		} 

		void* cos(void* value) //double
		{
			double val = *(double*)value;
			double _val = std::cos(val);
			double* v = &_val;
			return (void*)v;
		}

		void* scaler(void* value)
		{
			double _val = *(double*)value ;
			double* v = &_val;
			return (void*)v;
		}
	}
}