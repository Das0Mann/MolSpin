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
#include <stack>
#include <unordered_map>

#include "Function.h"

namespace SpinAPI
{
    //Function class stuff
	arma::cx_double Function::operator()(void* value)
	{
		arma::cx_double ReturnValue;
		if(m_funcType == ReturnType::d)
		{
			double val = m_factors[0] * *(double*)value; 
			double _val = *(double*)m_func((void*)(double*)&val);
			ReturnValue = arma::cx_double(_val, 0);
		}
		else if(m_funcType == ReturnType::cd)
		{
			double val = m_factors[0] * *(double*)value; 
			arma::cx_double _val = *(arma::cx_double*)m_func((void*)(arma::cx_double*)&val);
			ReturnValue = _val;
		}

		return ReturnValue;
	}

	arma::cx_double Function::operator()(std::vector<void*> value)
	{
		arma::cx_double ReturnValue;
		if(m_funcType == ReturnType::d)
		{
			double val = EvaluateFuncValue(value);
			double _val = *(double*)m_func((void*)(double*)&val);
			ReturnValue = arma::cx_double(_val, 0);
		}
		else if(m_funcType == ReturnType::cd)
		{
			double val = EvaluateFuncValue(value);
			arma::cx_double _val = *(arma::cx_double*)m_func((void*)(arma::cx_double*)&val);
			ReturnValue = _val;
		}

		return ReturnValue;
	}

	double Function::EvaluateFuncValue(std::vector<void*> values)
	{
		std::vector<double> VarValues;
		for(int i = 0; i < m_variables.size(); i++)
		{
			if(m_variables[i] == "")
				VarValues.push_back(1.0);
			else
				VarValues.push_back(0.0);
		}
		int index = 0;
		for(auto v = values.begin(); v != values.end(); v++)
		{
			if(VarValues[index] == 0.0)
			{
				VarValues[index] = *(double*)(*v);
			}
			else
			{
				v--;
			}
			index++;

		} 

		for(int i = 0; i < VarValues.size(); i++)
		{
			VarValues[i]= m_factors[i] * VarValues[i];
		}

		std::string func = "";
		int current_depth = 0;
		int number = 0;

		std::unordered_map<std::string, double> ValueMap;

		std::string CarriedForward = "";
		for(int i = 0; i < VarValues.size(); i++)
		{
			if(m_VarDepth[i] > current_depth)
			{
				while(current_depth != m_VarDepth[i])
				{
					current_depth++;
					func = func + "(";
				}
			}
			else if(m_VarDepth[i] < current_depth)
			{
				while(current_depth != m_VarDepth[i])
				{
					current_depth--;
					func = func + ")";
				}
			}
			
			func = func + CarriedForward;
			CarriedForward = "";
			std::string key = "x" + std::to_string(number);
			func = func + key;
			ValueMap[key] = VarValues[i];
			number++;
			
			if(i != VarValues.size()-1)
			{
				std::string temp;
				switch(m_op[i])
				{
				case InternalOperations::p:
					temp = "+";
					break;
				case InternalOperations::d:
					temp =  "/";
					break;
				case InternalOperations::mu:
					temp = "*";
					break;
				case InternalOperations::mi:
					temp = "-";
					break;
				}

				if(m_VarDepth[i+1] >= m_VarDepth[i])
				{
					func = func + temp;
				}
				else
				{
					CarriedForward = temp;
				}
			}
			
		}

		//std::cout << func << std::endl; //for debug purposes
		
		//using reverse polish notation 
		std::stack<char> st;
		std::string Postfix;

		auto prec = [](char c) {
			if(c == '/' || c == '*')
				return 2;
			else if(c == '+' || c == '-')
				return 1;
			else 
				return -1;
		};

		std::vector<char> chr = {'+', '-', '/', '*', '(', ')'};
		for(int i = 0; i < func.length(); i++)
		{
			char c = func[i];
			auto it = std::find(chr.begin(), chr.end(), c);
			if(it == chr.end())
			{
				Postfix += c;
			}
			else if(c == '(')
			{
				st.push('(');
			}
			else if(c == ')')
			{
				while (st.top() != '(')
				{
					Postfix += st.top();
					st.pop();
				}
				st.pop();
			}
			else
			{
				while(!st.empty() && prec(c) < prec(st.top()) || !st.empty() && prec(c) == prec(st.top()))
				{
					Postfix += st.top();
					st.pop();
				}
				st.push(c);
			}
		}

		while(!st.empty())
		{
			Postfix += st.top();
			st.pop();
		}

		//std::cout << Postfix << std::endl;

		//evaluate using the same method
		std::stack<double> ValueStack;
		for(int i = 0; i < Postfix.length(); i++)
		{
			char c = Postfix[i];
			auto it = std::find(chr.begin(), chr.end(), c);
			if(it == chr.end())
			{
				std::string VarName = "";
				VarName += c;
				i++;
				c = Postfix[i];
				while(c != 'x' && std::find(chr.begin(), chr.end(), c) == chr.end())
				{
					VarName += c;
					i++;
					c = Postfix[i];
				}
				i--;
				ValueStack.push(ValueMap[VarName]);
			}
			else
			{
				double Val2 = ValueStack.top();
				ValueStack.pop();
				double Val1 = ValueStack.top();
				ValueStack.pop();
				double Total;

				switch (c)
				{
				case '+':
					Total = Val1 + Val2;
					break;
				case '-':
					Total = Val1 - Val2;
					break;
				case '*':
					Total = Val1 * Val2;
					break;
				case '/':
					Total = Val1 / Val2;
					break;
				default:
					Total = 1.0;
					break;
				}

				ValueStack.push(Total);
			}
		}
		return ValueStack.top();
		
	}

	Function::Function(FuncPtr FunctionPtr, ReturnType type, std::string name, std::string var, double factor)
	{
		m_func = FunctionPtr;
		m_funcType = type;
		m_FunctionName = name;
		m_variables = {var};
		m_factors = {factor};
		m_op = {};
		m_VarDepth = {};
	}

	Function::Function(FuncPtr FunctionPtr, ReturnType type, std::string name, std::vector<std::string> vars, std::vector<double> factors)
	{
		m_func = FunctionPtr;
		m_funcType = type;
		m_FunctionName = name;
		m_variables = vars;
		m_factors = factors;
		m_op = {};
		m_VarDepth = {};
	}

	std::vector<std::string> Function::GetVariable()
	{
		return m_variables;
	}

	std::string Function::GetName()
	{
		return m_FunctionName;
	}

	std::shared_ptr<Function> FunctionParser(std::string& func, std::string& var)
	{
		bool value = true;

		std::vector<double> factors;
		std::vector<std::string> variables;
		std::vector<Function::InternalOperations> operations;
		std::vector<int> VarDepth;
		bool VarPreFactor = true;
		bool Initial = false;
		bool decimal = false;
		std::pair<double,double> nums = {0.0, 0.0};
		std::vector<char> SpecialCharacters = {'+', '-', '*', '/','(', ')'};
		int VarNum = 0;
		int depth = 0;

		std::string buffer = "";

		for(auto c = var.cbegin(); c != var.cend(); c++)
		{
			auto it = std::find(SpecialCharacters.begin(), SpecialCharacters.end(), (*c)); 
			if((*c) == '.' || (*c) == ',')
			{
				decimal = true;
				nums.first = std::stod(buffer);
				buffer = "";
			}
			else if(it != SpecialCharacters.end())
			{
				if(VarPreFactor && buffer == "")
				{
					VarNum++;
					variables.push_back("");
					VarDepth.push_back(depth);
					auto TempIt = it;

					if((*it) == '-')
					{
						factors.push_back(-1.0);
						it = it + 1;
					}
					else if((*it) == '+')
					{
						factors.push_back(1.0);
						it = it + 2;
					}
					else
					{
						factors.push_back(1.0);
					}

					if(*(c-1) == ')')
					{
						VarNum--;
						variables.pop_back();
						VarDepth.pop_back();
						factors.pop_back();	
						it = TempIt;
					}
				}
				else if(VarPreFactor && buffer != "")
				{
					VarNum++;
					variables.push_back("");
					VarDepth.push_back(depth);

					double PreFactor = std::stod(buffer);
					if(!decimal)
					{
						nums.first = PreFactor;
					}
					else
					{
						int DivideBy = buffer.length();
						PreFactor = PreFactor * std::pow(10.0, -1.0 * (double)DivideBy);
						nums.second = PreFactor;
					}

					factors.push_back(nums.first + nums.second);
					buffer = "";
					decimal = false;
					nums = {0,0};
				}
				else
				{
					VarNum++;
					variables.push_back(buffer);
					VarDepth.push_back(depth);
					buffer = "";
				}

				//if(VarNum == 1)
				//{
				//	if(operations.size() == 0)
				//	{
				//		operations.push_back(Function::InternalOperations::p);
				//	}
				//}
				//else if(VarDepth[VarNum-1] == VarDepth[VarNum-2] + 1)
				//{
				//	auto d = *(c-1);
				//	if(*(c-1) == '(')
				//	{
				//		VarDepth.pop_back();
				//		variables.pop_back();
				//		factors.pop_back();
				//		VarNum = VarNum-1;
				//		Initial = false;
				//	}
				//	else if(Initial)
				//	{
				//		operations.push_back(Function::InternalOperations::p);
				//		Initial = false;
				//	}
				//}
				if(VarDepth[VarNum-1] == VarDepth[VarNum-2] - 1)
				{
					VarDepth.pop_back();
					variables.pop_back();
					factors.pop_back();
					VarNum = VarNum-1;
				}

				if((*it) == '+')
				{
					operations.push_back(Function::InternalOperations::p);
				}
				else if((*it) == '-')
				{
					operations.push_back(Function::InternalOperations::mi);
				}
				else if((*it) == '*')
				{
					operations.push_back(Function::InternalOperations::mu);
				}
				else if((*it) == '/')
				{
					operations.push_back(Function::InternalOperations::d);
				}
				else if((*it) == '(')
				{
					Initial = true;
					depth = depth + 1;
					if(operations[VarNum - 2] != Function::InternalOperations::mu && operations[VarNum - 2] != Function::InternalOperations::d)
					{
						operations.push_back(Function::InternalOperations::mu);
					}
					else
					{
						VarDepth.pop_back();
						variables.pop_back();
						factors.pop_back();
						VarNum = VarNum-1;	
					}			
				}
				else if((*it) == ')')
				{
					depth = depth - 1;
				}

				VarPreFactor = true;
			}
			else if(std::isdigit((*c)) == 0 && VarPreFactor)
			{
				if(buffer == "")
				{
					factors.push_back(1.0);
				}
				else
				{
					double PreFactor = std::stod(buffer);
					if(!decimal)
					{
						nums.first = PreFactor;
					}
					else
					{
						int DivideBy = buffer.length();
						PreFactor = PreFactor * std::pow(10.0, -1.0 * (double)DivideBy);
						nums.second = PreFactor;
					}
					factors.push_back(nums.first + nums.second);
				}
				buffer = (*c);
				decimal = false;
				VarPreFactor = false;
				nums = {0,0};

			}
			else
			{
				buffer += (*c);
			}
		}
		
		if(buffer != "")
		{
			if(!VarPreFactor)
			{
				variables.push_back(buffer);
				VarDepth.push_back(depth);
				buffer = "";
			}
			else
			{
				double PreFactor = std::stod(buffer);
				if(!decimal)
				{
					nums.first = PreFactor;
				}
				else
				{
					int DivideBy = buffer.length();
					PreFactor = PreFactor * std::pow(10.0, -1.0 * (double)DivideBy);
					nums.second = PreFactor;
				}
				factors.push_back(nums.first + nums.second);

				variables.push_back("");
				VarDepth.push_back(depth);
				buffer = "";
			}

		}

		//for(auto c = var.cbegin(); c != var.cend(); c++)
		//{
		//	if((*c) == ' ')
		//	{
		//		continue;
		//	}
		//	else if((*c) == '.')
		//	{
		//		decimal = true;
		//		nums.first = std::stod(buffer);
		//		buffer = "";
		//	}
		//	else if(std::isdigit((*c)) == 0 && value)
		//	{
		//		value = false;
		//		double val = std::stod(buffer);
		//		if(!decimal)
		//		{
		//			nums.first = val;
		//		}
		//		else
		//		{
		//			int DivideBy = buffer.length();
		//			val = val * std::pow(10.0, -1.0 * (double)DivideBy);
		//			nums.second = val;
		//		}
		//		buffer = (*c);
		//	}
		//	else
		//	{
		//		buffer += (*c);
		//	}
		//}
//
		//std::string variable = buffer;
		//double factor = nums.first + nums.second;

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
			 _func = std::make_shared<Function>(MathematicalFunctions::sin, Function::ReturnType::d, FunctionName, variables, factors);
		}
		if(FunctionName.compare("cos") == 0)
		{
			_func = std::make_shared<Function>(MathematicalFunctions::cos, Function::ReturnType::d, FunctionName, variables, factors);
		}
		if(FunctionName.compare("scalar") == 0)
		{
			_func = std::make_shared<Function>(MathematicalFunctions::scalar, Function::ReturnType::d, FunctionName, variables, factors);
		}
		
		_func->SetOp(operations);
		_func->SetVarDepth(VarDepth);

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

		void* scalar(void* value) //scalar not scaler
		{
			double _val = *(double*)value ;
			double* v = &_val;
			return (void*)v;
		}
	}
}