/////////////////////////////////////////////////////////////////////////
// ActionFibonacciSphere implementation (RunSection module)
//
// Molecular Spin Dynamics Software - developed by Claus Nielsen and Luca Gerhards.
// (c) 2019 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
#include "ActionFibonacciSphere.h"
#include "ObjectParser.h"

#include<array>

namespace RunSection
{
    ActionFibonacciSphere::ActionFibonacciSphere(const MSDParser::ObjectParser & _parser, const std::map<std::string, ActionScalar> & _scalars, const std::map<std::string, ActionVector> & _vectors)
        :Action(_parser, _scalars, _vectors), actionVector(nullptr)
    {
        m_Points = nullptr;
        m_Step = 0;
        m_Magnitude = 1.0;
    }

    ActionFibonacciSphere::~ActionFibonacciSphere()
    {
        if(m_Points != nullptr)
        {
            free(m_Points);
        }
    }

    bool ActionFibonacciSphere::CalculatePoints(int n)
    {
        m_Points = (point*)malloc(n*sizeof(point));
        m_Num = n;

        if(m_Points == NULL)
        {
            std::cout << "Memory not allocated" << std::endl;
            return false;
        }

        double phi = M_PI * (3.0 - std::sqrt(5.0)); //golden angle in radians;
        for(int i = 0; i < n; i++)
        {
            double y = 1.0 - ((double)i / double(n - 1)) * 2;
            double theta = phi * (double)i;

            m_Points[i] = {y, theta};
        }
        return true;
    }

    bool ActionFibonacciSphere::GetPoint(std::array<double,3>& arr)
    {
        if(m_Step >= m_Num)
        {
            return false;
        }

        auto p = m_Points[m_Step];

        float y = p.first;
        float theta = p.second;

        double r = std::sqrt(1.0 - (y * y));
        double x = m_Magnitude * std::cos(theta) * r;
        double z = m_Magnitude * std::sin(theta) * r;

        m_Step++;

        y = y * m_Magnitude;

        arr = {x,y,z};
        return true;
    }

    bool ActionFibonacciSphere::DoStep()
    {

        // Make sure we have an ActionVector to act on
		if (actionVector == nullptr || !this->IsValid())
		{
			return false;
		}

        // Retrieve the vector we want to change
        std::array<double,3> points;
        GetPoint(points);
		arma::vec vec = {m_Magnitude*points[0], m_Magnitude*points[1], m_Magnitude*points[2]};

        return this->actionVector->Set(vec);
    }

    bool ActionFibonacciSphere::DoValidate()
    {
        std::string str;
		if (!this->Properties()->Get("actionvector", str) && !this->Properties()->Get("vector", str))
		{
			std::cout << "ERROR: No ActionVector specified for the FibonacciSphere action \"" << this->Name() << "\"!" << std::endl;
			return false;
		}

        int NumPoints = 0;
        if(!this->Properties()->Get("points", NumPoints))
        {
            std::cout << "ERROR: No Number of points specified for the FibonacciSphere action \"" << this->Name() << "\"!" << std::endl;
			return false;
        }

        m_Num = NumPoints;
        CalculatePoints(m_Num);

        // Attemp to set the ActionVector
        if (!this->Vector(str, &(this->actionVector)))
		{
			std::cout << "ERROR: Could not find ActionVector \"" << str << "\" specified for the FibonacciSphere action \"" << this->Name() << "\"!" << std::endl;
			return false;
		}

		// Readonly ActionTargets cannot be acted on
		if (this->actionVector->IsReadonly())
		{
			std::cout << "ERROR: Read only ActionVector \"" << str << "\" specified for the FibonacciSphere action \"" << this->Name() << "\"! Cannot act on this vector!" << std::endl;
			return false;
		}

        auto vec1 = this->actionVector->Get();
        m_Magnitude = 0;
        
        for(int i = 0; i < 3; i++)
        {
            m_Magnitude += std::pow(vec1[i],2);
        }

        m_Magnitude = std::sqrt(m_Magnitude);

        std::array<double,3> points;
        GetPoint(points);

		arma::vec vec = {m_Magnitude*points[0], m_Magnitude*points[1], m_Magnitude*points[2]};
        if(!this->actionVector->Set(vec))
        {
            return false;
        }

        return true;
    }

    bool ActionFibonacciSphere::Reset()
	{
		m_Step = 0;
        std::array<double,3> points;
        GetPoint(points);
		arma::vec vec = {points[0], points[1], points[2]};
        this->actionVector->Set(vec);
	}
}