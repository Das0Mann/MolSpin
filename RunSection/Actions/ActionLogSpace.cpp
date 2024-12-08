#include "ActionLogSpace.h"

bool RunSection::ActionLogSpace::CalculatePoints(int n, double start, double stop)
{
    m_Points = std::make_shared<arma::rowvec>(arma::logspace<arma::rowvec>(start,stop,n));
    return true;
}

bool RunSection::ActionLogSpace::DoStep()
{
    // Make sure we have an ActionScaler to act on
	if (actionScaler == nullptr || !this->IsValid())
	{
		return false;
	}

    double val = 0;
    GetPoint(val);
    
    if(!this->actionScaler->Set(val))
    {
        return false;
    }
    return true;
}

bool RunSection::ActionLogSpace::DoValidate()
{
    std::string str;
	if (!this->Properties()->Get("actionscaler", str) && !this->Properties()->Get("scaler", str))
	{
		std::cout << "ERROR: No ActionScaler specified for the LogSpace action \"" << this->Name() << "\"!" << std::endl;
		return false;
	}
    int NumPoints = 0;
    if(!this->Properties()->Get("points", NumPoints))
    {
        std::cout << "ERROR: No Number of points specified for the LogSpace action \"" << this->Name() << "\"!" << std::endl;
		return false;
    }
    double lowbounds = 0;
    if(!this->Properties()->Get("minvalue", lowbounds) && !this->Properties()->Get("min", lowbounds))
    {
        std::cout << "ERROR: No miniumum value specified for the LogSpace action \"" << this->Name() << "\"!" << std::endl;
		return false;
    }
    double upbounds = 0;
    if(!this->Properties()->Get("maxvalue", upbounds) && !this->Properties()->Get("max", upbounds))
    {
        std::cout << "ERROR: No maximum value specified for the LogSpace action \"" << this->Name() << "\"!" << std::endl;
		return false;
    }

    if(lowbounds >= upbounds)
    {
        std::cout << "ERROR: No incorrect bounds specified for the LogSpace action \"" << this->Name() << "\"!" << std::endl;
		return false;
    }

    m_Num = NumPoints;
    m_Bounds = {lowbounds,upbounds};

    // Attemp to set the ActionVector
    if (!this->Scalar(str, &(this->actionScaler)))
	{
		std::cout << "ERROR: Could not find ActionScaler \"" << str << "\" specified for the LogSpace action \"" << this->Name() << "\"!" << std::endl;
		return false;
	}

    if(!this->actionScaler->IsReadonly())
    {
		std::cout << "ERROR: Read only ActionScaler \"" << str << "\" specified for the LogSpace action \"" << this->Name() << "\"! Cannot act on this scaler!" << std::endl;
		return false;
	}


    if(!this->DoStep())
    {
        return false;
    }
    return true;

}

bool RunSection::ActionLogSpace::Reset()
{
    m_Step = 0
    double val = 0;
    GetPoint(val);
    this->actionScaler->Set(val);
    return true;
}

RunSection::ActionLogSpace::ActionLogSpace(const MSDParser::ObjectParser &_parser, const std::map<std::string, ActionScalar> &_scaler, const std::map<std::string, ActionVector> &_vector)
    :Action(_parser, _scaler, _vector), actionScaler(nullptr)
{
    m_Num = 0;
    m_Bounds = {0.0,0.0};
    m_Points = arma::rowvec(arma::zeros(m_Num));
}


bool RunSection::ActionLogSpace::GetPoint(double& val)
{
    val = m_Points->[m_Step];    
    m_Step++;
    return true;
}
