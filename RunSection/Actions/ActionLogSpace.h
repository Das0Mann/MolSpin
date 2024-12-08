/////////////////////////////////////////////////////////////////////////
// ActionLogSpace (RunSection module)
// ------------------
//
// Molecular Spin Dynamics Software - developed by Claus Nielsen and Luca Gerhards.
// (c) 2019 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
#ifndef MOD_RunSection_ActionLogSpace
#define MOD_RunSection_ActionLogSpace

#include "Action.h"
#include <armadillo>
#include <memory>

namespace RunSection
{
    class ActionLogSpace : public Action
    {
    private:
        //Data members
        ActionScalar* actionScaler;
        
        std::shared_ptr<arma::rowvec> m_Points;
        int m_Num;
        std::pair<double,double> m_Bounds;
    private:
        bool CalculatePoints(int, double, double);
        bool GetPoint(double&);
    protected:
        //overwritten protected methods
        bool DoStep() override;
		bool DoValidate() override;
		bool Reset() override;
    public:
        ActionLogSpace(const MSDParser::ObjectParser &, const std::map<std::string, ActionScalar> &, const std::map<std::string, ActionVector> &); // Normal constructor
    };
}