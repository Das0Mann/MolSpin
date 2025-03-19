/////////////////////////////////////////////////////////////////////////
// ActionLogSpace (RunSection module)
// ------------------
//
// Molecular Spin Dynamics Software - developed by Claus Nielsen and Luca Gerhards.
// (c) 2025 Quantum Biology and Computational Physics Group.
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
        // Data members
        ActionScalar *actionScaler;
        arma::rowvec m_Points;
        int m_Num; // number of points
        int m_Step;
        std::pair<double, double> m_Bounds; // upper and lower bounds for the log space
    private:
        bool CalculatePoints(int, double, double); // generates a list of numbers between log base 10 num1 and log base 10 num2
        bool GetPoint(double &);                   // gets the current value and assignes it to the double& variable
    protected:
        // overwritten protected methods
        bool DoStep() override;
        bool DoValidate() override;
        bool Reset() override;

    public:
        ActionLogSpace(const MSDParser::ObjectParser &, const std::map<std::string, ActionScalar> &, const std::map<std::string, ActionVector> &); // Normal constructor
    };
}

#endif