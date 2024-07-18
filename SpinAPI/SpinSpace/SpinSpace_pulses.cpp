/////////////////////////////////////////////////////////////////////////
// SpinSpace class (SpinAPI Module) - developed by Irina Anisimova and Luca Gerhards.
// ------------------
// This source file contains methods for managing pulses
// and having rotating functions.
// Molecular Spin Dynamics Software - developed by Claus Nielsen and Luca Gerhards..
// (c) 2024 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
namespace SpinAPI
{
    bool SpinSpace::PulseOperator(const pulse_ptr &_pulse, arma::cx_mat &_out) const
    {
        // Make sure the pulse is valid
        if (_pulse == nullptr)
            return false;

        // Create temporary matrix to hold the result
        arma::cx_mat tmp = arma::zeros<arma::cx_mat>(this->HilbertSpaceDimensions(), this->HilbertSpaceDimensions());

        if (_pulse->Type() == PulseType::InstantPulse)
        {
            // Create the rotation angle
            double angle;

            if (!this->CreateRotAngle(_pulse->Angle(), angle))
            {
                std::cout << "Failed to create a rotation angle for the pulse." << std::endl;
            }

            // Get rotation axis
            auto rotvec = _pulse->Rotationaxis();

            // this is should be fuction CreateRotDirection
            arma::cx_mat Sx;
            arma::cx_mat Sy;
            arma::cx_mat Sz;

            auto spinlist = _pulse->Group();

            for (auto i = spinlist.cbegin(); i < spinlist.cend(); i++)
            {
                this->CreateOperator(arma::conv_to<arma::cx_mat>::from((*i)->Sx()), (*i), Sx);
                this->CreateOperator(arma::conv_to<arma::cx_mat>::from((*i)->Sy()), (*i), Sy);
                this->CreateOperator(arma::conv_to<arma::cx_mat>::from((*i)->Sz()), (*i), Sz);

                tmp += Sx * rotvec[0] + Sy * rotvec[1] + Sz * rotvec[2];
            }

            // Check whether we want a superspace or Hilbert space result
            if (this->useSuperspace)
            {
                arma::cx_mat lhs;
                arma::cx_mat rhs;
                lhs = arma::expmat((arma::cx_double(0.0, -1.0) * angle * tmp));
                rhs = arma::expmat((arma::cx_double(0.0, 1.0) * angle * tmp));

                _out = kron(lhs, rhs.st());
            }
            else
            {
                // We already have the result in the Hilbert space
                _out = arma::expmat((arma::cx_double(0.0, -1.0) * angle * tmp));
            }
        }
        else if (_pulse->Type() == PulseType::LongPulse)
        {
            // this is should be fuction CreateRotDirection
            arma::cx_mat Sx;
            arma::cx_mat Sy;
            arma::cx_mat Sz;

            arma::vec field;
            field = _pulse->Field();

            auto spinlist = _pulse->Group();

            for (auto i = spinlist.cbegin(); i < spinlist.cend(); i++)
            {
                this->CreateOperator(arma::conv_to<arma::cx_mat>::from((*i)->Sx()), (*i), Sx);
                this->CreateOperator(arma::conv_to<arma::cx_mat>::from((*i)->Sy()), (*i), Sy);
                this->CreateOperator(arma::conv_to<arma::cx_mat>::from((*i)->Sz()), (*i), Sz);

                tmp += Sx * field(0) + Sy * field(1) + Sz * field(2);
            }

            // Multiply with the given prefactor (or 1.0 if none was specified)
            tmp *= _pulse->Prefactor();

            // Multiply with common prefactor
            if (_pulse->AddCommonPrefactor())
                tmp *= 8.794e+1;

            // Check whether we want a superspace or Hilbert space result
            if (this->useSuperspace)
            {
                arma::cx_mat lhs;
                arma::cx_mat rhs;
                auto result = this->SuperoperatorFromLeftOperator(tmp, lhs);
                result &= this->SuperoperatorFromRightOperator(tmp, rhs);
                if (result)
                    _out = lhs - rhs;
                else
                    return false;
            }
            else
            {
                // We already have the result in the Hilbert space
                _out = tmp;
            }
        }

        return true;
    }

    // Create a rotational angle
    bool SpinSpace::CreateRotAngle(double _angle_deg, double &result) const
    {
        if (_angle_deg <= 360.0 || _angle_deg >= -360.0)
        {
            result = _angle_deg * M_PI / 180.0;
        }
        else
        {
            std::cout << "Pulse angle is not recognised. Please choose in an interval [-360,360] degrees." << std::endl;
        }

        return true;
    }
}