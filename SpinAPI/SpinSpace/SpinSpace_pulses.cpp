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
    bool SpinSpace::PulseOperator(const pulse_ptr &_pulse, arma::sp_cx_mat &_out) const // Instant Pulse in SS
    {
        // Check whether we want a superspace result
        if (!this->useSuperspace)
        {
            std::cout << "Failed to create a rotation pulse in SS." << std::endl;
            return false;
        }

        // Make sure the pulse is valid
        if (_pulse == nullptr)
            return false;

        // Create temporary matrix to hold the result
        arma::cx_mat tmp = arma::zeros<arma::cx_mat>(this->HilbertSpaceDimensions(), this->HilbertSpaceDimensions());

        if (_pulse->Type() == PulseType::InstantPulse)
        {
            // Create the rotation angle
            double angle = 0;
            if (!this->CreateRotAngle(_pulse->Angle(), angle))
            {
                std::cout << "Failed to create a rotation angle for the pulse." << std::endl;
            }

            // Get rotation axis
            auto rotvec = _pulse->Rotationaxis();

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

            // Create pulse operator in SS
            arma::cx_mat lhs;
            arma::cx_mat rhs;

            this->SuperoperatorFromLeftOperator(tmp, lhs);
            this->SuperoperatorFromRightOperator(tmp, rhs);
            _out = arma::conv_to<arma::sp_cx_mat>::from(arma::expmat((arma::cx_double(0.0, -1.0) * angle * (lhs - rhs))));

            // /////Alternative way/////
            // lhs = arma::expmat((arma::cx_double(0.0, -1.0) * angle * tmp));
            // rhs = arma::expmat((arma::cx_double(0.0, 1.0) * angle * tmp));
            // _out = kron(lhs, rhs.st());
        }
        else if (_pulse->Type() == PulseType::LongPulse || _pulse->Type() == PulseType::LongPulseStaticField)
        {

            arma::cx_mat Sx;
            arma::cx_mat Sy;
            arma::cx_mat Sz;

            arma::vec field;
            field = _pulse->Field();

            arma::vec prefactorList;
            std::vector<bool> ignoreTensorsList;
            std::vector<bool> addCommonPrefactorList;
            prefactorList = _pulse->PrefactorList();
            ignoreTensorsList = _pulse->IgnoreTensorsList();
            addCommonPrefactorList = _pulse->AddCommonPrefactorList();

            auto spinlist = _pulse->Group();

            arma::cx_mat tmp_part = arma::zeros<arma::cx_mat>(this->HilbertSpaceDimensions(), this->HilbertSpaceDimensions());

            int n = 0;
            for (auto i = spinlist.cbegin(); i < spinlist.cend(); i++)
            {
                if (!ignoreTensorsList.empty())
                {
                    if (ignoreTensorsList[n])
                    {
                        this->CreateOperator(arma::conv_to<arma::cx_mat>::from((*i)->Sx()), (*i), Sx);
                        this->CreateOperator(arma::conv_to<arma::cx_mat>::from((*i)->Sy()), (*i), Sy);
                        this->CreateOperator(arma::conv_to<arma::cx_mat>::from((*i)->Sz()), (*i), Sz);
                    }
                    else
                    {
                        this->CreateOperator(arma::conv_to<arma::cx_mat>::from((*i)->Tx()), (*i), Sx);
                        this->CreateOperator(arma::conv_to<arma::cx_mat>::from((*i)->Ty()), (*i), Sy);
                        this->CreateOperator(arma::conv_to<arma::cx_mat>::from((*i)->Tz()), (*i), Sz);
                    }
                }
                else
                {
                    this->CreateOperator(arma::conv_to<arma::cx_mat>::from((*i)->Sx()), (*i), Sx);
                    this->CreateOperator(arma::conv_to<arma::cx_mat>::from((*i)->Sy()), (*i), Sy);
                    this->CreateOperator(arma::conv_to<arma::cx_mat>::from((*i)->Sz()), (*i), Sz);
                }

                tmp_part = Sx * field(0) + Sy * field(1) + Sz * field(2);

                if (!addCommonPrefactorList.empty())
                {
                    if (addCommonPrefactorList[n])
                    {
                        tmp_part *= 8.79410005e+1;
                    }
                }

                if (!prefactorList.is_empty())
                {
                    tmp += prefactorList(n) * tmp_part;
                }
                else
                {
                    tmp += tmp_part;
                }
                n++;
            }

            arma::cx_mat lhs;
            arma::cx_mat rhs;
            this->SuperoperatorFromLeftOperator(tmp, lhs);
            this->SuperoperatorFromRightOperator(tmp, rhs);
            _out = arma::conv_to<arma::sp_cx_mat>::from(lhs - rhs);
        }
        else
        {
            std::cout << "Not implemented yet. Sorry." << std::endl;
        }

        return true;
    }

    bool SpinSpace::PulseOperator(const pulse_ptr &_pulse, arma::sp_cx_mat &_left, arma::sp_cx_mat &_right) const // Instant Pulse in HS
    {
        // Check whether we want a Hilbert space result
        if (this->useSuperspace)
        {
            std::cout << "Failed to create a rotation pulse in HS." << std::endl;
            return false;
        }

        // Make sure the pulse is valid
        if (_pulse == nullptr)
            return false;

        // Create temporary matrix to hold the result
        arma::cx_mat tmp = arma::zeros<arma::cx_mat>(this->HilbertSpaceDimensions(), this->HilbertSpaceDimensions());

        if (_pulse->Type() == PulseType::InstantPulse)
        {
            // Create the rotation angle
            double angle = 0;

            if (!this->CreateRotAngle(_pulse->Angle(), angle))
            {
                std::cout << "Failed to create a rotation angle for the pulse." << std::endl;
            }

            // Get rotation axis
            auto rotvec = _pulse->Rotationaxis();

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

            // Create left and right pulse operators (use as _left*rho*_right)
            _left = arma::conv_to<arma::sp_cx_mat>::from(arma::expmat((arma::cx_double(0.0, -1.0) * angle * tmp)));
            _right = arma::conv_to<arma::sp_cx_mat>::from(arma::expmat((arma::cx_double(0.0, 1.0) * angle * tmp)));
        }
        else if (_pulse->Type() == PulseType::LongPulse || _pulse->Type() == PulseType::LongPulseStaticField)
        {
            arma::cx_mat Sx;
            arma::cx_mat Sy;
            arma::cx_mat Sz;

            arma::vec field;
            field = _pulse->Field();

            arma::vec prefactorList;
            std::vector<bool> ignoreTensorsList;
            std::vector<bool> addCommonPrefactorList;
            prefactorList = _pulse->PrefactorList();
            ignoreTensorsList = _pulse->IgnoreTensorsList();
            addCommonPrefactorList = _pulse->AddCommonPrefactorList();

            auto spinlist = _pulse->Group();

            arma::cx_mat tmp_part = arma::zeros<arma::cx_mat>(this->HilbertSpaceDimensions(), this->HilbertSpaceDimensions());

            int n = 0;
            for (auto i = spinlist.cbegin(); i < spinlist.cend(); i++)
            {
                if (!ignoreTensorsList.empty())
                {
                    if (ignoreTensorsList[n])
                    {
                        this->CreateOperator(arma::conv_to<arma::cx_mat>::from((*i)->Sx()), (*i), Sx);
                        this->CreateOperator(arma::conv_to<arma::cx_mat>::from((*i)->Sy()), (*i), Sy);
                        this->CreateOperator(arma::conv_to<arma::cx_mat>::from((*i)->Sz()), (*i), Sz);
                    }
                    else
                    {
                        this->CreateOperator(arma::conv_to<arma::cx_mat>::from((*i)->Tx()), (*i), Sx);
                        this->CreateOperator(arma::conv_to<arma::cx_mat>::from((*i)->Ty()), (*i), Sy);
                        this->CreateOperator(arma::conv_to<arma::cx_mat>::from((*i)->Tz()), (*i), Sz);
                    }
                }
                else
                {
                    this->CreateOperator(arma::conv_to<arma::cx_mat>::from((*i)->Sx()), (*i), Sx);
                    this->CreateOperator(arma::conv_to<arma::cx_mat>::from((*i)->Sy()), (*i), Sy);
                    this->CreateOperator(arma::conv_to<arma::cx_mat>::from((*i)->Sz()), (*i), Sz);
                }

                tmp_part = Sx * field(0) + Sy * field(1) + Sz * field(2);

                if (!addCommonPrefactorList.empty())
                {
                    if (addCommonPrefactorList[n])
                    {
                        tmp_part *= 8.79410005e+1;
                    }
                }

                if (!prefactorList.is_empty())
                {
                    tmp += prefactorList(n) * tmp_part;
                }
                else
                {
                    tmp += tmp_part;
                }
                n++;
            }

            // Create left and right pulse operators (use as exp(-i*_left*t)*rho*exp(-i*_right*t))
            _left = arma::conv_to<arma::sp_cx_mat>::from(tmp);
            _right = arma::conv_to<arma::sp_cx_mat>::from(tmp.st());
        }
        else
        {
            std::cout << "Not implemented yet. Sorry." << std::endl;
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
