/////////////////////////////////////////////////////////////////////////
// SubSystem class (SpinAPI Module)
// ------------------
// TODO: Write a description of the SubSystem class 
// 
//
// Molecular Spin Dynamics Software - developed by Claus Nielsen and Luca Gerhards.
// (c) 2019 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
#ifndef MOD_SpinAPI_SubSystem
#define MOD_SpinAPI_SubSystem

#include <memory>
#include "MSDParserfwd.h"
#include "SpinAPIfwd.h"

namespace SpinAPI
{
    class SubSystem
    {
    public:
        struct Transition
        {
            std::shared_ptr<transition_ptr> TransitionObject;
            int type; //0 transtion out of the system, 1 transiton into a other system, -1 unknown;
            std::string name;

            std::string source;
            std::string target;

            Transition(std::shared_ptr<transition_ptr> t, std::string n)
                :TransitionObject(t), type(-1), name(n), source(""), target("")
            {

            }
        };   
    private:
        // Implementation details
        std::string m_name;
		system_ptr m_system;
        std::vector<std::shared_ptr<spin_ptr>> m_spins;
        std::vector<std::shared_ptr<interaction_ptr>> m_interactions;
        std::vector<Transition> m_transtions;
        std::vector<std::shared_ptr<operator_ptr>> m_operators;
        std::vector<std::shared_ptr<pulse_ptr>> m_pulses;

        std::shared_ptr<MSDParser::ObjectParser> m_properties;

    public:
        SubSystem(std::string, std::string, system_ptr); //Normal constructor
        //SubSystem(const SubSystem& ) = delete;//Copy constructor
        ~SubSystem(); //Destructor

        //Operator overloading 
        const SubSystem &operator=(const SubSystem& ) = delete; //Copy Assignment 

        std::string Name() const;
        bool Validate(); //validates all the objects asssinged to the subsystem exist

        // Allow access to custom properties to be used for custom tasks
		//std::shared_ptr<const MSDParser::ObjectParser> Properties() const;

        inline std::vector<Transition> GetTransitions() {
            return m_transtions;
        }

        inline std::vector<Transition>& GetTransitions_Ref() {
            return m_transtions;
        }

        std::vector<interaction_ptr> GetInteractions();

        void ModifyTransitionVec(int index, Transition NewTransition)
        {
            m_transtions[index] = NewTransition;
        }

        std::vector<std::string> GetSpinNames();
        
    };

    using subsystem_ptr = std::shared_ptr<SubSystem>;

    bool LinkTransitions(system_ptr);
}

#endif