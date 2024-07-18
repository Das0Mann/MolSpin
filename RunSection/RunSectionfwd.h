/////////////////////////////////////////////////////////////////////////
// ForwardDeclarations (MSDParser Module)
// ------------------
// Forward declarations of classes etc. to decrease file inter-
// dependencies and speedup compilation time.
//
// NOTE: Make sure this file is included AFTER the other includes!
//
// Molecular Spin Dynamics Software - developed by Claus Nielsen and Luca Gerhards.
// (c) 2019 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
#ifndef MOD_RunSection_ForwardDeclarations
#define MOD_RunSection_ForwardDeclarations

namespace RunSection
{
#ifndef MOD_RunSection_RunSection
	class RunSection;
#endif

#ifndef MOD_RunSection_BasicTask
	class BasicTask;
#endif

#ifndef MOD_RunSection_Action
	class Action;
#endif

#ifndef MOD_RunSection_Settings
	class Settings;
#endif
}

#endif
