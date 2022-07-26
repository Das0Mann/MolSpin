/////////////////////////////////////////////////////////////////////////
// ForwardDeclarations (MSDParser Module)
// ------------------
// Forward declarations of classes etc. to decrease file inter-
// dependencies and speedup compilation time.
// 
// NOTE: Make sure this file is included AFTER the other includes!
// 
// Molecular Spin Dynamics Software - developed by Claus Nielsen.
// (c) 2019 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
#ifndef MOD_MSDParser_ForwardDeclarations
#define MOD_MSDParser_ForwardDeclarations

namespace MSDParser
{
	#ifndef MOD_MSDParser_MSDParser
	class MSDParser;
	#endif
	
	#ifndef MOD_MSDParser_ObjectParser
	class ObjectParser;
	#endif
}

#endif
