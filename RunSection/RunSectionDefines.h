/////////////////////////////////////////////////////////////////////////
// RunSection Defines (RunSection Module)
// ------------------
// Definitions used by the RunSection Module.
//
// Molecular Spin Dynamics Software - developed by Claus Nielsen and Luca Gerhards.
// (c) 2025 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
#ifndef MOD_RunSection_Defines
#define MOD_RunSection_Defines

namespace RunSection
{
	// Message type used by the OutputHandler class
	using MessageType = unsigned int;

	// A bitmask with all bits for normal message categories set to 1
	const MessageType NormalMessage = 7;

	// Normal message categories
	const MessageType MessageType_Details = 1;
	const MessageType MessageType_Normal = 2;
	const MessageType MessageType_Important = 3;
	const MessageType MessageType_Critical = 4;

	// Special message categories
	const MessageType MessageType_Warning = 8;
	const MessageType MessageType_Error = 16;

	// The notification level that is used if nothing else is specified
	const MessageType DefaultNotificationLevel = MessageType_Normal | MessageType_Warning | MessageType_Error;
}

#endif
