#include <iostream>

#include <getopt.h>
#include <cstring>

#include "TFile.h"
#include "TTree.h"

#include "XiaData.hpp"
#include "ChannelEvent.hpp"

// Local files
#include "TraceQdcScan.hpp"

#ifdef USE_HRIBF
#include "GetArguments.hpp"
#include "Scanor.hpp"
#include "ScanorInterface.hpp"
#endif

// Define the name of the program.
#ifndef PROG_NAME
#define PROG_NAME "TraceQdcScan"
#endif

///////////////////////////////////////////////////////////////////////////////
// class TraceQdcScanUnpacker
///////////////////////////////////////////////////////////////////////////////

/** Process all events in the event list.
  * \param[in]  addr_ Pointer to a location in memory. 
  * \return Nothing.
  */
void TraceQdcScanUnpacker::ProcessRawEvent(ScanInterface *addr_/*=NULL*/){
	if(!addr_){ return; }
	
	XiaData *current_event = NULL;
	
	// Fill the processor event deques with events
	while(!rawEvent.empty()){
		current_event = rawEvent.front();
		rawEvent.pop_front(); // Remove this event from the raw event deque.

#ifdef USE_HRIBF		
		// If using scanor, output to the generic histogram so we know that something is happening.
		count1cc_(8000, (current_event->modNum*16+current_event->chanNum), 1);
#endif	
	
		// Check that this channel event exists.
		if(!current_event){ continue; }

		// Send the event to the scan interface object for processing.
		if(addr_->AddEvent(current_event))
			addr_->ProcessEvents();
	}
	
	// Finish up with this raw event.
	addr_->ProcessEvents();
}

///////////////////////////////////////////////////////////////////////////////
// class TraceQdcScanScanner
///////////////////////////////////////////////////////////////////////////////

/// Default constructor.
TraceQdcScanScanner::TraceQdcScanScanner() : ScanInterface(), init(false), mult(0), setMod(0), setChan(0), lowMax(20), highMax(40), outFile(NULL), outTree(NULL) {
}

/// Destructor.
TraceQdcScanScanner::~TraceQdcScanScanner(){
	if(init){
		outFile->cd();
		outTree->Write();
		outFile->Close();
	
		delete outFile;

		/*for(size_t i = 0; i < lowMax; i++){
			delete[] traceQDC[i];
		}
		delete[] traceQDC;*/
	}
}

/** ExtraCommands is used to send command strings to classes derived
  * from ScanInterface. If ScanInterface receives an unrecognized
  * command from the user, it will pass it on to the derived class.
  * \param[in]  cmd_ The command to interpret.
  * \param[out] arg_ Vector or arguments to the user command.
  * \return True if the command was recognized and false otherwise.
  */
bool TraceQdcScanScanner::ExtraCommands(const std::string &cmd_, std::vector<std::string> &args_){
	return false;
}

/** ExtraArguments is used to send command line arguments to classes derived
  * from ScanInterface. This method should loop over the optionExt elements
  * in the vector userOpts and check for those options which have been flagged
  * as active by ::Setup(). This should be overloaded in the derived class.
  * \return Nothing.
  */
void TraceQdcScanScanner::ExtraArguments(){
	if(userOpts.at(0).active){
		lowMax = strtoul(userOpts.at(0).argument.c_str(), NULL, 0);
		std::cout << msgHeader << "Set low side integration limit to " << lowMax << " bins.\n";
	}
	if(userOpts.at(1).active){
		lowMax = strtoul(userOpts.at(1).argument.c_str(), NULL, 0);
		std::cout << msgHeader << "Set high side integration limit to " << lowMax << " bins.\n";
	}
	if(userOpts.at(2).active){
		setMod = strtoul(userOpts.at(2).argument.c_str(), NULL, 0);
		std::cout << msgHeader << "Set module to " << setMod << ".\n";
	}	
	if(userOpts.at(3).active){
		setChan = strtoul(userOpts.at(3).argument.c_str(), NULL, 0);
		std::cout << msgHeader << "Set channel to " << setChan << ".\n";
	}
}

/** CmdHelp is used to allow a derived class to print a help statement about
  * its own commands. This method is called whenever the user enters 'help'
  * or 'h' into the interactive terminal (if available).
  * \param[in]  prefix_ String to append at the start of any output. Not used by default.
  * \return Nothing.
  */
void TraceQdcScanScanner::CmdHelp(const std::string &prefix_/*=""*/){
	std::cout << "   mycmd1                   - A useful terminal command.\n";
	std::cout << "   mycmd2 [param]           - A useful terminal command with an optional argument.\n";
	std::cout << "   mycmd3 <param>           - A useful terminal command with one argument.\n";
	std::cout << "   mycmd4 <param1> <param2> - A useful terminal command with two arguments.\n";
}

/** ArgHelp is used to allow a derived class to add a command line option
  * to the main list of options. This method is called at the end of
  * from the ::Setup method.
  * Does nothing useful by default.
  * \return Nothing.
  */
void TraceQdcScanScanner::ArgHelp(){
	AddOption(optionExt("low", required_argument, NULL, 'L', "<lowLimit>", "Set the furthest bin from the peak maximum on the low side (default 20)."));
	AddOption(optionExt("high", required_argument, NULL, 'H', "<highLimit>", "Set the furthest bin from the peak maximum on the high side (default 20)."));
	AddOption(optionExt("mod", required_argument, NULL, 'M', "<module>", "Set the pixie module number."));
	AddOption(optionExt("chan", required_argument, NULL, 'C', "<channel>", "Set the pixie channel number."));
	
	// Note that the following single character options are reserved by ScanInterface
	//  b, h, i, o, q, s, and v
}

/** SyntaxStr is used to print a linux style usage message to the screen.
  * \param[in]  name_ The name of the program.
  * \return Nothing.
  */
void TraceQdcScanScanner::SyntaxStr(char *name_){ 
	std::cout << " usage: " << std::string(name_) << " [options]\n"; 
}

/** Initialize the map file, the config file, the processor handler, 
  * and add all of the required processors.
  * \param[in]  prefix_ String to append to the beginning of system output.
  * \return True upon successfully initializing and false otherwise.
  */
bool TraceQdcScanScanner::Initialize(std::string prefix_){
	if(init){ return false; }

	std::string fname = this->GetOutputFilename();
	if(fname.empty()){
		std::cout << msgHeader << "Warning! No output filename specified. Using \"traceQdcScan.root\".\n";
		fname = "traceQdcScan.root";
	}

	// Setup the output root file.
	outFile = new TFile(fname.c_str(), "RECREATE");

	if(!outFile || !outFile->IsOpen()){
		std::cout << msgHeader << "Error! Failed to open output root file \"" << fname << "\".\n";
		return false;	
	}

	// Initialize the TTree.
	outTree = new TTree("data", "QDC scanner tree");
	outTree->Branch("traceQDC[20][40]", traceQDC, "traceQDC[20][40]/f");

	// Initialize the QDC array.
	/*traceQDC = new float*[lowMax];
	for(size_t i = 0; i < lowMax; i++){
		traceQDC[i] = new float[highMax];	
	}*/

	return (init = true);
}

/** Peform any last minute initialization before processing data. 
  * /return Nothing.
  */
void TraceQdcScanScanner::FinalInitialization(){
	// Do some last minute initialization before the run starts.
}

/** Receive various status notifications from the scan.
  * \param[in] code_ The notification code passed from ScanInterface methods.
  * \return Nothing.
  */
void TraceQdcScanScanner::Notify(const std::string &code_/*=""*/){
	if(code_ == "START_SCAN"){  }
	else if(code_ == "STOP_SCAN"){  }
	else if(code_ == "SCAN_COMPLETE"){ std::cout << msgHeader << "Scan complete.\n"; }
	else if(code_ == "LOAD_FILE"){ std::cout << msgHeader << "File loaded.\n"; }
	else if(code_ == "REWIND_FILE"){  }
	else{ std::cout << msgHeader << "Unknown notification code '" << code_ << "'!\n"; }
}

/** Return a pointer to the Unpacker object to use for data unpacking.
  * If no object has been initialized, create a new one.
  * \return Pointer to an Unpacker object.
  */
Unpacker *TraceQdcScanScanner::GetCore(){ 
	if(!core){ core = (Unpacker*)(new TraceQdcScanUnpacker()); }
	return core;
}

/** Add a channel event to the deque of events to send to the processors.
  * This method should only be called from TraceQdcScanUnpacker::ProcessRawEvent().
  * \param[in]  event_ The raw XiaData to add to the channel event deque.
  * \return False.
  */
bool TraceQdcScanScanner::AddEvent(XiaData *event_){
	if(!event_){ return false; }

	unsigned int mod = event_->GetModuleNumber();
	unsigned int chan = event_->GetChannelNumber();

	// Check if this is the correct channel.
	if(mod == setMod && chan == setChan){ // Handle the individual XiaData.
		mult++;

		ChannelEvent *chEvent = new ChannelEvent(event_);

		/*for (size_t i = 0; i < lowMax; i++) {
			for (size_t j = 0; j < highMax; j++) {
				traceQDC[i][j] = chEvent->IntegratePulse(chEvent->max_index - i, chEvent->max_index + j);	
			}
		}*/

		for (size_t i = 0; i < 20; i++) {
			for (size_t j = 0; j < 40; j++) {
				traceQDC[i][j] = chEvent->IntegratePulse(chEvent->max_index - i, chEvent->max_index + j);	
			}
		}

		// Fill the TTree.
		outTree->Fill();

		// We're done with the event now.
		delete chEvent;
	}
	else{ // Do nothing with the event.
		delete event_;
	}
	
	return false;
}

/** Process all channel events read in from the rawEvent.
  * This method should only be called from TraceQdcScanUnpacker::ProcessRawEvent().
  * \return False.
  */
bool TraceQdcScanScanner::ProcessEvents(){
	// Process all of the events added so far.
	return false;
}

#ifndef USE_HRIBF
int main(int argc, char *argv[]){
	// Define a new unpacker object.
	TraceQdcScanScanner scanner;
	
	// Set the output message prefix.
	scanner.SetProgramName(std::string(PROG_NAME));	
	
	// Initialize the scanner.
	if(!scanner.Setup(argc, argv))
		return 1;

	// Run the main loop.
	int retval = scanner.Execute();
	
	scanner.Close();
	
	return retval;
}
#else
TraceQdcScanScanner *scanner = NULL;

// Do some startup stuff.
extern "C" void startup_()
{
	scanner = new TraceQdcScanScanner();	

	// Handle command line arguments.
	scanner->Setup(GetNumberArguments(), GetArguments());
	
	// Get a pointer to a class derived from Unpacker.
	ScanorInterface::get()->SetUnpacker(scanner->GetCore());
}

///@brief Defines the main interface with the SCANOR library, the program
/// essentially starts here.
///@param [in] iexist : unused paramter from SCANOR call
extern "C" void drrsub_(uint32_t &iexist) {
	drrmake_();
	hd1d_(8000, 2, 256, 256, 0, 255, "Run DAMM you!", strlen("Run DAMM you!"));
	endrr_();
}

// Catch the exit call from scanor and clean up c++ objects CRT
extern "C" void cleanup_()
{
	// Do some cleanup.
	std::cout << "\nCleaning up..\n";
	scanner->Close();
	delete scanner;
}
#endif
