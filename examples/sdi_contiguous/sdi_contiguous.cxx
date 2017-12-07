/*
  Simple test of algorithm for determining contiguous blocks from a vector
  of ints (i.e. an SDI region)
 */
#include <utility>
#include <vector>
#include <chrono>

#include <bout/physicsmodel.hxx>
#include "output.hxx"

class SDIContiguous : public PhysicsModel {
public:

  int init(bool restarting) {
    typedef std::chrono::time_point<std::chrono::steady_clock> SteadyClock;
    typedef std::chrono::duration<double> Duration;
    using namespace std::chrono;


    int maxBlockSize;
    int numLoops;
    Options *globalOptions = Options::getRoot();
    Options *modelOpts=globalOptions->getSection("contig");
    OPTION(modelOpts, maxBlockSize, 4);
    OPTION(modelOpts, numLoops, 100);
    
    mesh->coordinates()->geometry();

    Field3D f; f=0.; 

    std::vector<std::pair<int,int>> contiguousBlocks;
    //auto regionAll = f.sdi_region(RGN_NOBNDRY).region;
    auto regionAll = f.sdi_region(RGN_ALL).region;

    SteadyClock start1 = steady_clock::now();
    for (int x = 0; x < numLoops; x++) {
    
      int position = 0;
      const int lastPoint = regionAll.size();
      
      while (position < lastPoint){
	const auto startPair = position;
	int count = 1; //We will always have at least startPair in the block so count starts at 1
      
	//Consider if the next point should be added to this block
	for(position++; count<maxBlockSize; position++){
	  if((regionAll[position]-regionAll[position-1])==1){
	    count++;
	  }else{//Reached the end of this block so break
	    break;
	  }
	}

	output_debug<<"Building contiguous block from "<<startPair<<" up to but not including "<<position<<". This is a block of size "<<count<<"\n";

	//Add pair to vector, denotes start and size of blocks
	contiguousBlocks.push_back({startPair, count});
      }
    }
    Duration elapsed1 = steady_clock::now() - start1;
    
    int totalElements = 0;
    for (const auto& block: contiguousBlocks){
      totalElements += block.second;
    }
    totalElements/=numLoops;
    output_debug<<"Constructed "<<contiguousBlocks.size()<<" blocks corresponding to "<<totalElements<<" elements\n";
    
    output<<"This took "<<elapsed1.count() / numLoops<<"s per try and "<<elapsed1.count() / (numLoops*regionAll.size())<<"s per try per point.\n";

    output_warn<<"NumReg\tLoops\tblockSz\tBlocks\tinBlck\tTime/try\tTime/(try*npt)\n";
    output_warn<<regionAll.size()<<"\t";
    output_warn<<numLoops<<"\t";
    output_warn<<maxBlockSize<<"\t";
    output_warn<<contiguousBlocks.size() / numLoops<<"\t";
    output_warn<<totalElements<<"\t";
    output_warn<< elapsed1.count() / numLoops <<"\t";
    output_warn<< elapsed1.count() / (numLoops*regionAll.size()) <<"\t";
    output_warn<<"\n";
    return 0; //Exit now
   
  }
  int rhs(BoutReal time) {
    return 0;
  }
  
private:

};


BOUTMAIN(SDIContiguous);
