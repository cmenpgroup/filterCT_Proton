#include "Riostream.h"
#include "TApplication.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TClasTool.h"
#include "TIdentificator.h"
#include "TMath.h"
#include "TString.h"
#include "TMath.h"
#include "massConst.h"
using namespace std;

#define MAX_ELECTRONS 1
#define MAX_PROTONS 1

#define PDG_PHOTON 22
#define PDG_POSITRON -11
#define PDG_ELECTRON 11
#define PDG_PIPLUS 211
#define PDG_PIMINUS -211
#define PDG_KPLUS 321
#define PDG_KMINUS -321
#define PDG_NEUTRON 2112
#define PDG_PROTON 2212

#define GEANT3_PHOTON 1
#define GEANT3_POSITRON 2
#define GEANT3_ELECTRON 3
#define GEANT3_PIPLUS 8
#define GEANT3_PIMINUS 9
#define GEANT3_KPLUS 11
#define GEANT3_KMINUS 12
#define GEANT3_NEUTRON 13
#define GEANT3_PROTON 14

#define CUT_Q2 1.0
#define CUT_W 2.0
#define CUT_NU 0.85

#define EBEAM 5.015  // e- beam energy in GeV

//declarations of functions
void PrintAnalysisTime(float tStart, float tStop);
void PrintUsage(char *processName);
int GetPID(string partName, int kind);

typedef struct{
    Float_t EvtNum, ElecVertTarg, Q2, Nu, Xb, W;
    Float_t Xcorr, Ycorr, Zcorr;
    Int_t nElec, nPip, nPim, nGam, nProton, nNeutron, nKp, nKm, nPositron;
    Int_t PartComb;
} KINEVAR;

typedef struct{
    int Sector;
    float Charge, Pid, Beta;
    float Px, Py, Pz, Mom, Mass2;
    float X, Y, Z;
    float ECx, ECy, ECz, ECu, ECv, ECw, ECtot, ECin, ECout, ECtime, ECpath;
    float EChit_M2, EChit_M3, EChit_M4, Chi2EC;
    float SCpath, SCtime;
    float CCnphe;
    float T, Xf, Mx2, Pt, Zh, ThetaPQ, PhiPQ, TimeCorr4;
} PARTVAR;

PARTVAR SetPARTVAR(TIdentificator *t, int index, int kind, bool simul_key);

int main(int argc, char **argv)
{
    gROOT->SetBatch(true);
    extern char *optarg;
    int c;
    extern int optind;

    int i, j, k;
    int nRows, kind, tempPid;
    int candCtr = 0;
    int combCtr = 0;
    int ctr_nElec = 0;
    int dEvents = 1000; // increment of events for processing print statement
    int MaxEvents = 0; // max. number of events to process
    int nfiles = 0; // number of processed files
    int minRows = MAX_ELECTRONS + MAX_PROTONS; // number of particles in an event to start filtering

    TString catPid;

    bool bBatchMode = false;    // quiet mode
    bool simul_key = false;  // simulation flag (true = simulation, false = data)
    bool mflag = true; // cut flag for GetCategorization(k,tt,mflag)
    int cat_key = 0; // PID categorization 0 = EVNT (default), 1 = Full
    int tgt_key = 1;  // intitialize target flag 1 = Carbon (default), 2 = Iron, 3 = Lead
    string target; // solid target name

    char *inFile;
    string outFile = "protons.root";

    bool topology = false;
    vector<int> elecIndex;
    vector<int> protonIndex;

    TVector3 *vert;

    float timeStart = clock(); // start time

    TClasTool *input = new TClasTool();
    input->InitDSTReader("ROOTDSTR");

    TIdentificator *t = new TIdentificator(input);

    for (i = 0; i < argc; ++i) cerr << argv[i] << " "; cerr << endl;
    while ((c = getopt(argc,argv, "o:M:D:c:t:Sih")) != -1 ) {
        switch (c) {
            case 'o': outFile = optarg; break;
            case 'M': MaxEvents = atoi(optarg); break;
            case 'D': dEvents = atoi(optarg); break;
            case 'c': cat_key = atoi(optarg); break;
            case 't': tgt_key = atoi(optarg); break;
            case 'S': simul_key = true; break;
            case 'i': bBatchMode = true; break;
            case 'h':
                PrintUsage(argv[0]);
                exit(0);
                break;

            default:
                cerr << "Unrecognized argument: " << optarg << endl;
                PrintUsage(argv[0]);
                exit(0);
                break;
        }
    }

    TFile *output; // ROOT output file

    // check target selection
    switch(tgt_key){
        case 1: target = "C"; break;
        case 2: target = "Fe"; break;
        case 3: target = "Pb"; break;
        default: cout<<"Unknown target "<<target<<endl; exit(0); break;
    }
    cout<<"Analyzing " << target << " target data"<<endl;

    string kineList = "EvtNum/F:ElecVertTarg/F:Q2/F:Nu/F:Xb/F:W:Xcorr/F:Ycorr/F:Zcorr/F:nElec/I:nPip/I:nPim/I:nGam/I:nProton/I:nNeutron/I:nKp/I:nKm/I:nPositron/I:PartComb/I";

    string partList = "Sector/I:Charge/F:Pid/F:Beta/F:Px/F:Py/F:Pz/F:Mom/F:Mass2/F:X/F:Y/F:Z/F:ECx/F:ECy/F:ECz/F:ECu/F:ECv/F:ECw/F:ECtot/F:ECin/F:ECout/F:ECtime/F:ECpath/F:EChit_M2/F:EChit_M3/F:EChit_M4/F:Chi2EC/F:SCpath/F:SCtime/F:CCnphe/F:T/F:Xf/F:Mx2/F:Pt/F:Zh/F:ThetaPQ/F:PhiPQ/F:TimeCorr4/F";
    KINEVAR myKine;
    PARTVAR myElec;
    PARTVAR myProton;

    TTree *dataTree = new TTree("Data","Experimental Data Tree");
    dataTree->Branch("Kinematics",&myKine,kineList.c_str());
    dataTree->Branch("Electron",&myElec,partList.c_str());
    dataTree->Branch("Proton",&myProton,partList.c_str());

    output = new TFile(outFile.c_str(), "RECREATE", "Experimental Data");

    for (i = optind; i < argc; ++i) {
        inFile = argv[i]; // process all arguments on command line.
        if (*inFile != '-') { // we have a file to process
            cout << "Analyzing file " << inFile << endl; // let user know which file is being processed

            input->Add(inFile); // read file into ClasTool object

            nfiles++; // increment file counter
        }
    }

    Long_t nEntries = (Long_t) input->GetEntries(); // get total number of events

    cout<<"Analyzing "<<nEntries<<" from "<<nfiles<< " files."<<endl; // print out stats

    input->Next();

    k = 0; // event counter

    if(MaxEvents == 0) MaxEvents = nEntries; // if user does not set max. number of events, set to nEntries

    while (k < MaxEvents) {
    	if (!bBatchMode && ((k % dEvents) == 0)){
    		cerr << k << "\r";
    	}

        if(simul_key){
            kind = 1;
            nRows = input->GetNRows("GSIM");
        }else{
            kind = 0;
            nRows = input->GetNRows("EVNT");
        }

        if(nRows >= minRows){
      		myKine.nElec = 0;
      		myKine.nPip = 0;
	    	  myKine.nPim = 0;
	    	  myKine.nGam = 0;
	    	  myKine.nProton = 0;
          myKine.nNeutron = 0;
          myKine.nKp = 0;
          myKine.nKm = 0;
          myKine.nPositron = 0;
          elecIndex.clear(); // clear out the electron list
          protonIndex.clear(); // clear out the protons list

          topology = false; // init. the event topology cut
	        for (j = 0; j < nRows; j++) {

                // select the PID selection scheme
                if(simul_key){
                    catPid = t -> GetCategorizationGSIM(j);
                }else{
                    switch(cat_key){
                        case 0: catPid = t -> GetCategorizationEVNT(j); break;
                        case 1: catPid = t -> GetCategorization(j,target.c_str(),mflag); break;
                        default: cout<<"Incorrect PID categorization.  Try again."<<endl; exit(0); break;
                    }
                }
                tempPid = t -> Id(j,kind);

                if(catPid.EqualTo("electron")){
                    myKine.nElec++;
                    elecIndex.push_back(j);
                    ctr_nElec++;
                }
                if(catPid.EqualTo("proton")){
                    myKine.nProton++;
                    protonIndex.push_back(j);
                }

                // using PDG id numbers for both EVNT and GSIM
                if(tempPid == GetPID("Photon",0)) myKine.nGam++;
                if(tempPid == GetPID("PiPlus",0)) myKine.nPip++;
                if(tempPid == GetPID("PiMinus",0)) myKine.nPim++;
                if(tempPid == GetPID("Neutron",0)) myKine.nNeutron++;
                if(tempPid == GetPID("KPlus",0)) myKine.nKp++;
                if(tempPid == GetPID("KMinus",0)) myKine.nKm++;
                if(tempPid == GetPID("Positron",0)) myKine.nPositron++;
            }

            // check event topology
            topology = (myKine.nElec>=MAX_ELECTRONS && myKine.nProton>=MAX_PROTONS);

	    	    if(topology && t->Q2(kind) > CUT_Q2 && t->W(kind) > CUT_W && t->Nu(kind)/EBEAM < CUT_NU) {
                candCtr++;
                myKine.EvtNum = t -> NEvent();
                myKine.ElecVertTarg = t -> ElecVertTarg(kind);
                myKine.Q2 = t -> Q2(kind);
		     	      myKine.Nu = t -> Nu(kind);
	       		    myKine.Xb = t -> Xb(kind);
        		    myKine.W = t -> W(kind);

                if(simul_key){
                    myKine.Xcorr = t->X(0, kind);
                    myKine.Ycorr = t->Y(0, kind);
                    myKine.Zcorr = t->Z(0, kind);
                }else{
                    vert = t->GetCorrectedVert();
                    myKine.Xcorr = vert->X();
                    myKine.Ycorr = vert->Y();
                    myKine.Zcorr = vert->Z();
                }

                for(unsigned iElec=0; iElec<elecIndex.size(); iElec++){
                    myElec = SetPARTVAR(t, elecIndex.at(iElec), kind, simul_key);
                    for(unsigned iProton=0; iProton<protonIndex.size(); iProton++){
                        myProton = SetPARTVAR(t, protonIndex.at(iProton), kind, simul_key);
                        myKine.PartComb = 100*(iElec+1) + iProton + 1;
                        combCtr++;
                        dataTree->Fill();
                    }
                }
            }
        }
    	  k++; // increment event counter
        input->Next();
    }
//    dataTree->Print();
//    dataTree->Scan("Kinematics.EvtNum:Electron.Pid:PiPlus.Pid:PiMinus.Pid:Photon1.Pid:Photon2.Pid:Photon2.Beta");

    dataTree->Write();
    cout<<"Candidate data events = "<<candCtr<<endl;
    cout<<"Combination data events = "<<combCtr<<endl;
    cout<<"#(e-) = "<<ctr_nElec<<endl;

    output->Write();
    output->Close();

    float timeStop = clock();
    PrintAnalysisTime(timeStart,timeStop);

    return 0;
}

void PrintUsage(char *processName)
{
    cerr << processName << " <options> <filename>\n";
    cerr << "\toptions are:\n";
    cerr << "\t-o<filename>\tROOT output file (def. = protons.root).\n";
    cerr << "\t-M#\t\tprocess maximum # of events.\n";
    cerr << "\t-D#\t\tinform user when # of events have been processed (def. = 1000).\n";
    cerr << "\t-c#\t\tType in categoization # scheme of 0=EVNT or 1=Full (def. = 0).\n";
    cerr << "\t-t#\t\tTarget # of 1=C, 2=Fe, or 3=Pb (def. = 1).\n";
    cerr << "\t-S\t\tAnalyze simulation.\n";
    cerr << "\t-i\t\tquiet mode (no counter).\n";
    cerr << "\t-h\t\tprint the above" << endl;
}


void PrintAnalysisTime(float tStart, float tStop){
    //time to complete function
    float minutes = 0;
    float seconds = 0;
    minutes = (tStop - tStart)/1000000;
    minutes = (minutes)/60;
    seconds = fmod(minutes,1);
    minutes = minutes-seconds;
    seconds = seconds*60;

    if (minutes==0){
        cout<<endl<<"Completed in "<<seconds<<" seconds."<<endl<<endl;
    }
    else{
        cout<<endl<<"Completed in "<<minutes<<" minutes and "<<seconds<<" seconds."<<endl<<endl;
    }
}

int GetPID(string partName, int kind){

    int ret = 0;

    if(kind==0){
        if(partName.compare("Electron")==0){
            ret = PDG_ELECTRON;
        }else if(partName.compare("Positron")==0){
            ret = PDG_POSITRON;
        }else if(partName.compare("Photon")==0){
            ret = PDG_PHOTON;
        }else if(partName.compare("PiPlus")==0){
            ret = PDG_PIPLUS;
        }else if(partName.compare("PiMinus")==0){
            ret = PDG_PIMINUS;
        }else if(partName.compare("KPlus")==0){
            ret = PDG_KPLUS;
        }else if(partName.compare("KMinus")==0){
            ret = PDG_KMINUS;
        }else if(partName.compare("Neutron")==0){
            ret = PDG_NEUTRON;
        }else if(partName.compare("Proton")==0){
            ret = PDG_PROTON;
        }else{
            cerr<<"GetPid(): Unknown PDG particle "<<partName.c_str()<<endl; exit(0);
        }
    }else if(kind==1){
        if(partName.compare("Electron")==0){
            ret = GEANT3_ELECTRON;
        }else if(partName.compare("Positron")==0){
            ret = GEANT3_POSITRON;
        }else if(partName.compare("Photon")==0){
            ret = GEANT3_PHOTON;
        }else if(partName.compare("PiPlus")==0){
            ret = GEANT3_PIPLUS;
        }else if(partName.compare("PiMinus")==0){
            ret = GEANT3_PIMINUS;
        }else if(partName.compare("KPlus")==0){
            ret = GEANT3_KPLUS;
        }else if(partName.compare("KMinus")==0){
            ret = GEANT3_KMINUS;
        }else if(partName.compare("Neutron")==0){
            ret = GEANT3_NEUTRON;
        }else if(partName.compare("Proton")==0){
            ret = GEANT3_PROTON;
        }else{
            cerr<<"GetPid(): Unknown GEANT3 particle "<<partName.c_str()<<endl; exit(0);
        }
    }else{
        cerr<<"GetPID: Unknown analysis channel "<<kind<<endl;
    }
    return ret;
}

PARTVAR SetPARTVAR(TIdentificator *t, int index, int kind, bool simul_key){

    TVector3 *ECxyz = new TVector3(0.0,0.0,0.0);
    TVector3 *ECuvw;
    PARTVAR myPart;

    myPart.Sector = t->Sector(index,kind);
    myPart.Charge = t->Charge(index,kind);
    myPart.Beta = t->Betta(index,kind);
    myPart.Pid = t->Id(index,kind);
    myPart.Mom = t->Momentum(index,kind);
    myPart.Px = t->Px(index, kind);
    myPart.Py = t->Py(index, kind);
    myPart.Pz = t->Pz(index, kind);
    myPart.X = t->X(index, kind);
    myPart.Y = t->Y(index, kind);
    myPart.Z = t->Z(index, kind);
    myPart.Mass2 = t->Mass2(index,kind);

    myPart.ThetaPQ = t -> ThetaPQ(index, kind);
    myPart.PhiPQ = t -> PhiPQ(index, kind);
    myPart.Zh = t -> Zh(index, kind);
    myPart.Pt = TMath::Sqrt(t -> Pt2(index, kind));
    myPart.Mx2 = t -> Mx2(index, kind);
    myPart.Xf = t -> Xf(index, kind);
    myPart.T = t -> T(index, kind);

    // initialize detector info
    myPart.ECtot = 0;
    myPart.ECin = 0;
    myPart.ECout = 0;
    myPart.ECx = 0;
    myPart.ECy = 0;
    myPart.ECz = 0;
    myPart.ECu = 0;
    myPart.ECv = 0;
    myPart.ECw = 0;
    myPart.ECtime = 0;
    myPart.ECpath = 0;
    myPart.EChit_M2 = 0;
    myPart.EChit_M3 = 0;
    myPart.EChit_M4 = 0;
    myPart.Chi2EC = 0;
    myPart.SCtime = 0;
    myPart.SCpath = 0;
    myPart.CCnphe = 0;
    myPart.TimeCorr4 = 0;

    if(simul_key == false){
        myPart.ECtot = TMath::Max(t->Etot(index),t->Ein(index)+t->Eout(index));
        myPart.ECin = t->Ein(index);
        myPart.ECout = t->Eout(index);
        myPart.ECx = t->XEC(index);
        myPart.ECy = t->YEC(index);
        myPart.ECz = t->ZEC(index);
        ECxyz->SetXYZ(t->XEC(index),t->YEC(index),t->ZEC(index));
        ECuvw = t->XYZToUVW(ECxyz);
        myPart.ECu = ECuvw->X();
        myPart.ECv = ECuvw->Y();
        myPart.ECw = ECuvw->Z();
        myPart.ECtime = t->TimeEC(index);
        myPart.ECpath = t->PathEC(index);
        myPart.EChit_M2 = t->EChit_Moment2(index);
        myPart.EChit_M3 = t->EChit_Moment3(index);
        myPart.EChit_M4 = t->EChit_Moment4(index);
        myPart.Chi2EC = t->Chi2EC(index);

        myPart.SCtime = t->TimeSC(index);
        myPart.SCpath = t->PathSC(index);

        myPart.CCnphe = t->Nphe(index);

        if(myPart.Pid == GetPID("Electron",kind)) myPart.TimeCorr4 = t -> TimeCorr4(0.000511,index);
        if(myPart.Pid == GetPID("PiPlus",kind)) myPart.TimeCorr4 = t -> TimeCorr4(kMassPi_plus,index);
        if(myPart.Pid == GetPID("PiMinus",kind)) myPart.TimeCorr4 = t -> TimeCorr4(kMassPi_min,index);
        if(myPart.Pid == GetPID("Proton",kind)) myPart.TimeCorr4 = t -> TimeCorr4(0.938,index);
        if(myPart.Pid == GetPID("Photon",kind)) myPart.TimeCorr4 = t -> TimeCorr4(0.0,index);
    }

    return myPart;
}
