#include <string>

#include "lcio.h"
#include "EVENT/MCParticle.h"
#include "IMPL/MCParticleImpl.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TNtupleD.h"
#include "TSystem.h"
#include "TPad.h"
#include "TStyle.h"

#include "lcfiplus.h"
#include "process.h"
#include "testproc.h"
#include "VertexSelector.h"
#include "algoEtc.h"
#include "VertexFinderSuehara.h"
#include "VertexFitterSimple.h"

using namespace lcfiplus;

namespace lcfiplus {

const Jet* JetMCMatch(JetVec& jets, const MCParticle* mcp, vector<const Track*>& assignedTracks, vector<const Track*>& residualTracks) {
  const vector<const Track*>* pTracks;
  pTracks = &(Event::Instance()->getTracks());

  vector<const Track*> bTracks;

  vector<int> nTrackInJet;
  nTrackInJet.resize(jets.size());
  vector<int> nVertexTrackInJet;
  nVertexTrackInJet.resize(jets.size());

  int nvtx = 0;
  // get tracks
  for (unsigned int i=0; i<pTracks->size(); i++) {
    const MCParticle* mcpc = (*pTracks)[i]->getMcp();

    if (mcpc==0)continue;
    if (mcpc->isParent(mcp)) {
      bTracks.push_back((*pTracks)[i]);
      for (unsigned int j=0; j<jets.size(); j++) {
        for (unsigned int k=0; k<jets[j]->getVertices().size(); k++) {
          const vector<const Track*>& vtr = jets[j]->getVertices()[k]->getTracks();
          if (find(vtr.begin(), vtr.end(), (*pTracks)[i]) != vtr.end()) {
            // the track matched to this jet
            nTrackInJet[j] ++;
// 						if(nvtx == 0 && k > 0){cout << "CAUTION: vertices in the same jet might be from different semistables!" << endl;}
// 						if(nvtx > 0 && k == 0){cout << "CAUTION: vertices in different jets might be from the same samistable!" << endl;}
            nVertexTrackInJet[j] ++;
            nvtx ++;
          }
        }
        if (find(jets[j]->getTracks().begin(), jets[j]->getTracks().end(), (*pTracks)[i]) != jets[j]->getTracks().end()) {
          // the track matched to this jet
          nTrackInJet[j] ++;
        }
      }
    }
  }

  int ntijMax = 0;
  int ntijMaxIndex = -1;

  int ntrsum = 0;
  int ntrvtxsum = 0;
  // determine best-match jet
  for (unsigned int j=0; j<jets.size(); j++) {
    if (ntijMax < nTrackInJet[j]) {
      ntijMax = nTrackInJet[j];
      ntijMaxIndex = j;
    }
    ntrsum += nTrackInJet[j];
    ntrvtxsum += nVertexTrackInJet[j];
  }

  if (ntijMaxIndex == -1)
    return 0;

  vector<const Track*> jetTracks = jets[ntijMaxIndex]->getTracks();
  for (unsigned int i=0; i<jets[ntijMaxIndex]->getVertices().size(); i++) {
    const Vertex* vtx = jets[ntijMaxIndex]->getVertices()[i];
    jetTracks.insert(jetTracks.end(), vtx->getTracks().begin(), vtx->getTracks().end());
  }

  // obtain assignedtracks
  for (unsigned int i=0; i<bTracks.size(); i++) {
    if (find(jetTracks.begin(), jetTracks.end(), bTracks[i]) != jetTracks.end())
      assignedTracks.push_back(bTracks[i]);
    else
      residualTracks.push_back(bTracks[i]);
  }

  cout << "Assigned jet " << ntijMaxIndex << ", Vertex tracks: ";
  for (unsigned int j=0; j<jets.size(); j++)cout << nVertexTrackInJet[j] << ",";
  cout << "/" << ntrvtxsum << ", all tracks: ";
  for (unsigned int j=0; j<jets.size(); j++)cout << nTrackInJet[j] << ",";
  cout << "/" << bTracks.size() << ", PDG: " << mcp->getPDG() << endl;//"," ;
//	cout << " ( " << mcp->getVertex().x() << " " << mcp->getVertex().y() << " " << mcp->getVertex().z() << ")" << endl;

  return jets[ntijMaxIndex];
}

void TestAlgoV0::init(Parameters* param) {
  Algorithm::init(param);
  string filename = param->get("FileName",string("testv0.root"));
  _vtxname = param->get("VertexCollectionName",string("BuildUpVertex"));
  _file = new TFile(filename.c_str(),"RECREATE");
  _ntp = new TTree("v0","v0");
  _vertices = 0;

  VtxData& d = _data;
  _ntp->Branch("x",&d.x,"x/D");
  _ntp->Branch("y",&d.y,"y/D");
  _ntp->Branch("z",&d.z,"z/D");
  _ntp->Branch("r",&d.r,"r/D");
  _ntp->Branch("cs",&d.cs,"cs/D");
  _ntp->Branch("phi",&d.phi,"phi/D");
  _ntp->Branch("chrg",&d.chrg,"chrg/D");
  _ntp->Branch("dirdot",&d.dirdot,"dirdot/D");
  _ntp->Branch("dirdot2",&d.dirdot2,"dirdot2/D");
  _ntp->Branch("ntrk",&d.ntrk,"ntrk/I");
  _ntp->Branch("mks",&d.mks,"mks/D");
  _ntp->Branch("ml0",&d.ml0,"ml0/D");
  _ntp->Branch("mconv",&d.mconv,"mconv/D");
  _ntp->Branch("mks2",&d.mks2,"mks2/D");
  _ntp->Branch("ml02",&d.ml02,"ml02/D");
  _ntp->Branch("v0",&d.v0,"v0/I");
  _ntp->Branch("ks",&d.ks,"ks/I");
  _ntp->Branch("l0",&d.l0,"l0/I");
  _ntp->Branch("conv",&d.conv,"conv/I");
  _ntp->Branch("mcpdg1",&d.mcpdg1,"mcpdg1/I");
  _ntp->Branch("mcpdg2",&d.mcpdg2,"mcpdg2/I");
  _ntp->Branch("mcppdg1",&d.mcppdg1,"mcppdg1/I");
  _ntp->Branch("mcppdg2",&d.mcppdg2,"mcppdg2/I");
  _ntp->Branch("mcpp1",&d.mcpp1,"mcpp1/I");
  _ntp->Branch("mcpp2",&d.mcpp2,"mcpp2/I");
}

void TestAlgoV0::process() {
  if (!_vertices) {
    Event::Instance()->Get(_vtxname.c_str(), _vertices);
  }

  const VertexVec& vtx_list = *_vertices;
  for (unsigned int i=0; i < vtx_list.size(); ++i) {
    const Vertex* vtx = vtx_list[i];
    if (vtx->isPrimary()) continue;

    memset(&_data,0,sizeof(_data));

    _data.x = vtx->getX();
    _data.y = vtx->getY();
    _data.z = vtx->getZ();
    TVector3 pos = vtx->getPos();
    _data.r = pos.Mag();
    _data.cs = pos.CosTheta();
    _data.phi = pos.Phi();

    TVector3 mom;
    TVector3 mom2;
    for (unsigned int j=0; j<vtx->getTracks().size(); ++j) {
      mom += vtx->getTracks()[j]->Vect();
      mom2 += vtx->getTracks()[j]->momentumAtVertex(vtx);
      _data.chrg += vtx->getTracks()[j]->getCharge();
    }

    _data.dirdot = mom.Unit().Dot( pos.Unit() );
    _data.dirdot2 = mom2.Unit().Dot( pos.Unit() );
    _data.ntrk = vtx->getTracks().size();

    // compute ks mass
    if (_data.ntrk == 2 && _data.chrg == 0) {
      const Track* trk1 = vtx->getTracks()[0];
      const Track* trk2 = vtx->getTracks()[1];
      TVector3 mom1 = trk1->Vect();
      TVector3 mom2 = trk2->Vect();
      TVector3 mom1v = trk1->momentumAtVertex(vtx);
      TVector3 mom2v = trk2->momentumAtVertex(vtx);

      TLorentzVector lvec1;
      TLorentzVector lvec2;
      lvec1.SetVectM( mom1, 0.1396 );
      lvec2.SetVectM( mom2, 0.1396 );
      _data.mks = (lvec1+lvec2).M();

      lvec1.SetVectM( mom1v, 0.1396 );
      lvec2.SetVectM( mom2v, 0.1396 );
      _data.mks2 = (lvec1+lvec2).M();

      // compute l0 mass
      TLorentzVector protonForLambda;
      TLorentzVector pionForLambda;
      if (mom1.Mag() > mom2.Mag()) {
        protonForLambda.SetVectM( mom1, 0.9383 );
        pionForLambda.SetVectM( mom2, 0.1396 );
      } else {
        protonForLambda.SetVectM( mom2, 0.9383 );
        pionForLambda.SetVectM( mom1, 0.1396 );
      }
      _data.ml0 = (protonForLambda+pionForLambda).M();

      if (mom1v.Mag() > mom2v.Mag()) {
        protonForLambda.SetVectM( mom1v, 0.9383 );
        pionForLambda.SetVectM( mom2v, 0.1396 );
      } else {
        protonForLambda.SetVectM( mom2v, 0.9383 );
        pionForLambda.SetVectM( mom1v, 0.1396 );
      }
      _data.ml02 = (protonForLambda+pionForLambda).M();

      // compute photon mass
      double ang1 = atan( trk1->getTanLambda() );
      double ang2 = atan( trk2->getTanLambda() );
      _data.mconv = sqrt( mom1.Mag()*mom2.Mag()*(1-cos(ang1-ang2)) );


      const MCParticle* mcp1 = trk1->getMcp();
      const MCParticle* mcp2 = trk2->getMcp();

      if (mcp1 && mcp2) {
        _data.mcpdg1 = mcp1->getPDG();
        _data.mcpdg2 = mcp2->getPDG();
        _data.mcppdg1 = mcp1->getParent()->getPDG();
        _data.mcppdg2 = mcp2->getParent()->getPDG();
      }

      if (mcp1 && mcp2) {
        const MCParticle* parent1 = mcp1->getParent();
        const MCParticle* parent2 = mcp2->getParent();
        const MCParticle* parent = mcp1->getParent();

        _data.mcpp1 = (int)((long long) parent1 );
        _data.mcpp2 = (int)((long long) parent2 );

        if ( abs(mcp1->getPDG())==11 && abs(mcp2->getPDG())==11 && parent->getPDG()==22 ) _data.conv = 1;
        if ( abs(mcp1->getPDG())==211 && abs(mcp2->getPDG())==211 && parent->getPDG()==310 ) _data.ks = 1;
        if ( abs(mcp1->getPDG())==211 && abs(mcp2->getPDG())==2212 && abs(parent->getPDG())==3122 ) _data.l0 = 1;
        if ( abs(mcp2->getPDG())==211 && abs(mcp1->getPDG())==2212 && abs(parent->getPDG())==3122 ) _data.l0 = 1;
      }

      _data.v0 = _data.ks || _data.l0 || _data.conv;
    }
    _ntp->Fill();
  }
}

void TestAlgoV0::end() {
  _file->Write();
  _file->Close();
}

void ZHHAlgo::init(Parameters* param) {
  Algorithm::init(param);

  string filename = param->get("FileName",string("test.root"));
  _jetname8 = param->get("JetCollectionName8",string("RefinedJets_8"));
  _jetname7 = param->get("JetCollectionName7",string("RefinedJets_7"));
  _jetname  = param->get("JetCollectionName6",string("RefinedJets_6"));
  _jetname5 = param->get("JetCollectionName5",string("RefinedJets_5"));
  _jetname4 = param->get("JetCollectionName4",string("RefinedJets_4"));

  _jetnamenv8 = param->get("JetCollectionNameNV8",string("RefinedNJets_8"));
  _jetnamenv7 = param->get("JetCollectionNameNV7",string("RefinedNJets_7"));
  _jetnamenv6 = param->get("JetCollectionNameNV6",string("RefinedNJets_6"));
  _jetnamenv5 = param->get("JetCollectionNameNV5",string("RefinedNJets_5"));
  _jetnamenv4 = param->get("JetCollectionNameNV4",string("RefinedNJets_4"));

  _file = new TFile(filename.c_str(),"RECREATE");
  _tree = new TTree("tree","tree");

  // mc info
  _tree->Branch("mchdecaypdg",&_d.mchdecaypdg,"mchdecaypdg[2]/I");
  _tree->Branch("mchbb",&_d.mchbb,"mchbb/I");
  _tree->Branch("mcnb",&_d.mcnb,"mcnb/I");

  // non-jet variables
  _tree->Branch("ycuts",&_d.ycuts,"ycuts[10]/D");

  _tree->Branch("thrust",&_d.thrust,"thrust/D");
  _tree->Branch("thaxis",&_d.thaxis,"thaxis[3]/D");
  _tree->Branch("ntr", &_d.ntr, "ntr/I");
  _tree->Branch("npfo", &_d.npfo, "ntr/I");

  // combined variables for compatibility
  _tree->Branch("mass",&_d.mass,"mass[15]/D");
  _tree->Branch("ntrjetmin",&_d.ntrjetmin,"ntrjetmin/D");
  _tree->Branch("pmiss",&_d.pmiss,"pmiss[3]/D");
  _tree->Branch("emiss",&_d.emiss,"emiss/D");

  // 6-jet variables
  _tree->Branch("bcat",&_d.bcat,"bcat[6]/D");
  _tree->Branch("btag",&_d.btag,"btag[6]/D");
  _tree->Branch("ctag",&_d.ctag,"ctag[6]/D");
  _tree->Branch("ejet",&_d.ejet,"ejet[6]/D");
  _tree->Branch("pxjet",&_d.pxjet,"pxjet[6]/D");
  _tree->Branch("pyjet",&_d.pyjet,"pyjet[6]/D");
  _tree->Branch("pzjet",&_d.pzjet,"pzjet[6]/D");
  _tree->Branch("ntrjet",&_d.ntrjet,"ntrjet[6]/D");
  _tree->Branch("mcnb6",&_d.mcnb6,"mcnb6[6]/D");
  _tree->Branch("mcnc6",&_d.mcnc6,"mcnc6[6]/D");
  _tree->Branch("twovtxprobjet",&_d.twovtxprobjet,"twovtxprobjet[6]/D");
  _tree->Branch("vtxangle",&_d.vtxangle,"vtxangle[6]/D");

  // 4-jet variables
  _tree->Branch("bcat4",&_d.bcat4,"bcat4[4]/D");
  _tree->Branch("btag4",&_d.btag4,"btag4[4]/D");
  _tree->Branch("ctag4",&_d.ctag4,"ctag4[4]/D");
  _tree->Branch("ejet4",&_d.ejet4,"ejet4[4]/D");
  _tree->Branch("pxjet4",&_d.pxjet4,"pxjet4[4]/D");
  _tree->Branch("pyjet4",&_d.pyjet4,"pyjet4[4]/D");
  _tree->Branch("pzjet4",&_d.pzjet4,"pzjet4[4]/D");
  _tree->Branch("ntrjet4",&_d.ntrjet4,"ntrjet4[4]/D");
  _tree->Branch("twovtxprobjet4",&_d.twovtxprobjet4,"twovtxprobjet4[4]/D");
  _tree->Branch("vtxangle4",&_d.vtxangle4,"vtxangle4[4]/D");

  _tree->Branch("bcat5",&_d.bcat5,"bcat5[5]/D");
  _tree->Branch("btag5",&_d.btag5,"btag5[5]/D");
  _tree->Branch("ctag5",&_d.ctag5,"ctag5[5]/D");
  _tree->Branch("ejet5",&_d.ejet5,"ejet5[5]/D");
  _tree->Branch("pxjet5",&_d.pxjet5,"pxjet5[5]/D");
  _tree->Branch("pyjet5",&_d.pyjet5,"pyjet5[5]/D");
  _tree->Branch("pzjet5",&_d.pzjet5,"pzjet5[5]/D");
  _tree->Branch("ntrjet5",&_d.ntrjet5,"ntrjet5[5]/D");
  _tree->Branch("twovtxprobjet5",&_d.twovtxprobjet5,"twovtxprobjet5[5]/D");
  _tree->Branch("vtxangle5",&_d.vtxangle5,"vtxangle5[5]/D");

  _tree->Branch("bcat7",&_d.bcat7,"bcat7[7]/D");
  _tree->Branch("btag7",&_d.btag7,"btag7[7]/D");
  _tree->Branch("ctag7",&_d.ctag7,"ctag7[7]/D");
  _tree->Branch("ejet7",&_d.ejet7,"ejet7[7]/D");
  _tree->Branch("pxjet7",&_d.pxjet7,"pxjet7[7]/D");
  _tree->Branch("pyjet7",&_d.pyjet7,"pyjet7[7]/D");
  _tree->Branch("pzjet7",&_d.pzjet7,"pzjet7[7]/D");
  _tree->Branch("ntrjet7",&_d.ntrjet7,"ntrjet7[7]/D");
  _tree->Branch("twovtxprobjet7",&_d.twovtxprobjet7,"twovtxprobjet7[7]/D");
  _tree->Branch("vtxangle7",&_d.vtxangle7,"vtxangle7[7]/D");

  _tree->Branch("bcat8",&_d.bcat8,"bcat8[8]/D");
  _tree->Branch("btag8",&_d.btag8,"btag8[8]/D");
  _tree->Branch("ctag8",&_d.ctag8,"ctag8[8]/D");
  _tree->Branch("ejet8",&_d.ejet8,"ejet8[8]/D");
  _tree->Branch("pxjet8",&_d.pxjet8,"pxjet8[8]/D");
  _tree->Branch("pyjet8",&_d.pyjet8,"pyjet8[8]/D");
  _tree->Branch("pzjet8",&_d.pzjet8,"pzjet8[8]/D");
  _tree->Branch("ntrjet8",&_d.ntrjet8,"ntrjet8[8]/D");
  _tree->Branch("twovtxprobjet8",&_d.twovtxprobjet8,"twovtxprobjet8[8]/D");
  _tree->Branch("vtxangle8",&_d.vtxangle8,"vtxangle8[8]/D");

  // jet clustering with no vertex
  _tree->Branch("bcatnv4",&_d.bcatnv4,"bcatnv4[4]/D");
  _tree->Branch("btagnv4",&_d.btagnv4,"btagnv4[4]/D");
  _tree->Branch("ctagnv4",&_d.ctagnv4,"ctagnv4[4]/D");
  _tree->Branch("ejetnv4",&_d.ejetnv4,"ejetnv4[4]/D");
  _tree->Branch("pxjetnv4",&_d.pxjetnv4,"pxjetnv4[4]/D");
  _tree->Branch("pyjetnv4",&_d.pyjetnv4,"pyjetnv4[4]/D");
  _tree->Branch("pzjetnv4",&_d.pzjetnv4,"pzjetnv4[4]/D");
  _tree->Branch("ntrjetnv4",&_d.ntrjetnv4,"ntrjetnv4[4]/D");
  _tree->Branch("twovtxprobjetnv4",&_d.twovtxprobjetnv4,"twovtxprobjetnv4[4]/D");
  _tree->Branch("vtxanglenv4",&_d.vtxanglenv4,"vtxanglenv4[4]/D");

  _tree->Branch("bcatnv5",&_d.bcatnv5,"bcatnv5[5]/D");
  _tree->Branch("btagnv5",&_d.btagnv5,"btagnv5[5]/D");
  _tree->Branch("ctagnv5",&_d.ctagnv5,"ctagnv5[5]/D");
  _tree->Branch("ejetnv5",&_d.ejetnv5,"ejetnv5[5]/D");
  _tree->Branch("pxjetnv5",&_d.pxjetnv5,"pxjetnv5[5]/D");
  _tree->Branch("pyjetnv5",&_d.pyjetnv5,"pyjetnv5[5]/D");
  _tree->Branch("pzjetnv5",&_d.pzjetnv5,"pzjetnv5[5]/D");
  _tree->Branch("ntrjetnv5",&_d.ntrjetnv5,"ntrjetnv5[5]/D");
  _tree->Branch("twovtxprobjetnv5",&_d.twovtxprobjetnv5,"twovtxprobjetnv5[5]/D");
  _tree->Branch("vtxanglenv5",&_d.vtxanglenv5,"vtxanglenv5[5]/D");

  _tree->Branch("bcatnv6",&_d.bcatnv6,"bcatnv6[6]/D");
  _tree->Branch("btagnv6",&_d.btagnv6,"btagnv6[6]/D");
  _tree->Branch("ctagnv6",&_d.ctagnv6,"ctagnv6[6]/D");
  _tree->Branch("ejetnv6",&_d.ejetnv6,"ejetnv6[6]/D");
  _tree->Branch("pxjetnv6",&_d.pxjetnv6,"pxjetnv6[6]/D");
  _tree->Branch("pyjetnv6",&_d.pyjetnv6,"pyjetnv6[6]/D");
  _tree->Branch("pzjetnv6",&_d.pzjetnv6,"pzjetnv6[6]/D");
  _tree->Branch("ntrjetnv6",&_d.ntrjetnv6,"ntrjetnv6[6]/D");
  _tree->Branch("twovtxprobjetnv6",&_d.twovtxprobjetnv6,"twovtxprobjetnv6[6]/D");
  _tree->Branch("vtxanglenv6",&_d.vtxanglenv6,"vtxanglenv6[6]/D");

  _tree->Branch("bcatnv7",&_d.bcatnv7,"bcatnv7[7]/D");
  _tree->Branch("btagnv7",&_d.btagnv7,"btagnv7[7]/D");
  _tree->Branch("ctagnv7",&_d.ctagnv7,"ctagnv7[7]/D");
  _tree->Branch("ejetnv7",&_d.ejetnv7,"ejetnv7[7]/D");
  _tree->Branch("pxjetnv7",&_d.pxjetnv7,"pxjetnv7[7]/D");
  _tree->Branch("pyjetnv7",&_d.pyjetnv7,"pyjetnv7[7]/D");
  _tree->Branch("pzjetnv7",&_d.pzjetnv7,"pzjetnv7[7]/D");
  _tree->Branch("ntrjetnv7",&_d.ntrjetnv7,"ntrjetnv7[7]/D");
  _tree->Branch("twovtxprobjetnv7",&_d.twovtxprobjetnv7,"twovtxprobjetnv7[7]/D");
  _tree->Branch("vtxanglenv7",&_d.vtxanglenv7,"vtxanglenv7[7]/D");

  _tree->Branch("bcatnv8",&_d.bcatnv8,"bcatnv8[8]/D");
  _tree->Branch("btagnv8",&_d.btagnv8,"btagnv8[8]/D");
  _tree->Branch("ctagnv8",&_d.ctagnv8,"ctagnv8[8]/D");
  _tree->Branch("ejetnv8",&_d.ejetnv8,"ejetnv8[8]/D");
  _tree->Branch("pxjetnv8",&_d.pxjetnv8,"pxjetnv8[8]/D");
  _tree->Branch("pyjetnv8",&_d.pyjetnv8,"pyjetnv8[8]/D");
  _tree->Branch("pzjetnv8",&_d.pzjetnv8,"pzjetnv8[8]/D");
  _tree->Branch("ntrjetnv8",&_d.ntrjetnv8,"ntrjetnv8[8]/D");
  _tree->Branch("twovtxprobjetnv8",&_d.twovtxprobjetnv8,"twovtxprobjetnv8[8]/D");
  _tree->Branch("vtxanglenv8",&_d.vtxanglenv8,"vtxanglenv8[8]/D");

  _jets = 0;
  _jets4 = 0;
  _jets5 = 0;
  _jets7 = 0;
  _jets8 = 0;

  _jetsnv4 = 0;
  _jetsnv5 = 0;
  _jetsnv6 = 0;
  _jetsnv7 = 0;
  _jetsnv8 = 0;
}



void ZHHAlgo::process() {
  if (!_jets4) {
    Event::Instance()->Get(_jetname4.c_str(), _jets4);
  }
  if (!_jets5) {
    Event::Instance()->Get(_jetname5.c_str(), _jets5);
  }
  if (!_jets) {
    Event::Instance()->Get(_jetname.c_str(), _jets);
  }
  if (!_jets7) {
    Event::Instance()->Get(_jetname7.c_str(), _jets7);
  }
  if (!_jets8) {
    Event::Instance()->Get(_jetname8.c_str(), _jets8);
  }
  if (!_jetsnv4) {
    Event::Instance()->Get(_jetnamenv4.c_str(), _jetsnv4);
  }
  if (!_jetsnv5) {
    Event::Instance()->Get(_jetnamenv5.c_str(), _jetsnv5);
  }
  if (!_jetsnv6) {
    Event::Instance()->Get(_jetnamenv6.c_str(), _jetsnv6);
  }
  if (!_jetsnv7) {
    Event::Instance()->Get(_jetnamenv7.c_str(), _jetsnv7);
  }
  if (!_jetsnv8) {
    Event::Instance()->Get(_jetnamenv8.c_str(), _jetsnv8);
  }

  const Vertex* privtx = Event::Instance()->getPrimaryVertex();

  // check higgs decay & nbs
  const MCParticleVec& mcps = Event::Instance()->getMCParticles();

  _d.mcnb = 0;
  _d.mchdecaypdg[0] = _d.mchdecaypdg[1] = 0;
  _d.mchbb = 0;

  int hcount = 0;
  for (unsigned int i=0; i<mcps.size(); i++) {
    int abspdg = abs(mcps[i]->getPDG());
    int parpdg = 0;
    if (mcps[i]->getParent())parpdg = abs(mcps[i]->getParent()->getPDG());
    if (((abspdg > 500 && abspdg < 600) || (abspdg > 5000 && abspdg < 6000)) && parpdg < 100)
      _d.mcnb ++;

    if (mcps[i]->getPDG() == 25) {
      // higgs
      if (mcps[i]->getDaughters().size() != 2) {
        cout << "ERR: # of higgs daughters = " << mcps[i]->getDaughters().size() << endl;
        break;
      }
      if (hcount == 2) {
        cout << "Too many higgs found!, ignore decay" << endl;
        break;
      }
      int apdg = abs(mcps[i]->getDaughters()[0]->getPDG());
      _d.mchdecaypdg[hcount++] = apdg;
      if (apdg == 4)_d.mchbb ++;
    }
    if (mcps[i]->getPDG() == 5 && parpdg == 21 && mcps[i]->getParent()->getDaughters().size() == 2) {
      TVector3 v1 = mcps[i]->getParent()->getDaughters()[0]->Vect();
      TVector3 v2 = mcps[i]->getParent()->getDaughters()[1]->Vect();
      cout << "g->bb angle: " << v1.Angle(v2) << endl;
    }
  }

  // thrust
  vector<TVector3> v;
  const TrackVec& tracks = Event::Instance()->getTracks();
  const NeutralVec& neutrals = Event::Instance()->getNeutrals();

  for (unsigned int n=0; n<tracks.size(); n++) {
    v.push_back(tracks[n]->Vect());
  }
  for (unsigned int n=0; n<neutrals.size(); n++) {
    v.push_back(neutrals[n]->Vect());
  }

  TVector3 taxis;
  _d.thrust = algoEtc::calcThrust(v, taxis);
  _d.thaxis[0] = taxis.x();
  _d.thaxis[1] = taxis.y();
  _d.thaxis[2] = taxis.z();

  _d.ntr = tracks.size();
  _d.npfo = tracks.size() + neutrals.size();

  // sorting btag
  vector<const Jet*> jets;
  jets = *_jets;
  sort(jets.begin(), jets.end(), sortBtag);

  int nmass = 0;
  TLorentzVector totp;
  _d.ntrjetmin = 10000.;

  for (unsigned int nj = 0; nj < 6; nj ++) {
    Jet* j = const_cast<Jet*>(jets[nj]);  // TODO: bad boy...
    j->recalcFourMomentum();

    _d.btag[nj] = j->getParam("lcfiplus")->get<double>("BTag");
    _d.bcat[nj] = j->getParam("lcfiplus")->get<double>("Category");
    _d.ctag[nj] = j->getParam("lcfiplus")->get<double>("CTag");
    _d.ejet[nj] = j->E();
    _d.pxjet[nj] = j->Px();
    _d.pyjet[nj] = j->Py();
    _d.pzjet[nj] = j->Pz();

    _d.mcnb6[nj] = j->getParam("lcfiplus")->get<double>("MCnb");
    _d.mcnc6[nj] = j->getParam("lcfiplus")->get<double>("MCnc");

    totp += *j;

    unsigned int ntr = j->getAllTracks().size();
    _d.ntrjet[nj] = ntr;
    if (_d.ntrjetmin > ntr)_d.ntrjetmin = ntr;

    unsigned int nv = j->getVertices().size();
    _d.twovtxprobjet[nj] = 1.;
    _d.vtxangle[nj] = 0.;
    if (nv == 2) {
      // two vtx prob
      _d.twovtxprobjet[nj] *= j->getVertices()[0]->getProb();
      _d.twovtxprobjet[nj] *= j->getVertices()[1]->getProb();

      // vertex angle
      TVector3 pripos = privtx->getPos();
      TVector3 pos1 = j->getVertices()[0]->getPos() - pripos;
      TVector3 pos2 = j->getVertices()[1]->getPos() - pripos;

      _d.vtxangle[nj] = pos1.Angle(pos2);
    }

    // ycut values
    if (nj == 0) {
      for (int i=0; i<10; i++) {
        TString s;
        s.Form("y%d%d",i,i+1);
        _d.ycuts[i] = j->getParam("yth")->get<double>(s);
      }
    }

    // masses
    if (nj == 5)continue;
    for (unsigned int nj2 = nj + 1; nj2 < 6; nj2 ++) {
      const Jet* j2 = jets[nj2];
      TLorentzVector v = *j;
      v += *j2;
      _d.mass[nmass++] = v.M();
    }

  }
  _d.emiss = 500 - totp.E();
  _d.pmiss[0] = -totp.Px();
  _d.pmiss[1] = -totp.Py();
  _d.pmiss[2] = -totp.Pz();

  // sorting btag for
  vector<const Jet*> jets4;
  jets4 = *_jets4;
  sort(jets4.begin(), jets4.end(), sortBtag);

  for (unsigned int nj = 0; nj < 4; nj ++) {
    Jet* j = const_cast<Jet*>(jets4[nj]);  // TODO: bad boy...
    j->recalcFourMomentum();

    _d.btag4[nj] = j->getParam("lcfiplus")->get<double>("BTag");
    _d.bcat4[nj] = j->getParam("lcfiplus")->get<double>("Category");
    _d.ctag4[nj] = j->getParam("lcfiplus")->get<double>("CTag");
    _d.ejet4[nj] = j->E();
    _d.pxjet4[nj] = j->Px();
    _d.pyjet4[nj] = j->Py();
    _d.pzjet4[nj] = j->Pz();

    unsigned int ntr = j->getAllTracks().size();
    _d.ntrjet4[nj] = ntr;

    unsigned int nv = j->getVertices().size();
    _d.twovtxprobjet4[nj] = 1.;
    _d.vtxangle4[nj] = 0.;
    if (nv == 2) {
      // two vtx prob
      _d.twovtxprobjet4[nj] *= j->getVertices()[0]->getProb();
      _d.twovtxprobjet4[nj] *= j->getVertices()[1]->getProb();

      // vertex angle
      TVector3 pripos = privtx->getPos();
      TVector3 pos1 = j->getVertices()[0]->getPos() - pripos;
      TVector3 pos2 = j->getVertices()[1]->getPos() - pripos;

      _d.vtxangle4[nj] = pos1.Angle(pos2);
    }
  }

  // sorting btag for
  vector<const Jet*> jets5;
  jets5 = *_jets5;
  sort(jets5.begin(), jets5.end(), sortBtag);

  for (unsigned int nj = 0; nj < 5; nj ++) {
    Jet* j = const_cast<Jet*>(jets5[nj]);  // TODO: bad boy...
    j->recalcFourMomentum();

    _d.btag5[nj] = j->getParam("lcfiplus")->get<double>("BTag");
    _d.bcat5[nj] = j->getParam("lcfiplus")->get<double>("Category");
    _d.ctag5[nj] = j->getParam("lcfiplus")->get<double>("CTag");
    _d.ejet5[nj] = j->E();
    _d.pxjet5[nj] = j->Px();
    _d.pyjet5[nj] = j->Py();
    _d.pzjet5[nj] = j->Pz();

    unsigned int ntr = j->getAllTracks().size();
    _d.ntrjet5[nj] = ntr;

    unsigned int nv = j->getVertices().size();
    _d.twovtxprobjet5[nj] = 1.;
    _d.vtxangle5[nj] = 0.;
    if (nv == 2) {
      // two vtx prob
      _d.twovtxprobjet5[nj] *= j->getVertices()[0]->getProb();
      _d.twovtxprobjet5[nj] *= j->getVertices()[1]->getProb();

      // vertex angle
      TVector3 pripos = privtx->getPos();
      TVector3 pos1 = j->getVertices()[0]->getPos() - pripos;
      TVector3 pos2 = j->getVertices()[1]->getPos() - pripos;

      _d.vtxangle5[nj] = pos1.Angle(pos2);
    }
  }

  // sorting btag for
  vector<const Jet*> jets7;
  jets7 = *_jets7;
  sort(jets7.begin(), jets7.end(), sortBtag);

  for (unsigned int nj = 0; nj < 7; nj ++) {
    Jet* j = const_cast<Jet*>(jets7[nj]);  // TODO: bad boy...
    j->recalcFourMomentum();

    _d.btag7[nj] = j->getParam("lcfiplus")->get<double>("BTag");
    _d.bcat7[nj] = j->getParam("lcfiplus")->get<double>("Category");
    _d.ctag7[nj] = j->getParam("lcfiplus")->get<double>("CTag");
    _d.ejet7[nj] = j->E();
    _d.pxjet7[nj] = j->Px();
    _d.pyjet7[nj] = j->Py();
    _d.pzjet7[nj] = j->Pz();

    unsigned int ntr = j->getAllTracks().size();
    _d.ntrjet7[nj] = ntr;

    unsigned int nv = j->getVertices().size();
    _d.twovtxprobjet7[nj] = 1.;
    _d.vtxangle7[nj] = 0.;
    if (nv == 2) {
      // two vtx prob
      _d.twovtxprobjet7[nj] *= j->getVertices()[0]->getProb();
      _d.twovtxprobjet7[nj] *= j->getVertices()[1]->getProb();

      // vertex angle
      TVector3 pripos = privtx->getPos();
      TVector3 pos1 = j->getVertices()[0]->getPos() - pripos;
      TVector3 pos2 = j->getVertices()[1]->getPos() - pripos;

      _d.vtxangle7[nj] = pos1.Angle(pos2);
    }
  }

  // sorting btag for
  vector<const Jet*> jets8;
  jets8 = *_jets8;
  sort(jets8.begin(), jets8.end(), sortBtag);

  for (unsigned int nj = 0; nj < 8; nj ++) {
    Jet* j = const_cast<Jet*>(jets8[nj]);  // TODO: bad boy...
    j->recalcFourMomentum();

    _d.btag8[nj] = j->getParam("lcfiplus")->get<double>("BTag");
    _d.bcat8[nj] = j->getParam("lcfiplus")->get<double>("Category");
    _d.ctag8[nj] = j->getParam("lcfiplus")->get<double>("CTag");
    _d.ejet8[nj] = j->E();
    _d.pxjet8[nj] = j->Px();
    _d.pyjet8[nj] = j->Py();
    _d.pzjet8[nj] = j->Pz();

    unsigned int ntr = j->getAllTracks().size();
    _d.ntrjet8[nj] = ntr;

    unsigned int nv = j->getVertices().size();
    _d.twovtxprobjet8[nj] = 1.;
    _d.vtxangle8[nj] = 0.;
    if (nv == 2) {
      // two vtx prob
      _d.twovtxprobjet8[nj] *= j->getVertices()[0]->getProb();
      _d.twovtxprobjet8[nj] *= j->getVertices()[1]->getProb();

      // vertex angle
      TVector3 pripos = privtx->getPos();
      TVector3 pos1 = j->getVertices()[0]->getPos() - pripos;
      TVector3 pos2 = j->getVertices()[1]->getPos() - pripos;

      _d.vtxangle8[nj] = pos1.Angle(pos2);
    }
  }

  vector<const Jet*> jetsnv4;
  jetsnv4 = *_jetsnv4;
  sort(jetsnv4.begin(), jetsnv4.end(), sortBtag);
  for (unsigned int nj = 0; nj < 4; nj ++) {
    Jet* j = const_cast<Jet*>(jetsnv4[nj]);  // TODO: bad boy...
    j->recalcFourMomentum();

    _d.btagnv4[nj] = j->getParam("lcfiplus")->get<double>("BTag");
    _d.bcatnv4[nj] = j->getParam("lcfiplus")->get<double>("Category");
    _d.ctagnv4[nj] = j->getParam("lcfiplus")->get<double>("CTag");
    _d.ejetnv4[nj] = j->E();
    _d.pxjetnv4[nj] = j->Px();
    _d.pyjetnv4[nj] = j->Py();
    _d.pzjetnv4[nj] = j->Pz();

    unsigned int ntr = j->getAllTracks().size();
    _d.ntrjetnv4[nj] = ntr;

    unsigned int nv = j->getVertices().size();
    _d.twovtxprobjetnv4[nj] = 1.;
    _d.vtxanglenv4[nj] = 0.;
    if (nv == 2) {
      // two vtx prob
      _d.twovtxprobjetnv4[nj] *= j->getVertices()[0]->getProb();
      _d.twovtxprobjetnv4[nj] *= j->getVertices()[1]->getProb();

      // vertex angle
      TVector3 pripos = privtx->getPos();
      TVector3 pos1 = j->getVertices()[0]->getPos() - pripos;
      TVector3 pos2 = j->getVertices()[1]->getPos() - pripos;

      _d.vtxanglenv4[nj] = pos1.Angle(pos2);
    }
  }

  vector<const Jet*> jetsnv5;
  jetsnv5 = *_jetsnv5;
  sort(jetsnv5.begin(), jetsnv5.end(), sortBtag);
  for (unsigned int nj = 0; nj < 5; nj ++) {
    Jet* j = const_cast<Jet*>(jetsnv5[nj]);  // TODO: bad boy...
    j->recalcFourMomentum();

    _d.btagnv5[nj] = j->getParam("lcfiplus")->get<double>("BTag");
    _d.bcatnv5[nj] = j->getParam("lcfiplus")->get<double>("Category");
    _d.ctagnv5[nj] = j->getParam("lcfiplus")->get<double>("CTag");
    _d.ejetnv5[nj] = j->E();
    _d.pxjetnv5[nj] = j->Px();
    _d.pyjetnv5[nj] = j->Py();
    _d.pzjetnv5[nj] = j->Pz();

    unsigned int ntr = j->getAllTracks().size();
    _d.ntrjetnv5[nj] = ntr;

    unsigned int nv = j->getVertices().size();
    _d.twovtxprobjetnv5[nj] = 1.;
    _d.vtxanglenv5[nj] = 0.;
    if (nv == 2) {
      // two vtx prob
      _d.twovtxprobjetnv5[nj] *= j->getVertices()[0]->getProb();
      _d.twovtxprobjetnv5[nj] *= j->getVertices()[1]->getProb();

      // vertex angle
      TVector3 pripos = privtx->getPos();
      TVector3 pos1 = j->getVertices()[0]->getPos() - pripos;
      TVector3 pos2 = j->getVertices()[1]->getPos() - pripos;

      _d.vtxanglenv5[nj] = pos1.Angle(pos2);
    }
  }

  vector<const Jet*> jetsnv6;
  jetsnv6 = *_jetsnv6;
  sort(jetsnv6.begin(), jetsnv6.end(), sortBtag);
  for (unsigned int nj = 0; nj < 6; nj ++) {
    Jet* j = const_cast<Jet*>(jetsnv6[nj]);  // TODO: bad boy...
    j->recalcFourMomentum();

    _d.btagnv6[nj] = j->getParam("lcfiplus")->get<double>("BTag");
    _d.bcatnv6[nj] = j->getParam("lcfiplus")->get<double>("Category");
    _d.ctagnv6[nj] = j->getParam("lcfiplus")->get<double>("CTag");
    _d.ejetnv6[nj] = j->E();
    _d.pxjetnv6[nj] = j->Px();
    _d.pyjetnv6[nj] = j->Py();
    _d.pzjetnv6[nj] = j->Pz();

    unsigned int ntr = j->getAllTracks().size();
    _d.ntrjetnv6[nj] = ntr;

    unsigned int nv = j->getVertices().size();
    _d.twovtxprobjetnv6[nj] = 1.;
    _d.vtxanglenv6[nj] = 0.;
    if (nv == 2) {
      // two vtx prob
      _d.twovtxprobjetnv6[nj] *= j->getVertices()[0]->getProb();
      _d.twovtxprobjetnv6[nj] *= j->getVertices()[1]->getProb();

      // vertex angle
      TVector3 pripos = privtx->getPos();
      TVector3 pos1 = j->getVertices()[0]->getPos() - pripos;
      TVector3 pos2 = j->getVertices()[1]->getPos() - pripos;

      _d.vtxanglenv6[nj] = pos1.Angle(pos2);
    }
  }

  vector<const Jet*> jetsnv7;
  jetsnv7 = *_jetsnv7;
  sort(jetsnv7.begin(), jetsnv7.end(), sortBtag);
  for (unsigned int nj = 0; nj < 7; nj ++) {
    Jet* j = const_cast<Jet*>(jetsnv7[nj]);  // TODO: bad boy...
    j->recalcFourMomentum();

    _d.btagnv7[nj] = j->getParam("lcfiplus")->get<double>("BTag");
    _d.bcatnv7[nj] = j->getParam("lcfiplus")->get<double>("Category");
    _d.ctagnv7[nj] = j->getParam("lcfiplus")->get<double>("CTag");
    _d.ejetnv7[nj] = j->E();
    _d.pxjetnv7[nj] = j->Px();
    _d.pyjetnv7[nj] = j->Py();
    _d.pzjetnv7[nj] = j->Pz();

    unsigned int ntr = j->getAllTracks().size();
    _d.ntrjetnv7[nj] = ntr;

    unsigned int nv = j->getVertices().size();
    _d.twovtxprobjetnv7[nj] = 1.;
    _d.vtxanglenv7[nj] = 0.;
    if (nv == 2) {
      // two vtx prob
      _d.twovtxprobjetnv7[nj] *= j->getVertices()[0]->getProb();
      _d.twovtxprobjetnv7[nj] *= j->getVertices()[1]->getProb();

      // vertex angle
      TVector3 pripos = privtx->getPos();
      TVector3 pos1 = j->getVertices()[0]->getPos() - pripos;
      TVector3 pos2 = j->getVertices()[1]->getPos() - pripos;

      _d.vtxanglenv7[nj] = pos1.Angle(pos2);
    }
  }

  vector<const Jet*> jetsnv8;
  jetsnv8 = *_jetsnv8;
  sort(jetsnv8.begin(), jetsnv8.end(), sortBtag);
  for (unsigned int nj = 0; nj < 8; nj ++) {
    Jet* j = const_cast<Jet*>(jetsnv8[nj]);  // TODO: bad boy...
    j->recalcFourMomentum();

    _d.btagnv8[nj] = j->getParam("lcfiplus")->get<double>("BTag");
    _d.bcatnv8[nj] = j->getParam("lcfiplus")->get<double>("Category");
    _d.ctagnv8[nj] = j->getParam("lcfiplus")->get<double>("CTag");
    _d.ejetnv8[nj] = j->E();
    _d.pxjetnv8[nj] = j->Px();
    _d.pyjetnv8[nj] = j->Py();
    _d.pzjetnv8[nj] = j->Pz();

    unsigned int ntr = j->getAllTracks().size();
    _d.ntrjetnv8[nj] = ntr;

    unsigned int nv = j->getVertices().size();
    _d.twovtxprobjetnv8[nj] = 1.;
    _d.vtxanglenv8[nj] = 0.;
    if (nv == 2) {
      // two vtx prob
      _d.twovtxprobjetnv8[nj] *= j->getVertices()[0]->getProb();
      _d.twovtxprobjetnv8[nj] *= j->getVertices()[1]->getProb();

      // vertex angle
      TVector3 pripos = privtx->getPos();
      TVector3 pos1 = j->getVertices()[0]->getPos() - pripos;
      TVector3 pos2 = j->getVertices()[1]->getPos() - pripos;

      _d.vtxanglenv8[nj] = pos1.Angle(pos2);
    }
  }

  _tree->Fill();
}

void ZHHAlgo::end() {
  _file->Write();
  _file->Close();
}

void TestAlgo::init(Parameters* param) {
  Algorithm::init(param);

  string filename = param->get("FileName",string("test.root"));
  _jetname = param->get("JetCollectionName",string("Durham_2Jets"));
  string primvtxcolname = param->get("PrimaryVertexCollectionName",string("PrimaryVertex"));
  string secvtxcolname = param->get("SecondaryVertexCollectionName",string("BuildUpVertex"));
  Event::Instance()->setDefaultPrimaryVertex(primvtxcolname.c_str());
  Event::Instance()->setDefaultSecondaryVertices(secvtxcolname.c_str());

  _file = new TFile(filename.c_str(),"RECREATE");
//		_nt = new TNtupleD("nt","nt","nev:pdg:d0:d0sig:z0:z0sig:e:pt:pz:chi2:sd0:sz0:ecaldep:hcaldep");
  _nt = new TNtupleD("nt","nt","nev:pdg:parpdg:d0:d0sig:z0:z0sig:e:pt:pz:chi2:invtx:minchi2:nvtx");

  _jets = 0;
  _nev = 0;
}

void TestAlgo::process() {

  if (!_jets) {
    Event::Instance()->Get(_jetname.c_str(), _jets);
  }
//		TrackVec &tracks = Event::Instance()->getTracks();
  const Vertex* privtx = Event::Instance()->getPrimaryVertex();
//		VertexVec &vtcs = Event::Instance()->getSecondaryVertices();

  for (unsigned int nj = 0; nj < _jets->size(); nj++) {
    const Jet* j = (*_jets)[nj];
    TrackVec tracks = j->getAllTracks(true);
    VertexVec& vtcs = j->getVertices();

    int nvtx = 0;
    for (unsigned int nv=0; nv<vtcs.size(); nv++)
      if (vtcs[nv]->getTracks().size() >=2) nvtx ++;

    for (unsigned int n=0; n<tracks.size(); n++) {
      const Track* tr = tracks[n];

      //double sd0 = signedD0(tr, j, privtx, true);
      //double sd0sig = signedD0Significance(tr, j, privtx, true);
      //double sz0 = signedZ0(tr, j, privtx, true);
      //double sz0sig = signedZ0Significance(tr, j, privtx, true);

      // vertex-track association
      int invtx = 0;
      double minchi2 = 1e+300;

      if (find(privtx->getTracks().begin(), privtx->getTracks().end(), tr) != privtx->getTracks().end())invtx = 1;
      else {
        for (unsigned int nv = 0; nv < vtcs.size(); nv ++) {
          double chi2 = 1e+300;
          if (find(vtcs[nv]->getTracks().begin(), vtcs[nv]->getTracks().end(), tr) != vtcs[nv]->getTracks().end()) {
            invtx = 2;
            chi2 = vtcs[nv]->getChi2Track(tr);
          } else { // if(vtcs[nv]->getTracks().size() >= 2){
            chi2 = vtcs[nv]->getChi2TrackFit(tr,3);
          }

          if (minchi2 > chi2)minchi2 = chi2;
        }
      }

      const MCParticle* mcp = tracks[n]->getMcp();
      const MCParticle* pmcp = (mcp ? mcp->getSemiStableBParent() : 0);
      _nt->Fill(_nev, mcp ? mcp->getPDG() : 0, pmcp ? pmcp->getPDG() : 0, fabs(tr->getD0()), fabs(tr->getD0() / sqrt(tr->getCovMatrix()[tpar::d0d0])),
                fabs(tr->getZ0()), fabs(tr->getZ0() / sqrt(tr->getCovMatrix()[tpar::z0z0])),
                tr->E(), tr->Pt(), tr->Pz(), tr->getChi2(), invtx, minchi2, nvtx * 10 + vtcs.size());
      //sd0, sz0, tr->getCaloEdep()[tpar::ecal], tr->getCaloEdep()[tpar::hcal]);

    }
  }

  _nev ++;
}

void TestAlgo::end() {
  _file->Write();
  _file->Close();
}

void FlavtagReader::init(Parameters* param) {
  Algorithm::init(param);

  string filename = param->get("FileName",string("test.root"));
  _jetname = param->get("JetCollectionName",string("Durham_2Jets"));
  string primvtxcolname = param->get("PrimaryVertexCollectionName",string("PrimaryVertex"));
  Event::Instance()->setDefaultPrimaryVertex(primvtxcolname.c_str());

  _file = new TFile(filename.c_str(),"RECREATE");
  _nt = new TNtupleD("nt","nt","nev:nj:e:px:py:pz:btag:ctag:otag:bbtag:bctag:cctag");
  _ntev = new TNtupleD("ntev","ntev","nev:btag1:btag2:btag3:btag4:btag5:btag6:ctag1:ctag2:ctag3:ctag4:ctag5:ctag6");

  _jets = 0;
  _nev = 0;
}

void FlavtagReader::process() {
  if (!_jets) {
    Event::Instance()->Get(_jetname.c_str(), _jets);
  }
  //const Vertex * privtx = Event::Instance()->getPrimaryVertex();

  vector<double> btags, ctags;

  for (unsigned int nj = 0; nj < _jets->size(); nj++) {
    const Jet* j = (*_jets)[nj];

    const Parameters* para = j->getParam("lcfiplus");
    _nt->Fill(_nev, nj, j->E(), j->Px(), j->Py(), j->Pz(),
              para->get<double>("BTag"), para->get<double>("CTag"), para->get<double>("OTag"),  para->get<double>("BBTag"),  para->get<double>("CCTag"),  para->get<double>("BCTag"));
    btags.push_back(para->get<double>("BTag"));
    ctags.push_back(para->get<double>("CTag"));

    cout << "nvtx = " << para->get<double>("nvtx") << ", nvtxall = " << para->get<double>("nvtxall") << endl;
  }

  std::sort(btags.begin(), btags.end());
  std::sort(ctags.begin(), ctags.end());
  if (_jets->size() >= 6)
    _ntev->Fill(_nev,  btags[0],  btags[1],  btags[2],  btags[3],  btags[4],  btags[5],  ctags[0],  ctags[1],  ctags[2],  ctags[3],  ctags[4],  ctags[5]);

  _nev ++;
}

void FlavtagReader::end() {
  _file->Write();
  _file->Close();
}

#if 0

void TestAlgo::init(Parameters* param) {
  Algorithm::init(param);

  gStyle->SetPalette(1);
  _h = new TH2D("h","h",200,-2,2,200,-2,2);
  _he = new TH2D("he","he",200,-2,2,200,-2,2);
}

void TestAlgo::process() {
  // check bbbbbb (reject H->WW etc.)
  const MCParticleVec& mcps = Event::Instance()->getMCParticles();
  const MCColorSingletVec& mccss = Event::Instance()->getMCColorSinglets();
  cout << "# mccs = " << mccss.size() << endl;
  if (mccss.size() != 3) return;

  int nq[3];
  for (unsigned int i=0; i<3; i++) {
    nq[i] = mccss[i]->_initials.size();
  }
  cout << "# qs = " << nq[0] << " " << nq[1] << " " << nq[2] << endl;
  if (nq[0]!=2 ||nq[1]!=2 ||nq[2]!=2)return;

  for (unsigned int i=0; i<mcps.size(); i++) {
    const MCColorSinglet* mccs = mcps[i]->getColorSinglet(&mccss);
    const MCParticle* p1 = mccs->_initials[0];
    const MCParticle* p2 = mccs->_initials[1];

    TVector3 normal = p1->Vect().Cross(p2->Vect()).Unit();
    double ndp = (normal.Dot(mcps[i]->Vect().Unit()));
    TVector3 pplane = mcps[i]->Vect().Unit() - (normal * ndp);
    double nxp1 = p1->Vect().Unit().Dot(pplane);
    double nxp2 = p2->Vect().Unit().Dot(pplane);
    double nxp12 = fabs(p1->Vect().Unit().Dot(p2->Vect().Unit()));
    double nxp = (nxp1 - nxp2) / (1 - nxp12);

    _h->Fill(nxp, ndp);
    cout << nxp << " " << ndp << endl;
  }
}

void TestAlgo::end() {
  _h->Draw("colz");
  gPad->Update();
  gSystem->ProcessEvents();
}


void TestAlgo::init(Parameters* param) {
  Algorithm::init(param);

  string filename = param->get("FileName",string("test.root"));
  _privtxname = param->get("PrimaryVertexCollectionName",string("PrimaryVertex"));
  _vtxname = param->get("BuildUpVertexCollectionName",string("BuildUpVertex"));
  _v0vtxname = param->get("V0VertexCollectionName",string("BuildUpVertex_V0"));
  _jetname = param->get("JetCollectionName",string("Durham_2Jets"));
  _vtxsel = param->get("VertexSelection",(int)0);
  _refine = param->get("PerformRefining",(int)0);
  _bbhh = param->get("IsBBHH",int(0));

  _file = new TFile(filename.c_str(),"RECREATE");
  _ntJet = new TNtupleD("ntJet","ntJet","nvtx:1vtxprob:2vtxprob:cflt:ecvtx:ejet:vangle:vmass:esingle");

  _vertices = 0;
  _jets = 0;
}



void TestAlgo::process() {
  if (!_vertices) {
    Event::Instance()->Get(_vtxname.c_str(), _vertices);
    //cout << "vtx name: " << _vtxname << ", pointer = " << (unsigned int)_vertices << endl;
  }
  if (!_v0vertices) {
    Event::Instance()->Get(_v0vtxname.c_str(), _v0vertices);
    //cout << "vtx name: " << _vtxname << ", pointer = " << (unsigned int)_vertices << endl;
  }

  if (!_jets) {
    Event::Instance()->Get(_jetname.c_str(), _jets);
    //cout << "jet name: " << _jetname << ", pointer = " << (unsigned int)_jets << endl;
  }

  // check bbbbbb (reject H->WW etc.)
  const MCParticleVec& mcps = Event::Instance()->getMCParticles();
  if (_bbhh) {
    int hcount = 0;
    for (unsigned int i=0; i<mcps.size(); i++) {
      if (mcps[i]->getPDG() != 25)continue;

      // higgs
      if (mcps[i]->getDaughters().size() != 2) {
        cout << "ERR: # of higgs daughters = " << mcps[i]->getDaughters().size() << endl;
        break;
      }
      if (abs(mcps[i]->getDaughters()[0]->getPDG()) != 5)break;
      if (abs(mcps[i]->getDaughters()[1]->getPDG()) != 5)break;
      hcount ++;
    }
    if (hcount < 2)return;
  }

  // select vertices
  vector<const Track*> residualTracks;
  vector<const Vertex*> selectedVertices;
  const vector<const Vertex*>* pVertices;
  if (_vertices && _vtxsel) {
    VertexSelectorConfig vscfg;
    vscfg.rejectdist = true;
    vscfg.minpos = .3;
    vscfg.maxpos = 30.;
    vscfg.rejectk0 = true;
    vscfg.k0width = .01;

    selectedVertices = VertexSelector()(*_vertices, vscfg, residualTracks,false);
    pVertices = &selectedVertices;
  } else {
    pVertices = _vertices;
  }

  cout << "# jet = " << _jets->size() << endl;

  // copy vertices
  vector<Vertex*> vtcs2;
  const Vertex* ip = Event::Instance()->getPrimaryVertex(_privtxname.c_str());
  if (ip == 0)throw(Exception("IP not found!"));

  for (unsigned int v=0; v<pVertices->size(); v++) {
    if (!(*pVertices)[v]->isPrimary())
      vtcs2.push_back(new Vertex(*(*pVertices)[v]));
  }

  cout << "# sec vtx = " << vtcs2.size() << endl;

  vector<vector<Vertex*> > jetVertices;
  vector<vector<const Track*> > jetResidualTracks;
  // jet-vtx association
  lcfiplus::algoEtc::connectVerticesToJets(*_jets, vtcs2, jetVertices, jetResidualTracks,ip);

  VertexFinderSuehara::VertexFinderSueharaConfig cfg;

  for (unsigned int j=0; j<_jets->size(); j++) {

    // single track probability
    double singleprob = 0;
    double twoprob = 0;
    double cflt = 0;
    double ecvtx = 0;
    double vangle = 0;
    double vmass = 0;
    double esingle = 0;

    if (_refine) {
      vector<Vertex*> singleVtcs = VertexFinderSuehara::makeSingleTrackVertices(constVector(jetVertices[j]), jetResidualTracks[j], *_v0vertices, ip, cfg);

      if (jetVertices[j].size() + singleVtcs.size() >= 2) {
        cout << "Before recombination:" << endl;
        for (unsigned int k=0; k<jetVertices[j].size(); k++)
          jetVertices[j][k]->Print();
        for (unsigned int k=0; k<singleVtcs.size(); k++)
          singleVtcs[k]->Print();
      }

      VertexFinderSuehara::recombineVertices(jetVertices[j], singleVtcs);

      // v0 selection again
      VertexSelector()(jetVertices[j], cfg.v0selVertex);

      vector<const Track*> singletracklist;
      if (jetVertices[j].size() > 1) {
        twoprob = jetVertices[j][0]->getProb() * jetVertices[j][1]->getProb();
        singletracklist.resize(jetVertices[j][0]->getTracks().size() + jetVertices[j][1]->getTracks().size());
        std::copy(jetVertices[j][0]->getTracks().begin(), jetVertices[j][0]->getTracks().end(), singletracklist.begin());
        std::copy(jetVertices[j][1]->getTracks().begin(), jetVertices[j][1]->getTracks().end(), singletracklist.begin() + jetVertices[j][0]->getTracks().size());

        Vertex* single = VertexFitterSimple_V()(singletracklist.begin(), singletracklist.end());

        singleprob = single->getProb();

        cout << "twoprob: " << jetVertices[j][0]->getProb() << " " << jetVertices[j][1]->getProb() << " oneprob: " << singleprob << endl;

        delete single;

        cflt = (jetVertices[j][1]->getPos() - jetVertices[j][0]->getPos()).Mag();

        // looking for near vertex
        int nnear = (jetVertices[j][0]->getPos().Mag() > jetVertices[j][1]->getPos().Mag() ? 0 : 1);

        for (unsigned int ntr=0; ntr<jetVertices[j][nnear]->getTracks().size(); ntr++) {
          const Track* tr = jetVertices[j][nnear]->getTracks()[ntr];
          ecvtx += tr->E();
        }
        cout << "cflt = " << cflt << ", ecvtx = " << ecvtx << ", ejet = " << (*_jets)[j]->E() << endl;

        // single track investigation
        int idx = -1;
        if (jetVertices[j][0]->getTracks().size() == 1) idx = 0;
        if (jetVertices[j][1]->getTracks().size() == 1) idx = 1;

        if (idx >= 0) {
          const Track* tr = jetVertices[j][idx]->getTracks()[0];
          TVector3 vpos = jetVertices[j][idx]->getPos();
          vangle = vpos.Angle(tr->Vect());
          vmass = 2 * tr->E() * tr->E() * (1 - cos(vpos.Angle(tr->Vect())));
          esingle = tr->E();
        }

      } else if (jetVertices[j].size() == 1) {
        singleprob = jetVertices[j][0]->getProb();
      }
    }

    if (jetVertices[j].size() >= 2) {
      for (unsigned int k=0; k<jetVertices[j].size(); k++)
        jetVertices[j][k]->Print();
    }

    _file->cd();
    _ntJet->Fill((int)jetVertices[j].size(),singleprob, twoprob, cflt, ecvtx, (*_jets)[j]->E(), vangle, vmass, esingle);
  }
}

void TestAlgo::init(Parameters* param) {
  Algorithm::init(param);


  string filename = param->get("FileName",string("test.root"));
  _jetname = param->get("JetCollectionName",string("Durham_6Jets2"));
  _bbhh = param->get("IsBBHH",int(0));

  _file = new TFile(filename.c_str(),"RECREATE");

  _ntJet2 = new TNtuple("ntJet2","ntJet2","nev:njet:nbjetmc:nvtx:nvtxjet:ngoodvtx:fracgoodvtxtrack:ycut:nbjet:fracgoodtrack");
  _nbJet = new TNtuple("nbJet", "number of b tracks in eachjet", "nev:nb1:nb2:nb3:nb4:nb5:nb6:nb11:nb12:nb13:nb14:nb15:nb16");

  _jets = 0;
}

void TestAlgo::process() {
  if (!_jets) {
    Event::Instance()->Get(_jetname.c_str(), _jets);
    //cout << "jet name: " << _jetname << ", pointer = " << (unsigned int)_jets << endl;
  }
  JetVec& jets = *_jets;
  unsigned int nj = jets.size();

  Event* event = Event::Instance();
  MCParticleVec& mcps = event->getMCParticles();

  // check bbbbbb (reject H->WW etc.)
  if (_bbhh) {
    int hcount = 0;
    for (unsigned int i=0; i<mcps.size(); i++) {
      if (mcps[i]->getPDG() != 25)continue;

      // higgs
      if (mcps[i]->getDaughters().size() != 2) {
        cout << "ERR: # of higgs daughters = " << mcps[i]->getDaughters().size() << endl;
        break;
      }
      if (abs(mcps[i]->getDaughters()[0]->getPDG()) != 5)break;
      if (abs(mcps[i]->getDaughters()[1]->getPDG()) != 5)break;
      hcount ++;
    }
    if (hcount < 2)return;
  }

  // semistable B
  vector<const MCParticle*> blist;
  for (unsigned int i=0; i<mcps.size(); i++) {
    if (mcps[i]->isSemiStableB()) {
      blist.push_back(mcps[i]);
    }
  }
  cout << "Number of semistable B: " << blist.size() << endl;

//	TNtuple *ntResidual = new TNtuple("ntResidual","ResidualTracks","nev:bid:btracks:mcvx:mcvy:mcvz:d0:d0err:z0:z0err:tre");
  // calculate btracks
  /*		int *btracks = new int[blist.size()];
  		memset(btracks,0,sizeof(int)*blist.size());
  		for(unsigned int i=0;i<tracks.size();i++){
  			for(unsigned int k=0;k<blist.size();k++){
  				if(tracks[i]->getMcp()->isParent(blist[k]))btracks[k] ++;
  			}
  		}
  */
  vector<const Track*> assignedTracks;
  vector<const Track*> residualTracks;

  map<const Jet*, int > nbmap;
  int nvtxjet = 0;
  for (unsigned int nj=0; nj<jets.size(); nj++) {
    nbmap[jets[nj]] = 0;
    if (jets[nj]->getVertices().size()>0)nvtxjet ++;
  }

  // nbjet fill
  vector<int> nbjet0(max(int(jets.size()),6));
  vector<int> nbjet1(max(int(jets.size()),6));
  for (unsigned int nj=0; nj<jets.size(); nj++) {
    nbjet0[nj] = nbjet1[nj] = 0;

    for (unsigned int i=0; i<jets[nj]->getTracks().size(); i++) {
      const Track* tr = jets[nj]->getTracks()[i];
      if (tr->getMcp() == 0)continue;
      if (tr->getMcp()->getSemiStableParent() == 0)continue;
      int pdg = tr->getMcp()->getSemiStableParent()->getPDG();
      if ((abs(pdg)>400&&abs(pdg)<600) || (abs(pdg)>4000&&abs(pdg)<6000)) {
        for (unsigned int k=0; k<blist.size(); k++) {
          if (tr->getMcp()->isParent(blist[k])) {
            nbjet0[nj] ++;
            if (tr->E()>1.)nbjet1[nj] ++;
            break;
          }
        }
      }
    }
    for (unsigned int nv=0; nv<jets[nj]->getVertices().size(); nv++) {
      for (unsigned int i=0; i<jets[nj]->getVertices()[nv]->getTracks().size(); i++) {
        const Track* tr = jets[nj]->getVertices()[nv]->getTracks()[i];
        if (tr->getMcp() == 0)continue;
        if (tr->getMcp()->getSemiStableParent() == 0)continue;
        int pdg = tr->getMcp()->getSemiStableParent()->getPDG();
        if ((abs(pdg)>400&&abs(pdg)<600) || (abs(pdg)>4000&&abs(pdg)<6000)) {
          for (unsigned int k=0; k<blist.size(); k++) {
            if (tr->getMcp()->isParent(blist[k])) {
              nbjet0[nj] ++;
              if (tr->E()>1.)nbjet1[nj] ++;
              break;
            }
          }
        }
      }
    }
  }

  // sort by decending order
  sort(nbjet0.begin(), nbjet0.end(), greater<int>());
  sort(nbjet1.begin(), nbjet1.end(), greater<int>());
  _nbJet->Fill(0,nbjet0[0], nbjet0[1], nbjet0[2], nbjet0[3], nbjet0[4], nbjet0[5],
               nbjet1[0], nbjet1[1], nbjet1[2], nbjet1[3], nbjet0[4], nbjet0[5]);

  for (unsigned int ib=0; ib<blist.size(); ib++) {
    vector<const Track*> aTracks, rTracks;
//Jet * JetMCMatch(vector<Jet *> &jets, MCParticle *mcp, vector<Track *> &assignedTracks, vector<Track *> &residualTracks)
    const Jet* jet = JetMCMatch(jets, blist[ib], aTracks, rTracks);
    if (jet) {
      nbmap[jet] ++;
      assignedTracks.insert(assignedTracks.end(), aTracks.begin(), aTracks.end());
      residualTracks.insert(residualTracks.end(), rTracks.begin(), rTracks.end());
    }
  }

  int nbjet = 0;
  for (unsigned int nj=0; nj<jets.size(); nj++) {
    if (nbmap[jets[nj]] > 0)nbjet ++;
  }
//	TNtuple *ntJet2 = new TNtuple("ntJet2","ntJet2","nev:njet:nbjetmc:nvtx:nvtxjet:ngoodvtx:ngoodvtxtrack:ycut:nbjet:fracgoodtrack");
//	TNtuple *ntJet2 = new TNtuple("ntJet2","ntJet2","nev:njet:ycut:nbjet:fracgoodtrack");
  double fracgoodtrack = double(residualTracks.size()) / double(assignedTracks.size()+residualTracks.size());
  _ntJet2->Fill(0, jets.size(), blist.size(), 0, nvtxjet, 0, 0, 0., nbjet, fracgoodtrack);
}


void TestAlgo::end() {
  _file->Write();
  _file->Close();
}
#endif

void VertexAnalysis::init(Parameters* param) {
  Algorithm::init(param);
  _privtxname = param->get("PrimaryVertexCollectionName",string("PrimaryVertex"));

  string filename = param->get("VertexAnalysis.FileName",string("VertexAnalysis.root"));
  _secvtxname = param->get("VertexAnalysis.SecondaryVertexCollectionName",string("BuildUpVertex"));
  _file = new TFile(filename.c_str(),"RECREATE");
  _nt = new TNtupleD("vtxtree","vtxtree","track:invtx:d0:d0err:z0:z0err:e:pt:chi2:ndf:vtxftdhits");
}

void VertexAnalysis::process() {
  const Vertex* privtx = Event::Instance()->getPrimaryVertex(_privtxname.c_str());
  TrackVec& tracks = Event::Instance()->getTracks();
  MCParticleVec& mcps = Event::Instance()->getMCParticles();
  VertexVec* psecvtx = 0;
  if (_secvtxname != "")
    Event::Instance()->Get(_secvtxname.c_str(), psecvtx);

  if (mcps.size() == 0) {
    cout << "MCParticle collection not specified. We need the MCParticle collection for vertex analysis." << endl;
    return;
  }

  for (unsigned int i=0; i < tracks.size(); ++i) {
    const Track* tr = tracks[i];
    const MCParticle* mcp = tr->getMcp();

    double trackseed, invtx, d0, d0err, z0, z0err, e, pt, chi2, ndf, vtxftdhits;

    if (mcp->getSemiStableParent() == 0) trackseed = 0.;
    else if (mcp->getSemiStableBParent() != 0) trackseed = 1.;
    else if (mcp->getSemiStableCParent() != 0) trackseed = 2.;
    else trackseed = 3.;

    invtx = (find(privtx->getTracks().begin(), privtx->getTracks().end(), tr) != privtx->getTracks().end());
    if (invtx == 0. && psecvtx) {
      // looking for secondary vertices
      for (unsigned int j=0; j<psecvtx->size(); j++) {
        TrackVec& vtr = (*psecvtx)[j]->getTracks();
        if (find(vtr.begin(), vtr.end(), tr) != vtr.end())
          invtx = 2.;
      }
    }

    d0 = tr->getD0();
    d0err = sqrt(tr->getCovMatrix()[tpar::d0d0]);
    z0 = tr->getZ0();
    z0err = sqrt(tr->getCovMatrix()[tpar::z0z0]);

    e = tr->E();
    pt = tr->Pt();

    chi2 = tr->getChi2();
    ndf = tr->getNdf();
    vtxftdhits = tr->getVtxHits() + tr->getFtdHits();

    _nt->Fill(trackseed, invtx, d0, d0err, z0, z0err, e, pt, chi2, ndf, vtxftdhits);
  }
}

void VertexAnalysis::end() {
  _file->Write();
  _file->Close();
}

#if 1
//const int MyAnalysis::Blist[] = { 511, 521, 531, 541, 5122, 5132, 5232, 5332};
//const int MyAnalysis::Clist[] = { 411, 421, 431, 4122, 4132, 4232, 4332};
//const int MyAnalysis::Olist[] = { 11, 13, 15, 22, 130, 211, 310, 321, 2112, 2212, 3112, 3122, 3212, 3222, 3312, 3322, 3334 };
//const int MyAnalysis::V0list[] = { 22, 310, 3122 };
std::vector<int> MyAnalysis::Blist = { 511, 521, 531, 541, 5122, 5132, 5232, 5332};
std::vector<int> MyAnalysis::Clist = { 411, 421, 431, 4122, 4132, 4232, 4332};
std::vector<int> MyAnalysis::Olist = { 11, 13, 15, 22, 130, 211, 310, 321, 2112, 2212, 3112, 3122, 3212, 3222, 3312, 3322, 3334 };
std::vector<int> MyAnalysis::V0list = { 22, 310, 3122 };
const float MyAnalysis::aSmallNumber = 0.001;

void MyAnalysis::init(Parameters* param) {
  Algorithm::init(param);
  _privtxname = param->get("PrimaryVertexCollectionName",string("PrimaryVertex"));
  _secvtxname = param->get("SecondaryVertexCollectionName",string("BuildUpVertex"));
  _v0vtxname = param->get("V0VertexCollectionName",string("BuildUpVertex_V0"));
  _nEvt = 1;
  string filename = param->get("OutputRootFile",string("MyAnalysis_output.root"));
  _file = new TFile(filename.c_str(),"RECREATE");
          _trktree = new TTree("VertexTrack","Vertex track tree");
	  _trktree->Branch("nevt"              , &_evdata.nevt               , "nevt/I"              );
	  _trktree->Branch("ipTruthx"          , &_evdata.ipTruthx           , "ipTruthx"            );
	  _trktree->Branch("ipTruthy"          , &_evdata.ipTruthy           , "ipTruthy"            );
	  _trktree->Branch("ipTruthz"          , &_evdata.ipTruthz           , "ipTruthz"            );
	  _trktree->Branch("mcx_vtx"           , &_trkdata.mcx_vtx           , "mcx_vtx"             );
	  _trktree->Branch("mcy_vtx"           , &_trkdata.mcy_vtx           , "mcy_vtx"             );
	  _trktree->Branch("mcz_vtx"           , &_trkdata.mcz_vtx           , "mcz_vtx"             );
	  _trktree->Branch("mcx_evtx"          , &_trkdata.mcx_evtx          , "mcx_evtx"            );
	  _trktree->Branch("mcy_evtx"          , &_trkdata.mcy_evtx          , "mcy_evtx"            );
	  _trktree->Branch("mcz_evtx"          , &_trkdata.mcz_evtx          , "mcz_evtx"            );
	  _trktree->Branch("mce"               , &_trkdata.mce               , "mce/D"               );
	  _trktree->Branch("mcp"               , &_trkdata.mcp               , "mcp/D"               );
	  _trktree->Branch("mcpt"              , &_trkdata.mcpt              , "mcpt/D"              );
	  _trktree->Branch("mcpx"              , &_trkdata.mcpx              , "mcpx/D"              );
	  _trktree->Branch("mcpy"              , &_trkdata.mcpy              , "mcpy/D"              );
	  _trktree->Branch("mcpz"              , &_trkdata.mcpz              , "mcpz/D"              );
	  _trktree->Branch("mcpdg"             , &_trkdata.mcpdg             , "mcpdg/I"             );
	  _trktree->Branch("mccostheta"        , &_trkdata.mccostheta        , "mccostheta/D"        );
	  _trktree->Branch("mcvpdg"            , &_trkdata.mcvpdg            , "mcvpdg/I"            );
	  _trktree->Branch("mcve"              , &_trkdata.mcve              , "mcve/D"              );
	  _trktree->Branch("mcvp"              , &_trkdata.mcvp              , "mcvp/D"              );
	  _trktree->Branch("mcvpt"             , &_trkdata.mcvpt             , "mcvpt/D"             );
	  _trktree->Branch("mcvpx"             , &_trkdata.mcvpx             , "mcvpx/D"             );
	  _trktree->Branch("mcvpy"             , &_trkdata.mcvpy             , "mcvpy/D"             );
	  _trktree->Branch("mcvpz"             , &_trkdata.mcvpz             , "mcvpz/D"             );
	  _trktree->Branch("mcvcostheta"       , &_trkdata.mcvcostheta       , "mcvcostheta/D"       );
	  _trktree->Branch("isFromP"           , &_trkdata.isFromP           , "isFromP/I"           );
	  _trktree->Branch("isFromB"           , &_trkdata.isFromB           , "isFromB/I"           );
	  _trktree->Branch("isFromC"           , &_trkdata.isFromC           , "isFromC/I"           );
	  _trktree->Branch("isFromO"           , &_trkdata.isFromO           , "isFromO/I"           );
	  _trktree->Branch("isFromBmbkg"       , &_trkdata.isFromBmbkg       , "isFromBmbkg/O"       );
	  _trktree->Branch("rce"               , &_trkdata.rce               , "rce/D"               );
	  _trktree->Branch("rcp"               , &_trkdata.rcp               , "rcp/D"               );
	  _trktree->Branch("rcpt"              , &_trkdata.rcpt              , "rcpt/D"              );
	  _trktree->Branch("rcpx"              , &_trkdata.rcpx              , "rcpx/D"              );
	  _trktree->Branch("rcpy"              , &_trkdata.rcpy              , "rcpy/D"              );
	  _trktree->Branch("rcpz"              , &_trkdata.rcpz              , "rcpz/D"              );
	  _trktree->Branch("rccostheta"        , &_trkdata.rccostheta        , "rccostheta/D"        );
	  _trktree->Branch("d0"                , &_trkdata.d0                , "d0/D"                );
	  _trktree->Branch("d0err"             , &_trkdata.d0err             , "d0err/D"             );
	  _trktree->Branch("phi"               , &_trkdata.phi               , "phi/D"               );
	  _trktree->Branch("phierr"            , &_trkdata.phierr            , "phierr/D"            );
	  _trktree->Branch("omg"               , &_trkdata.omg               , "omg/D"               );
	  _trktree->Branch("omgerr"            , &_trkdata.omgerr            , "omgerr/D"            );
	  _trktree->Branch("z0"                , &_trkdata.z0                , "z0/D"                );
	  _trktree->Branch("z0err"             , &_trkdata.z0err             , "z0err/D"             );
	  _trktree->Branch("tanl"              , &_trkdata.tanl              , "tanl/D"              );
	  _trktree->Branch("tanlerr"           , &_trkdata.tanlerr           , "tanlerr/D"           );
	  _trktree->Branch("isAssociatedToPvtx", &_trkdata.isAssociatedToPvtx, "isAssociatedToPvtx/O");
	  _trktree->Branch("isAssociatedToSvtx", &_trkdata.isAssociatedToSvtx, "isAssociatedToSvtx/O");
	  _trktree->Branch("isCorrectVertex"   , &_trkdata.isCorrectVertex   , "isCorrectVertex/O"   );
	  _trktree->Branch("isCorrectDChain"   , &_trkdata.isCorrectDChain   , "isCorrectDChain/O"   );
          _trktree->Branch("distancePvtxToSvtx", &_trkdata.distancePvtxToSvtx, "distancePvtxToSvtx"  );

          _pvtxtree = new TTree("PrimaryVertex","Primary Vertex tree");
	  _pvtxtree->Branch("nevt"       , &_evdata.nevt         , "nevt/I"      );
	  _pvtxtree->Branch("rcx"        , &_vtxdata.rcx         , "rcx"         );
	  _pvtxtree->Branch("rcy"        , &_vtxdata.rcy         , "rcy"         );
	  _pvtxtree->Branch("rcz"        , &_vtxdata.rcz         , "rcz"         );
	  _pvtxtree->Branch("rcpt"       , &_vtxdata.rcpt        , "rcpt"        );
	  _pvtxtree->Branch("rcpx"       , &_vtxdata.rcpx        , "rcpx"        );
	  _pvtxtree->Branch("rcpy"       , &_vtxdata.rcpy        , "rcpy"        );
	  _pvtxtree->Branch("rcpz"       , &_vtxdata.rcpz        , "rcpz"        );
	  _pvtxtree->Branch("rce"        , &_vtxdata.rce         , "rce"         );
	  _pvtxtree->Branch("rcp"        , &_vtxdata.rcp         , "rcp"         );
	  _pvtxtree->Branch("rccostheta" , &_vtxdata.rccostheta  , "rccostheata" );
	  _pvtxtree->Branch("rcmass"     , &_vtxdata.rcmass      , "rcmass"      );
	  _pvtxtree->Branch("xerr"       , &_vtxdata.xerr        , "xerr"        );
	  _pvtxtree->Branch("yerr"       , &_vtxdata.yerr        , "yerr"        );
	  _pvtxtree->Branch("zerr"       , &_vtxdata.zerr        , "zerr"        );
	  _pvtxtree->Branch("chi2"       , &_vtxdata.chi2        , "chi2/D"      );
	  _pvtxtree->Branch("rcntrk"     , &_vtxdata.rcntrk      , "rcntrk/I"    );
	  _pvtxtree->Branch("ipTruthx"   , &_vtxdata.ipTruthx    , "ipTruthx"    );
	  _pvtxtree->Branch("ipTruthy"   , &_vtxdata.ipTruthy    , "ipTruthy"    );
	  _pvtxtree->Branch("ipTruthz"   , &_vtxdata.ipTruthz    , "ipTruthz"    );
	  _pvtxtree->Branch("mcx"        , &_vtxdata.mcx         , "mcx"         );
	  _pvtxtree->Branch("mcy"        , &_vtxdata.mcy         , "mcy"         );
	  _pvtxtree->Branch("mcz"        , &_vtxdata.mcz         , "mcz"         );
	  _pvtxtree->Branch("type"       , &_vtxdata.type        , "type/I"      );
	  _pvtxtree->Branch("wegiht"     , &_vtxdata.weight      , "weight"      );
	  _pvtxtree->Branch("mcntrk"     , &_vtxdata.mcntrk      , "mcntrk/I"    );
	  _pvtxtree->Branch("mcpt"       , &_vtxdata.mcpt        , "mcpt"        );
	  _pvtxtree->Branch("mcpx"       , &_vtxdata.mcpx        , "mcpx"        );
	  _pvtxtree->Branch("mcpy"       , &_vtxdata.mcpy        , "mcpy"        );
	  _pvtxtree->Branch("mcpz"       , &_vtxdata.mcpz        , "mcpz"        );
	  _pvtxtree->Branch("mcp"        , &_vtxdata.mcp         , "mcp"         );
	  _pvtxtree->Branch("mce"        , &_vtxdata.mce         , "mce"         );
	  _pvtxtree->Branch("mccostheta" , &_vtxdata.mccostheta  , "mccostheata" );
	  _pvtxtree->Branch("mcmass"     , &_vtxdata.mcmass      , "mcmass"      );
	  _pvtxtree->Branch("mcpdg"      , &_vtxdata.mcpdg       , "mcpdg/I"     );

          _svtxtree = new TTree("SecondaryVertex","Secondary Vertex tree");
	  _svtxtree->Branch("nevt"       , &_evdata.nevt         , "nevt/I"      );
	  _svtxtree->Branch("rcx"        , &_vtxdata.rcx         , "rcx"         );
	  _svtxtree->Branch("rcy"        , &_vtxdata.rcy         , "rcy"         );
	  _svtxtree->Branch("rcz"        , &_vtxdata.rcz         , "rcz"         );
	  _svtxtree->Branch("rcpt"       , &_vtxdata.rcpt        , "rcpt"        );
	  _svtxtree->Branch("rcpx"       , &_vtxdata.rcpx        , "rcpx"        );
	  _svtxtree->Branch("rcpy"       , &_vtxdata.rcpy        , "rcpy"        );
	  _svtxtree->Branch("rcpz"       , &_vtxdata.rcpz        , "rcpz"        );
	  _svtxtree->Branch("rce"        , &_vtxdata.rce         , "rce"         );
	  _svtxtree->Branch("rcp"        , &_vtxdata.rcp         , "rcp"         );
	  _svtxtree->Branch("rccostheta" , &_vtxdata.rccostheta  , "rccostheata" );
	  _svtxtree->Branch("rcmass"     , &_vtxdata.rcmass      , "rcmass"      );
	  _svtxtree->Branch("xerr"       , &_vtxdata.xerr        , "xerr"        );
	  _svtxtree->Branch("yerr"       , &_vtxdata.yerr        , "yerr"        );
	  _svtxtree->Branch("zerr"       , &_vtxdata.zerr        , "zerr"        );
	  _svtxtree->Branch("chi2"       , &_vtxdata.chi2        , "chi2/D"      );
	  _svtxtree->Branch("rcntrk"     , &_vtxdata.rcntrk      , "rcntrk/I"    );
	  _svtxtree->Branch("ipTruthx"   , &_vtxdata.ipTruthx    , "ipTruthx"    );
	  _svtxtree->Branch("ipTruthy"   , &_vtxdata.ipTruthy    , "ipTruthy"    );
	  _svtxtree->Branch("ipTruthz"   , &_vtxdata.ipTruthz    , "ipTruthz"    );
	  _svtxtree->Branch("mcx"        , &_vtxdata.mcx         , "mcx"         );
	  _svtxtree->Branch("mcy"        , &_vtxdata.mcy         , "mcy"         );
	  _svtxtree->Branch("mcz"        , &_vtxdata.mcz         , "mcz"         );
	  _svtxtree->Branch("type"       , &_vtxdata.type        , "type/I"      );
	  _svtxtree->Branch("wegiht"     , &_vtxdata.weight      , "weight"      );
	  _svtxtree->Branch("mcntrk"     , &_vtxdata.mcntrk      , "mcntrk/I"    );
	  _svtxtree->Branch("mcpt"       , &_vtxdata.mcpt        , "mcpt"        );
	  _svtxtree->Branch("mcpx"       , &_vtxdata.mcpx        , "mcpx"        );
	  _svtxtree->Branch("mcpy"       , &_vtxdata.mcpy        , "mcpy"        );
	  _svtxtree->Branch("mcpz"       , &_vtxdata.mcpz        , "mcpz"        );
	  _svtxtree->Branch("mcp"        , &_vtxdata.mcp         , "mcp"         );
	  _svtxtree->Branch("mce"        , &_vtxdata.mce         , "mce"         );
	  _svtxtree->Branch("mccostheta" , &_vtxdata.mccostheta  , "mccostheata" );
	  _svtxtree->Branch("mcmass"     , &_vtxdata.mcmass      , "mcmass"      );
	  _svtxtree->Branch("mcpdg"      , &_vtxdata.mcpdg       , "mcpdg/I"     );
          _svtxtree->Branch("dl"         , &_vtxdata.dl          , "dl/D"        );
          _svtxtree->Branch("dt1"        , &_vtxdata.dt1         , "dt1/D"       );
          _svtxtree->Branch("dt2"        , &_vtxdata.dt2         , "dt2/D"       );
          _svtxtree->Branch("distancePvtxToSvtx", &_vtxdata.distancePvtxToSvtx, "distancePvtxToSvtx");
  std::cerr << "init called." << std::endl;
}

void MyAnalysis::process() {
  std::cerr << "### ev = " << _nEvt << std::endl;
  _evdata.nevt = _nEvt;
  //std::cerr << "process called." << std::endl;
  const Vertex* pvtx = Event::Instance()->getPrimaryVertex(_privtxname.c_str());
  std::vector<const Vertex*> svtxs = Event::Instance()->getSecondaryVertices(_secvtxname.c_str());
  std::vector<const Vertex*> v0vtx = Event::Instance()->getSecondaryVertices(_v0vtxname.c_str());
  const TrackVec& tracks = Event::Instance()->getTracks();
  const NeutralVec& neutrals = Event::Instance()->getNeutrals();
  const MCParticleVec& mcps = Event::Instance()->getMCParticles();
  _mcps = &mcps;
  //std::cerr << pvtx << std::endl;
  //std::cerr << svtx.size() << std::endl;
  std::cerr << "v0vtx.size() = " << v0vtx.size() << std::endl;
  //std::cerr << tracks.size() << std::endl;
  //std::cerr << neutrals.size() << std::endl;
  //std::cerr << mcps.size() << std::endl;

  std::vector<MCVertex> BCs; // B or C hadrons to be found in this event.
  std::vector<MCVertex> V0s; // V0 to be found in this event.

  //for (unsigned int itrk = 0; itrk < tracks.size(); itrk++) {
  //  const MCParticle* mcp = tracks[itrk]->getMcp();
  //}

  _mcpIndex.clear(); // clean up 
  for (unsigned int imcp = 0; imcp < mcps.size(); imcp++) {
    const MCParticle* mcp = mcps[imcp];
    _mcpIndex.insert(std::map<const MCParticle*,int>::value_type(mcp,imcp)); 
  }

  // Define a virtual MC particle for IP.
  //TVector3 ipTruth = getIPTruth(mcps[0]);
  TVector3 ipTruth = getIPTruth();
  //MCParticle* ipmc = new MCParticle;
//std::cerr << "### createVtxMCParticle for PrimaryVertex" << std::endl;
  MCParticle* ipmc = createVtxMCParticle(ipTruth);
  //ipmc->setVertex(ipTruth);
  //setTracksOf(ipmc,mcps);
  //setTracksOf(ipmc);
  ipmc->setPDG(100000000);
  _evdata.ipTruthx = ipTruth.X();
  _evdata.ipTruthy = ipTruth.Y();
  _evdata.ipTruthz = ipTruth.Z();

  std::set<MCParticle*> BmbkgIPMCParticles;

  std::map<const MCParticle*,MCParticleExt> mcpMap;
  std::set<const MCParticle*> BCRegisteredList; // This will be used to check if the BC has already been registered.
  std::set<const MCParticle*> V0RegisteredList; // This will be used to check if the V0 has already been registered.
  for (unsigned int imcp = 0; imcp < mcps.size(); imcp++) {

    const MCParticle* mcp = mcps[imcp];

    MCParticleExt mcpe;
    mcpe.setMCParticle(mcp);

    int mcpdg = mcp->getPDG();
    MCParticleExt Bpart, Cpart, Opart, V0part;
//if (imcp==125) std::cerr << " ############# Check this part ! pdg " << mcpdg << " ndaughters = " << mcp->getDaughters().size() << std::endl;
    if (mcp->getParent()&&mcp->getCharge()!=0&&
        abs(mcpdg)>6&&  // remove bare quarks
        mcpdg!=92&&mcpdg!=94 // remove intermediate states
       //!(abs(mcpdg)==11&&mcp->getGeneratorStatus()==2) // remove beam particles
       // && mcp->getDaughters().size()>0
       ) { // corresponds to obserbed tracks
    
      //std::cerr << itrk << " " <<  mcp->getPDG() << std::endl;
      Bpart.setMCParticle(mcp);    
      Cpart.setMCParticle(mcp);    
      Opart.setMCParticle(mcp);    
      V0part.setMCParticle(mcp);    

      updateSemistableAncestorOf(&Bpart, Blist);
      if (Bpart.getAncestor()) { // if found
         std::set<const MCParticle*>::iterator mcpitr = BCRegisteredList.find(Bpart.getAncestor());
         if (mcpitr == BCRegisteredList.end()) {
           BCRegisteredList.insert(Bpart.getAncestor()); 
           MCVertex mcvtx;
           mcvtx.mcpe.setMCParticle(Bpart.getAncestor()); 
           mcvtx.mcpe.setAncestor(Bpart.getAncestor());
           mcvtx.mcpe.setType(2);
           BCs.push_back(mcvtx);
         }
      }
//if (imcp==125) std::cerr << " ############# Check this part !" << std::endl;
      updateSemistableAncestorOf(&Cpart, Clist);
      if (Cpart.getAncestor()) {
         std::set<const MCParticle*>::iterator mcpitr = BCRegisteredList.find(Cpart.getAncestor());
         if (mcpitr == BCRegisteredList.end()) {
           BCRegisteredList.insert(Cpart.getAncestor()); 
           MCVertex mcvtx;
           mcvtx.mcpe.setMCParticle(Cpart.getAncestor());
           mcvtx.mcpe.setAncestor(Cpart.getAncestor());
           mcvtx.mcpe.setType(3);
           BCs.push_back(mcvtx);
         }
      }

      updateSemistableAncestorOf(&Opart, Olist); 

      updateSemistableAncestorOf(&V0part, V0list); 
      // check if the daughtor has two charged particles that were actually found.
      if (V0part.getAncestor()) {
         bool isV0 = false;
         if (V0part.getAncestor()->getDaughters().size()==2) {
           const MCParticle* dau0 = V0part.getAncestor()->getDaughters()[0];
           const MCParticle* dau1 = V0part.getAncestor()->getDaughters()[1];
           if (dau0->getCharge() !=0 && dau1->getCharge() !=0 ) {
               bool isTrackFound0 = false;
               bool isTrackFound1 = false;
               for (unsigned int itrk = 0; itrk < tracks.size(); itrk++) {
                 const Track* trk = tracks[itrk];
                 const MCParticle* mcp = trk->getMcp();
                 if (dau0 == mcp) isTrackFound0 = true; 
                 if (dau1 == mcp) isTrackFound1 = true; 
               }            
               if (isTrackFound0&&isTrackFound1) isV0 = true; 
           }
 

         }
         if (!isV0) V0part.setAncestor(0); // Reset origin if it doesn't have two charged particles.
      }
    }

    if (V0part.getAncestor()) {
      std::map<const MCParticle*,int>::iterator mcindexitr;
      mcindexitr = _mcpIndex.find(V0part.getAncestor());
      int index = -1;
      if (mcindexitr != _mcpIndex.end()) index = mcindexitr->second;
      std::set<const MCParticle*>::iterator mcpitr = V0RegisteredList.find(V0part.getAncestor());
      if (mcpitr == V0RegisteredList.end()) {
        V0RegisteredList.insert(V0part.getAncestor()); 
        MCVertex mcvtx;
        mcvtx.mcpe.setMCParticle(V0part.getAncestor());
        mcvtx.mcpe.setAncestor(V0part.getAncestor());
        mcvtx.mcpe.setType(9); // V0 
        V0s.push_back(mcvtx);
//TVector3 pos(V0part.getAncestor()->getVertex());
//std::cerr << "V0 particle found : index = " << index << " " << V0part.getAncestor()->getPDG() << " (" << pos.X() << ", " << pos.Y() << ", " << pos.Z() << ") "<< V0part.getAncestor()->getDaughters()[0]->getPDG() << " " << V0part.getAncestor()->getDaughters()[1]->getPDG() << std::endl;
//std::cerr << "Myself = " << mcp->getPDG() << std::endl;
      }
    }
    if (isFromBeambackground(mcp)) {
      MCParticle* ipbeambkg = 0;
      const MCParticle* mcporig = getOriginMCParticle(mcp);
      TVector3 ipBeamBkg(mcporig->getVertex());
      //std::map<const MCParticle*, MCParticle*>::iterator bmbkgipitr;
      std::set<MCParticle*>::iterator bmbkgipitr;
      //bmbkgipitr = BmbkgIPMCParticles.find(mcporig);
      for (bmbkgipitr=BmbkgIPMCParticles.begin();bmbkgipitr!=BmbkgIPMCParticles.end();bmbkgipitr++) {
        TVector3 refpos((*bmbkgipitr)->getVertex()); 
        if ((refpos-ipBeamBkg).Mag()<aSmallNumber) { 
          ipbeambkg = *bmbkgipitr;
        }
      }
      if (!ipbeambkg) {
//std::cerr << "### createVtxMCParticle for Beam background. (" << ipBeamBkg.X() << ", " << ipBeamBkg.Y() << ", " << ipBeamBkg.Z() << ") mcporig = " << mcporig << std::endl;
         ipbeambkg = createVtxMCParticle(ipBeamBkg);
         ipbeambkg->setPDG(100000001);
         BmbkgIPMCParticles.insert(ipbeambkg);
      }
      mcpe.setType(0);
      mcpe.setMCParticle(ipbeambkg);
      mcpe.setAncestor(ipbeambkg);
    } else if (isFromPrimaryVertex(mcp)) {
      mcpe.setType(1);
      mcpe.setMCParticle(ipmc);
      mcpe.setAncestor(ipmc);
    } else {

      // Find the nearest ancester, e.g. B->C->O-> case, we define this particle as O.
      std::vector<MCParticleExt> ancestry;
      if (Bpart.getAncestor()) ancestry.push_back(Bpart); 
      if (Cpart.getAncestor()) ancestry.push_back(Cpart); 
      if (Opart.getAncestor()) ancestry.push_back(Opart); 

      if (ancestry.size()) {
         sort(ancestry.begin(),ancestry.end(),MCParticleExt::DescendingGeneration);
         if      (ancestry[0].getAncestor()==Bpart.getAncestor()) {
           mcpe.setType(2);
           mcpe.setAncestor(Bpart.getAncestor());
         }
         else if (ancestry[0].getAncestor()==Cpart.getAncestor()) {
           mcpe.setType(3);
           mcpe.setAncestor(Cpart.getAncestor());
         }
         else if (ancestry[0].getAncestor()==Opart.getAncestor()) {
            mcpe.setType(4);
            mcpe.setAncestor(Opart.getAncestor());
         } 
      } else { // non of the above, e.g. 0 daughter particles.
//std::cerr << "imcp = " << imcp << std::endl;
//std::abort();
         mcpe.setType(5);
         mcpe.setAncestor(getParentIfSamePDG(mcp)); // if the particle radiates any particle, it shouldn't be recognized as its decay.
      } 
    }
    mcpMap.insert(map<const MCParticle*,MCParticleExt>::value_type(mcp,mcpe));
  }
  std::cerr << "map size : " << mcpMap.size() << std::endl;
  // MCParticle loop end.
  

  std::cerr << "Primary Vertex" << std::endl;
  std::vector<const Vertex*> pvtxs;
  pvtxs.push_back(pvtx);
#if 0 // for single particle sample
  checkVertices(pvtxs,BCs,mcpMap,_pvtxtree);
  std::cerr << "Secondary Vertex" << std::endl;
  checkVertices(svtxs,BCs,mcpMap,_svtxtree);
#endif

  // Track loop.
  for (unsigned int itrk = 0; itrk < tracks.size(); itrk++) {
    const Track* trk = tracks[itrk];
    const MCParticle* mcp = trk->getMcp();
    MCParticleExt trkmcpe;
    trkmcpe.setMCParticle(mcp); // default MCParticle. This will be replaced if necessary. 

    _trkdata.mcpdg  = mcp->getPDG();
    _trkdata.mcx_vtx  = getParentIfSamePDG(mcp)->getVertex().X();
    _trkdata.mcy_vtx  = getParentIfSamePDG(mcp)->getVertex().Y();
    _trkdata.mcz_vtx  = getParentIfSamePDG(mcp)->getVertex().Z();
    _trkdata.mcx_evtx = mcp->getEndVertex().X();
    _trkdata.mcy_evtx = mcp->getEndVertex().Y();
    _trkdata.mcz_evtx = mcp->getEndVertex().Z();
    _trkdata.mcpt   = mcp->Pt();
    _trkdata.mcpx   = mcp->Px();
    _trkdata.mcpy   = mcp->Py();
    _trkdata.mcpz   = mcp->Pz();
    _trkdata.mce    = mcp->E();
    _trkdata.mcp    = mcp->Vect().Mag();
    _trkdata.mccostheta = _trkdata.mcpz/_trkdata.mcp;
    _trkdata.rcpt   = trk->Pt();
    _trkdata.rcpx   = trk->Px();
    _trkdata.rcpy   = trk->Py();
    _trkdata.rcpz   = trk->Pz();
    _trkdata.rce    = trk->E();
    _trkdata.rcp    = trk->Vect().Mag();
    _trkdata.rccostheta = _trkdata.rcpz/_trkdata.rcp;
    _trkdata.d0     = trk->getD0();
    _trkdata.z0     = trk->getZ0();
    _trkdata.phi    = trk->getPhi();
    _trkdata.omg    = trk->getOmega();
    _trkdata.tanl   = trk->getTanLambda();
    _trkdata.d0err     = TMath::Sqrt(trk->getCovMatrix()[tpar::d0d0]);
    _trkdata.z0err     = TMath::Sqrt(trk->getCovMatrix()[tpar::z0z0]);
    _trkdata.phierr    = TMath::Sqrt(trk->getCovMatrix()[tpar::phph]);
    _trkdata.omgerr    = TMath::Sqrt(trk->getCovMatrix()[tpar::omom]);
    _trkdata.tanlerr   = TMath::Sqrt(trk->getCovMatrix()[tpar::tdtd]);
    _trkdata.isFromP     = isFromPrimaryVertex(mcp);
    _trkdata.isFromBmbkg = isFromBeambackground(mcp);

    // Do same as in MCParticle loop.
    MCParticleExt Bpart, Cpart, Opart;
    Bpart.setMCParticle(mcp);    
    Cpart.setMCParticle(mcp);    
    Opart.setMCParticle(mcp);    
    updateSemistableAncestorOf(&Bpart, Blist); 
    if (Bpart.getAncestor()) _trkdata.isFromB = Bpart.getAncestor()->getPDG(); 
    else                   _trkdata.isFromB = 0;
    updateSemistableAncestorOf(&Cpart, Clist); 
    if (Cpart.getAncestor()) _trkdata.isFromC = Cpart.getAncestor()->getPDG(); 
    else                   _trkdata.isFromC = 0;
    updateSemistableAncestorOf(&Opart, Olist); 
    if (Opart.getAncestor()) _trkdata.isFromO = Opart.getAncestor()->getPDG(); 
    else                   _trkdata.isFromO = 0;
    // Find the nearest ancester, e.g. B->C case, we define this particle as C.
    std::vector<MCParticleExt> ancestry;
    if (Bpart.getAncestor()) ancestry.push_back(Bpart); 
    if (Cpart.getAncestor()) ancestry.push_back(Cpart); 
    if (Opart.getAncestor()) ancestry.push_back(Opart); 
    if (ancestry.size()) {
       sort(ancestry.begin(),ancestry.end(),MCParticleExt::DescendingGeneration);
       if      (ancestry[0].getAncestor()==Bpart.getAncestor()) {
         trkmcpe.setMCParticle(Bpart.getAncestor());
         trkmcpe.setAncestor(Bpart.getAncestor());
       }
       else if (ancestry[0].getAncestor()==Cpart.getAncestor()) {
         trkmcpe.setMCParticle(Cpart.getAncestor());
         trkmcpe.setAncestor(Cpart.getAncestor());
       }
       else if (ancestry[0].getAncestor()==Opart.getAncestor()) {
         trkmcpe.setMCParticle(Opart.getAncestor());
         trkmcpe.setAncestor(Opart.getAncestor());
       }
    }

    _trkdata.isAssociatedToPvtx = false;
    //check if this track is associated to secondary vertieces.
    for (int ipv = 0; ipv < pvtxs.size(); ipv++) {
      const Vertex* pvtx = pvtxs[ipv];
      const std::vector<const Track*> pvtracks = pvtx->getTracks();
      std::vector<const Track*>::const_iterator pvtrackitr; 
      pvtrackitr = std::find(pvtracks.begin(),pvtracks.end(),trk);
      if (pvtrackitr != pvtracks.end()) { // found the pvtx that is associated with this track.
        _trkdata.isAssociatedToPvtx = true;
      }
    }
//std::cerr << "  itrk = " << itrk;
//if (trkmcpe.getAncestor()) std::cerr << " orig : " << trkmcpe.getAncestor()->getPDG() << " ("<< trkmcpe.getAncestor() << ") "<< std::endl;
//else std::cerr << " no orig info." << std::endl;
    _trkdata.isAssociatedToSvtx = false;
    _trkdata.isCorrectVertex = false; 
    _trkdata.isCorrectDChain = false; // valid only for B,C vertices.
    bool isSearchEnd = false;
    //check if this track is associated to secondary vertieces.
    for (int isv = 0; isv < svtxs.size(); isv++) {
      const Vertex* svtx = svtxs[isv];
      const std::vector<const Track*> svtracks = svtx->getTracks();
      std::vector<const Track*>::const_iterator svtrackitr; 
      svtrackitr = std::find(svtracks.begin(),svtracks.end(),trk);
      if (svtrackitr != svtracks.end()) { // found the svtx that is associated with this track.
        _trkdata.isAssociatedToSvtx = true;
        for (int ibc = 0; ibc < BCs.size(); ibc++) {
          MCVertex& mcv = BCs[ibc];
//std::cerr << "ibc = " << ibc << " mcv.mcpe.getMCParticle()->getPDG() = " << mcv.mcpe.getMCParticle()->getPDG() << std::endl;
          if (mcv.matchedvtxrec==svtx) {  // found the mcvtx that is matched to this svtx.
//std::cerr << "### mcvtx matched" << std::endl;
             // check if this track origin has a C-hadron descendant.
             if(updateSemistableDescendantOf(&trkmcpe, Clist)){       // confirm trkmcpe.getDescendant() is found.
               if (mcv.mcpe.getAncestor()==trkmcpe.getDescendant()) { // track originated from B (trk) but associated to C (mcv).
                 _trkdata.isCorrectDChain = true;
//std::cerr << "track orig : " << trkmcpe.getAncestor()->getPDG() << " mcvtx : " << mcv.mcpe.getAncestor()->getPDG() << " track orig has " << trkmcpe.getDescendant()->getPDG() << " in descendants." << std::endl;
               }
             }
             // check if MC origin has a C-hadron descendant.
//std::cerr << "trkmcpe.getAncestor()->getPDG() = " << trkmcpe.getAncestor()->getPDG() << " mcv.mcpe.getAncestor()->getPDG() = " << mcv.mcpe.getAncestor()->getPDG() << " mcv.mcpe.getMCParticle()->getPDG() = " << mcv.mcpe.getMCParticle()->getPDG() << " type: " << mcv.mcpe.getType() << std::endl;
             mcv.mcpe.setGeneration(0);
             if (updateSemistableDescendantOf(&mcv.mcpe, Clist)){     // confirm mvc.mcpe.getDescendant() is found.
               if (mcv.mcpe.getDescendant()==trkmcpe.getAncestor()) { // track originated from C (trk) but associated to B (mcv).
                 _trkdata.isCorrectDChain = true;
//std::cerr << "track orig : " << trkmcpe.getAncestor()->getPDG() << " mcvtx : " << BCs[ibc].mcpe.getAncestor()->getPDG() << " mcvtx has " << BCs[ibc].mcpe.getDescendant()->getPDG() << " in descendants." << std::endl;
               }
             } 
             if (mcv.mcpe.getAncestor()==trkmcpe.getAncestor() && trkmcpe.getAncestor()) {
//std::cerr << "Correctly found : pdg " <<  trkmcpe.getAncestor()->getPDG() << std::endl;
                _trkdata.isCorrectDChain = true;
                _trkdata.isCorrectVertex = true;
             }
             isSearchEnd = true;

             // Before end, fill info related to vertex.
             TVector3 pvpos(pvtx->getPos());
             TVector3 svpos(svtx->getPos());
             _trkdata.distancePvtxToSvtx = (svpos-pvpos).Mag();
             _trkdata.mcvpdg = mcv.mcpe.getMCParticle()->getPDG();
             _trkdata.mcve   = mcv.mcpe.getMCParticle()->E();
             _trkdata.mcvp   = mcv.mcpe.getMCParticle()->Vect().Mag();
             _trkdata.mcvpt  = mcv.mcpe.getMCParticle()->Pt();
             _trkdata.mcvpx  = mcv.mcpe.getMCParticle()->Px();
             _trkdata.mcvpy  = mcv.mcpe.getMCParticle()->Py();
             _trkdata.mcvpz  = mcv.mcpe.getMCParticle()->Pz();
             _trkdata.mcvcostheta = _trkdata.mcvpz/_trkdata.mcvp;
//std::cerr << "_trkdata.mcvpdg = " << _trkdata.mcvpdg << " " 
//          << "_trkdata.mcve = "   << _trkdata.mcve   << " " 
//          << "_trkdata.mcvp = "   << _trkdata.mcvp   << " " 
//          << "_trkdata.mcvpt = "  << _trkdata.mcvpt  << " " 
//          << "_trkdata.mcvpx = "  << _trkdata.mcvpx  << " " 
//          << "_trkdata.mcvpy = "  << _trkdata.mcvpy  << " " 
//          << "_trkdata.mcvpz = "  << _trkdata.mcvpz  << " " 
//          << "_trkdata.mcvcostheta = "  << _trkdata.mcvcostheta  << std::endl; 

             break;
          }
        } 
      }
      if (isSearchEnd) break;
    }
    _trktree->Fill();
  }
#if 1
  const bool dump = true;
  if (dump) {
    std::cerr << "#### Summary (Secondary vertices) ####" << std::endl;
    std::cerr << "# of vertices to be found = " << BCs.size() << std::endl;
    int nbcvtxrec = 0;
    int nbcvtxmc = 0;
    int nbcvtxrec_b = 0;
    int nbcvtxmc_b = 0;
    int nbcvtxrec_c = 0;
    int nbcvtxmc_c = 0;
    std::vector<MCVertex>::iterator mcvtxitr;
    for (mcvtxitr = BCs.begin(); mcvtxitr!=BCs.end(); mcvtxitr++) {
      nbcvtxmc++;
      if (mcvtxitr->mcpe.getType()==2) nbcvtxmc_b++;
      if (mcvtxitr->mcpe.getType()==3) nbcvtxmc_c++;
      if (dump) std::cerr << "  " << setw(4) << nbcvtxmc << ") " << setw(5) << mcvtxitr->mcpe.getAncestor()->getPDG() 
                          << " (" << (*mcvtxitr).mcpe.getAncestor() << ") ";
        if (mcvtxitr->matchedvtxrec) {
          nbcvtxrec++;
          if ((*mcvtxitr).mcpe.getType()==2) nbcvtxrec_b++;
          if ((*mcvtxitr).mcpe.getType()==3) nbcvtxrec_c++;
          std::cerr << " Found.";
          //TVector3 vpos(mcvtxitr->mcpe.getAncestor()->getEndVertex());
          //std::cerr << " TEST : (" << vpos.X() << ", " << vpos.Y() << ", " << vpos.Z() << ") "; 
        } else {
          std::cerr << " Not found.";
          TVector3 vpos(mcvtxitr->mcpe.getAncestor()->getEndVertex());
          std::cerr << " Vertex : (" << vpos.X() << ", " << vpos.Y() << ", " << vpos.Z() << ") "; 
          //int nfoundtracks = 0;
	  //for(int ipfo=0; ipfo < colPFO->getNumberOfElements(); ipfo++) {
          //   ReconstructedParticle *part = dynamic_cast<ReconstructedParticle*>( colPFO->getElementAt( ipfo ) );
          //   if (part->getTracks().size()) {
          //     MCParticle *mcp = getBestMCParticleOf(part,relnav);
          //     if ((*mcvtxitr).mcpe.isAncestorOf(mcp)) {
          //       nfoundtracks++;
          //     }
          //   }
          //}
          //if (nfoundtracks==1)     {
          //   std::cerr << " (1 pfo track for this vertex exists but was not associated.)";
          //   TVector3 vpos(mcvtxitr->mcpe.getAncestor()->getEndVertex());
          //   std::cerr << " TEST : (" << vpos.X() << ", " << vpos.Y() << ", " << vpos.Z() << ") "; 
          //} else if (nfoundtracks>1) {
          //   std::cerr << " (" << nfoundtracks << " pfo tracks for this vertex exist but were not associated.)";
          //   TVector3 vpos(mcvtxitr->mcpe.getAncestor()->getEndVertex());
          //   std::cerr << " TEST : (" << vpos.X() << ", " << vpos.Y() << ", " << vpos.Z() << ") "; 
          //} else                     std::cerr << " (No track exist.)";
        }
        std::cerr << std::endl;
    }
    std::cerr << "Vertex reconstruction efficiency = " << float(nbcvtxrec)/float(nbcvtxmc) 
              << " (" << float(nbcvtxrec) << "/" <<float(nbcvtxmc) << ")" << std::endl;
    std::cerr << "Vertex reconstruction efficiency (b) = " << float(nbcvtxrec_b)/float(nbcvtxmc_b) << std::endl; 
    std::cerr << "Vertex reconstruction efficiency (c) = " << float(nbcvtxrec_c)/float(nbcvtxmc_c) << std::endl;
  }
#endif
  _nEvt++;

  // clean up
  if (ipmc) delete ipmc;
  std::set<MCParticle*>::reverse_iterator ipbmbkgmcsitr;
  for (ipbmbkgmcsitr=BmbkgIPMCParticles.rbegin(); ipbmbkgmcsitr!=BmbkgIPMCParticles.rend();ipbmbkgmcsitr++) delete *ipbmbkgmcsitr; // delete only obejct newly created in this process.
}

void MyAnalysis::end() {
  _trktree->Write();
  _pvtxtree->Write();
  _svtxtree->Write();
  _file->Close();
  std::cerr << "end called." << std::endl;
}

bool MyAnalysis::updateSemistableAncestorOf(MCParticleExt* p, std::vector<int>& plist)
{
   if (p->getGeneration()<0) p->setGeneration(0);

   int nElements = plist.size();

   int pdg = p->getMCParticle()->getPDG();

   if (p->getGeneration()>0) {
     
     std::vector<int>::iterator itr;
     itr = find(plist.begin(), plist.end(), abs(pdg));
     if (itr!=plist.end()) {
       if (p->getMCParticle()->getDaughters().size()>0) {
         p->setAncestor(p->getMCParticle());
         return true;
       }
     }
   }

   p->incrementGeneration();
   if (p->getMCParticle()->getParent()) {
     const MCParticle* orig = p->getMCParticle();
     p->setMCParticle(p->getMCParticle()->getParent());
     bool retval = updateSemistableAncestorOf(p,plist);
     p->setMCParticle(orig); // set it back.
     if (retval) return retval;
   }
   return false;
}

bool MyAnalysis::updateSemistableDescendantOf(MCParticleExt* p, std::vector<int>& plist)
{
   if (p->getGeneration()>0) p->setGeneration(0);

   int nElements = plist.size();
   int pdg = p->getMCParticle()->getPDG();
//std::cerr << "    updateSemistableDescendantOf :: pdg = " << pdg << "  gen = " << p->getGeneration() << std::endl;

   if (p->getGeneration()<0) {
     
     std::vector<int>::iterator itr;
     itr = find(plist.begin(), plist.end(), abs(pdg));
     if (itr!=plist.end()) {
       p->setDescendant(p->getMCParticle());
       return true;
     }
   }

   p->decrementGeneration();
   for ( int i = 0; i < p->getMCParticle()->getDaughters().size(); i++) {
     const MCParticle* orig = p->getMCParticle();
     p->setMCParticle(p->getMCParticle()->getDaughters()[i]);
     bool retval = updateSemistableDescendantOf(p,plist);
     p->setMCParticle(orig); // set it back.
//if (retval) std::cerr << "    updateSemistableDescendantOf :: return pdg = " << p->getDescendant()->getPDG() << std::endl;
     if (retval) return retval;
   }
   return false;
}

const MCParticle* MyAnalysis::getParentIfSamePDG(const MCParticle* p)
{
    if (p->getParent()){
      const MCParticle* parent = p->getParent();
      if (parent->getPDG()==p->getPDG()) {
        return getParentIfSamePDG(parent);
      }
    }
    return p; 
}

//TVector3 MyAnalysis::getIPTruth(const MCParticle* p) {
TVector3 MyAnalysis::getIPTruth() {
  //const MCParticle* mcp = p;
  //while (mcp->getParent()) {
  //  mcp = mcp->getParent();
  //}
  //return mcp->getVertex();
  return (*_mcps)[0]->getVertex();
}

bool MyAnalysis::isFromPrimaryVertex(const MCParticle* p)
{
  const MCParticle* newp = getParentIfSamePDG(p); // if the particle radiates any particle, it shouldn't be recognized as its decay.
  TVector3 vtxvec(newp->getVertex());
  //TVector3 ipTruth = getIPTruth(p);
  TVector3 ipTruth = getIPTruth();
  return ((vtxvec-ipTruth).Mag() < aSmallNumber); 
}

bool MyAnalysis::isFromBeambackground(const MCParticle* p)
{
#if 0

#if 1
   const MCParticleVec& mcps = Event::Instance()->getMCParticles();
   if (p==mcps[0] || p==mcps[1]) return false;
#endif
   return true;
#else // flag added to MCParticle. See simulation flag in LCIO::MCParticle.
      // 23: overlay
   return (p->getFlag() & (1<<23));
#endif
}

const MCParticle* MyAnalysis::getOriginMCParticle(const MCParticle* p)
{
   if (p->getParent()) {
     const MCParticle *parent = p->getParent();
     return getOriginMCParticle(parent);
   }
   return p;
}

void MyAnalysis::checkVertices(std::vector<const Vertex*>& vtxs, std::vector<MCVertex>& mcvtxs, 
                                                                 std::map<const MCParticle*,MCParticleExt>& mcpMap,
                                                                 TTree* outtree)
{
   const bool dump = true;
   for (unsigned int iv = 0; iv < vtxs.size(); iv++) {
      std::vector<MCParticleExt> vtxmcps; // Candidates of the truth vertex for this vertex.
      std::vector<MCParticleExt>::iterator vtxmcpitr; // Candidates of the truth vertex for this vertex.

      std::cerr << "New Vertex :" << std::endl;
      const Vertex* vtx = vtxs[iv];
      TVector3 vpos(vtx->getPos());
      _vtxdata.rcx = vpos.X();
      _vtxdata.rcy = vpos.Y();
      _vtxdata.rcz = vpos.Z();
      const vector<const Track*> tracks = vtx->getTracks();
      _vtxdata.rcntrk = tracks.size();
      //std::cerr << "tracks.size() = " << tracks.size() << std::endl;
      std::map<const MCParticle*,MCParticleExt>::iterator mcpeitr;
      TLorentzVector vtx4p(0.,0.,0.,0.);
      TVector3 ipTruth;
      for (unsigned int itrk = 0; itrk < tracks.size(); itrk++) {
         const Track* trk = tracks[itrk];
         vtx4p += (*trk);
//std::cerr << "Test : getChi2Track() = " << vtx->getChi2Track(trk) << std::endl;
//Helix hel(&(*trk),PointBase::SECVTX);
//double tmin = -9999999;
//std::cerr << "loglikelihood = " << hel.LogLikelihood(vpos, tmin) << std::endl;

         const MCParticle* mcp = trk->getMcp();
         mcpeitr = mcpMap.find(mcp);
         if (mcpeitr != mcpMap.end()) {
            //MCParticleExt& mcpe = mcpeitr->second;
            MCParticleExt mcpe = mcpeitr->second;
            mcpe.addAssociatedTrack(trk);
#if 0
std::map<const MCParticle*,int>::iterator mcindexitr;
mcindexitr = _mcpIndex.find(mcp);
int index = -1;
if (mcindexitr != _mcpIndex.end()) index = mcindexitr->second;
std::cerr << "itrk = " << itrk << " mcp index = " << index << std::endl;
#endif
         
               //std::cerr <<  mcpe.getMCParticle()->getPDG() << " mcpe.getAncestor() = " << mcpe.getAncestor() << std::endl;

            if (itrk==0) {
               //ipTruth = getIPTruth(mcp);
               ipTruth = getIPTruth();
               _vtxdata.ipTruthx = ipTruth.X();
               _vtxdata.ipTruthy = ipTruth.Y();
               _vtxdata.ipTruthz = ipTruth.Z();
               if (dump) {
                 std::cerr << "      IP = (" << ipTruth.X() << ", " << ipTruth.Y() << ", " << ipTruth.Z() << ")" << std::endl;
                 std::cerr << "      Reconstructed vertex has " << tracks.size() << " track";
                 if (tracks.size()>1) std::cerr << "s";
                 std::cerr << ". ";
                 std::cerr << "(" << vpos.X() << ", " << vpos.Y() << ", " << vpos.Z() << ")"<< std::endl;
               }
            }

            //const MCParticle* origin = mcpe.getAncestor();
            //if (origin) std::cerr << "origin : " << origin->getPDG() << std::endl;
            //else        std::cerr << "origin 0, type : " << mcpe.getType() << std::endl;


            bool isFirstMCTruth = true; // check if the correspoinding truth candidate for this vertex track has been already registered in vtxmcps or not.
            // vtxmcps is a list of MC vtx truth candidates found so far. 
            // If same MC truth already found, add this track (trk) to the MC vtx truth (vtxmcptr).
            for (vtxmcpitr = vtxmcps.begin(); vtxmcpitr!=vtxmcps.end(); vtxmcpitr++) {

              //if (!mcpe.getAncestor()) { // tracks from beam bkg or primary

                if (mcpe.getType()==1 && vtxmcpitr->getType()==1) { //primary
                     vtxmcpitr->addAssociatedTrack(trk);
                     isFirstMCTruth = false;
                     break;
                } 
            
                if (mcpe.getType()==0 && vtxmcpitr->getType()==0) { //beam bkg
                    vtxmcpitr->addAssociatedTrack(trk);
                    isFirstMCTruth = false;
                    break;
                }
            
              //} else { // tracks from other than beam bkg or primary vtx.
                if (vtxmcpitr->getAncestor()==mcpe.getAncestor() ) { // if the same origin has already been found.
                   vtxmcpitr->addAssociatedTrack(trk);
                   isFirstMCTruth = false;
                   break;
                } 
              //}
            }

            // if this MC truth has not yet been registered, just add it to vtxmcps. 
            if (isFirstMCTruth) { 
#if 0
//std::cerr << "pdg = " << mcpe.getMCParticle()->getPDG() << " origin PDG = " << mcpe.getAncestor()->getPDG() << "(" << mcpe.getAncestor() << ")"<< std::endl;
std::cerr << "pdg = " << mcp->getPDG() << " origin PDG = " << mcpe.getAncestor()->getPDG() << "(" << mcpe.getAncestor() << ")"<< std::endl;
#endif
                // if it is not any semistable particle.
                //if (mcpe.getType()!=5) vtxmcps.push_back(mcpe);
                vtxmcps.push_back(mcpe);
            }
         }
      } //track loop
      _vtxdata.rcmass = vtx4p.M();
      _vtxdata.rce    = vtx4p.E();
      _vtxdata.rcp    = vtx4p.Vect().Mag();
      _vtxdata.rcpt   = vtx4p.Pt();
      _vtxdata.rcpx   = vtx4p.Px();
      _vtxdata.rcpy   = vtx4p.Py();
      _vtxdata.rcpz   = vtx4p.Pz();
      _vtxdata.rccostheta = _vtxdata.rcpz/_vtxdata.rcp;
      _vtxdata.xerr = TMath::Sqrt(vtx->getCov()[Vertex::xx]);
      _vtxdata.yerr = TMath::Sqrt(vtx->getCov()[Vertex::yy]);
      _vtxdata.zerr = TMath::Sqrt(vtx->getCov()[Vertex::zz]);
      _vtxdata.chi2 = vtx->getChi2();
      //for (int i = 0; i < vtxmcps.size(); i++) {
      //  std::cerr << " i = " << i << " vtx->getPos().X() = " << vtx->getPos().X() << std::endl;
      //  vtxmcps[i].setRecPos(vtx->getPos()); // This will be used when # of NAssociatedTracks are identical. 
      //}
      sort(vtxmcps.begin(),vtxmcps.end(),MCParticleExt::AscendingNAssociatedTracks);

      int ntotaltracks = 0;
      for (int i = 0; i < vtxmcps.size(); i++) {
        if (dump) std::cerr << "      Candidate origin " << i << ") ";
        TVector3 mcvpos;
        string origin;
        if (vtxmcps[i].getType()==5) {
          std::cerr << "This is an unexpected case. Please check." << std::endl;
          std::abort();
        } else if (vtxmcps[i].getAncestor()->getPDG()==100000000) {
          origin = "IP";
          mcvpos = TVector3(vtxmcps[i].getAncestor()->getVertex());
        } else if (vtxmcps[i].getAncestor()->getPDG()==100000001) {
          origin = "Beam bkg.";
          mcvpos = TVector3(vtxmcps[i].getAncestor()->getVertex());
        } else {
          origin = std::to_string(vtxmcps[i].getAncestor()->getPDG());
          mcvpos = TVector3(vtxmcps[i].getAncestor()->getEndVertex()); // Note that this getEndVertex is different from LCIO::MCParticle::getEndVertex()!
          //mcvpos = TVector3(vtxmcps[i].getAncestor()->getVertex());
        }
        std::map<const MCParticle*,int>::iterator mcindexitr;
        mcindexitr = _mcpIndex.find(vtxmcps[i].getAncestor());
        int index = -1;
        if (mcindexitr != _mcpIndex.end()) index = mcindexitr->second;
        if (dump) std::cerr << origin << " ";
        if (dump) std::cerr << "(" << vtxmcps[i].getAncestor() << "), index = " << index << ", ntrcks assigned = " << vtxmcps[i].getNAssociatedTracks() 
                            << ", Vertex (" << mcvpos.X() << ", " << mcvpos.Y() << ", " << mcvpos.Z() << ")"<< std::endl;
//std::cerr << "nDau = " << vtxmcps[i].getAncestor()->getDaughters().size() << " " << vtxmcps[i].getGeneration() << " type = " << vtxmcps[i].getType() << std::endl;
#if 0 // check track fit
          for (int itrk = 0; itrk < vtxmcps[i].getNAssociatedTracks(); itrk++) {
             const Track* trk = vtxmcps[i].getAssociatedTrack(itrk);
             Helix hel(trk,PointBase::SECVTX);
             double tmin;
             std::cerr << "track loglikelihood = " << hel.LogLikelihood(vpos,tmin) << " " << tmin<< std::endl;
          }
#endif
        ntotaltracks += vtxmcps[i].getNAssociatedTracks();
      }
 
      // initialization
      _vtxdata.weight = -1.;
      _vtxdata.mcpdg  = 0;

      TVector3 mcvpos;
      if (vtxmcps.size() ) {
        if (vtxmcps.size()>1) {
          int index_w_sameNtracks;
          for (index_w_sameNtracks = 1; index_w_sameNtracks < vtxmcps.size(); index_w_sameNtracks++) {
            if (vtxmcps[0].getNAssociatedTracks()!=vtxmcps[index_w_sameNtracks].getNAssociatedTracks()) break; 
          }
          index_w_sameNtracks;
        
//if (index_w_sameNtracks>0) std::cerr << "[0-" << index_w_sameNtracks << "] must be tested." << std::endl;
//else std::cerr << "only [0] must be tested." << std::endl;
//std::cerr << "vtxmcps[0].getAncestor()->getDaughters().size() = " << std::endl;
//std::cerr << vtxmcps[0].getAncestor()->getDaughters().size() << std::endl;
          if (index_w_sameNtracks>1) {
             // When the # of tracks are identical, take the one that is nearer to the Rconstructed vertex position.
             int index_best = 0;
             float distance_best = (vtxmcps[0].getAncestor()->getEndVertex() - vtx->getPos()).Mag();
//std::cerr << "original distance_best = " << distance_best << std::endl;
             for (int i = 1; i < index_w_sameNtracks; i++) {
               float distance_test = (vtxmcps[i].getAncestor()->getEndVertex() - vtx->getPos()).Mag(); 
//std::cerr << "test " << i << " distance_test = " << distance_test << std::endl;
               if (distance_test < distance_best) {
                 index_best = i;
                 distance_best = distance_test; 
               } 
             }
//if (index_best>0) {
//std::cerr << "Must be swapped. index_best = " << index_best << std::endl;
//std::cerr << "before swap : (vtxmcps[0].getAncestor()->getEndVertex() - vtx->getPos()).Mag() = " << (vtxmcps[0].getAncestor()->getEndVertex() - vtx->getPos()).Mag() << std::endl;
             if (index_best>0) std::iter_swap(vtxmcps.begin(),vtxmcps.begin()+index_best);  
//std::cerr << "after swap : (vtxmcps[0].getAncestor()->getEndVertex() - vtx->getPos()).Mag() = " << (vtxmcps[0].getAncestor()->getEndVertex() - vtx->getPos()).Mag() << std::endl;
//}
          }

        }
        string origin;
        if (vtxmcps[0].getAncestor()->getPDG()==100000000) {
          origin = "IP";
          mcvpos = TVector3(vtxmcps[0].getAncestor()->getVertex());
        } else if (vtxmcps[0].getAncestor()->getPDG()==100000001) {
          origin = "Beam bkg.";
          mcvpos = TVector3(vtxmcps[0].getAncestor()->getVertex());
        } else {
          origin = std::to_string(vtxmcps[0].getAncestor()->getPDG());
          mcvpos = TVector3(vtxmcps[0].getAncestor()->getEndVertex()); // Note that this getEndVertex is different from LCIO::MCParticle::getEndVertex()!
        } 
        _vtxdata.weight = float(vtxmcps[0].getNAssociatedTracks())/float(ntotaltracks);
        _vtxdata.mcpdg = vtxmcps[0].getAncestor()->getPDG();
        if (dump) cerr << "      Best estimation :   " << origin << " (" << vtxmcps[0].getAncestor() << "), w = " << _vtxdata.weight << endl; 
      }

      // loop in MCVertex to be found. vtxmcps[0] is the best candidate.
      for (int imcv = 0; imcv < mcvtxs.size(); imcv++) {
//std::cerr << "##### mcvtxs[imcv].mcpe.getMCParticle()->getPDG() = " << mcvtxs[imcv].mcpe.getMCParticle()->getPDG() << std::endl; 
        if (vtxmcps[0].getAncestor()->getDaughters().size()) { // secondary vertex
          //if (mcvtxs[imcv].mcpe->getMCParticle()==vtxmcps[0]->getMCParticle()) {
          if (mcvtxs[imcv].mcpe.getAncestor()==vtxmcps[0].getAncestor()) {
            if (mcvtxs[imcv].purity < float(vtxmcps[0].getNAssociatedTracks())/float(ntotaltracks)) {
              mcvtxs[imcv].matchedvtxrec = vtx; // correctly found 
              mcvtxs[imcv].ntrk = ntotaltracks; 
              mcvtxs[imcv].purity = float(vtxmcps[0].getNAssociatedTracks())/float(ntotaltracks); 
            }
          }
        } else { // primary vertex, beam bkg.
            //cerr << "##### primary vertex!" << std::endl;
           if (mcvtxs[imcv].mcpe.getType() == vtxmcps[0].getType()) {
              if (mcvtxs[imcv].purity < float(vtxmcps[0].getNAssociatedTracks())/float(ntotaltracks)) {
                mcvtxs[imcv].mcpe = vtxmcps[0];
                mcvtxs[imcv].matchedvtxrec = vtx; // correctly found 
                mcvtxs[imcv].ntrk = ntotaltracks; 
                mcvtxs[imcv].purity = float(vtxmcps[0].getNAssociatedTracks())/float(ntotaltracks); 
              }
           }
        }
      }

      _vtxdata.type  = vtxmcps[0].getType();
      _vtxdata.mcx   = mcvpos.X();
      _vtxdata.mcy   = mcvpos.Y();
      _vtxdata.mcz   = mcvpos.Z();
      _vtxdata.mcpt  = vtxmcps[0].getAncestor()->Pt();
      _vtxdata.mcpx  = vtxmcps[0].getAncestor()->Px();
      _vtxdata.mcpy  = vtxmcps[0].getAncestor()->Py();
      _vtxdata.mcpz  = vtxmcps[0].getAncestor()->Pz();
      _vtxdata.mcp   = vtxmcps[0].getAncestor()->Vect().Mag();
      _vtxdata.mce   = vtxmcps[0].getAncestor()->E();
      _vtxdata.mccostheta = _vtxdata.mcpz/_vtxdata.mcp;
      _vtxdata.mcpdg = vtxmcps[0].getAncestor()->getPDG();
      _vtxdata.mcmass= vtxmcps[0].getAncestor()->M();
      _vtxdata.mcntrk=0; 
      for ( int idau = 0; idau < vtxmcps[0].getAncestor()->getDaughters().size(); idau++ ) {
         const MCParticle* dau = vtxmcps[0].getAncestor()->getDaughters()[idau];
         if (fabs(dau->getCharge())>0) {
           TVector3 start(dau->getVertex());
           TVector3 end(dau->getEndVertex());
           if ((end-start).Mag()>30) {
           //if ((end-start).Mag()>0.01) {
              _vtxdata.mcntrk++; // Assume 3cm-long track can be found (~ 1cm^3 voxel * 3 hits). 
              //std::cerr << "dau->getPDG() (" << dau << ") = " << dau->getPDG() << std::endl;
           }
         }
      }

      TVector3 axisL = (mcvpos - ipTruth).Unit(); // Axis along pvtx - svtx.
//if (axisL.Mag()<1) std::cerr << "axisL.Mag() = " << axisL.Mag() << std::endl;
      TVector3 axisZ(0.,0.,1.);
      TVector3 axisT1 = axisZ.Cross(axisL);
      TVector3 axisT2 = axisL.Cross(axisT1);

      TVector3 residualvec = vpos - mcvpos;
      float dl = residualvec.Dot(axisL);
      //TVector3 residualL = dl*axisL;
      //TVector3 residualT = residualvec - residualL;
      //float dt1 = residualT.Dot(axisT1);
      //float dt2 = residualT.Dot(axisT2);
      float dt1 = residualvec.Dot(axisT1);
      float dt2 = residualvec.Dot(axisT2);
      _vtxdata.dl = dl;
      _vtxdata.dt1 = dt1; // z x L
      _vtxdata.dt2 = dt2; // L x T1
      _vtxdata.distancePvtxToSvtx = (mcvpos - ipTruth).Mag();

      outtree->Fill();
   }//vertex loop
}

//void MyAnalysis::setTracksOf(MCParticle* vtx, const MCParticleVec& mcps) {
//void MyAnalysis::setTracksOf(MCParticle* vtx) {
//  TVector3 vpos(vtx->getVertex());
//  for (unsigned int i=0; i<_mcps->size(); i++) {
//    TVector3 pos((*_mcps)[i]->getVertex());
//    if ((vpos-pos).Mag() < aSmallNumber) vtx->addDaughter((*_mcps)[i]); // set tracks if the origin is near enough.
//  }
//}

MCParticle* MyAnalysis::createVtxMCParticle(TVector3& vpos) {

  MCParticle* vtx = new MCParticle;
  vtx->setVertex(vpos);
  std::set<const MCParticle*> dau_candidates; // avoid double-counting.
  for (unsigned int i=0; i<_mcps->size(); i++) {
    TVector3 pos(getParentIfSamePDG((*_mcps)[i])->getVertex());
    TVector3 epos((*_mcps)[i]->getEndVertex()); // Note that if the particle has no daughter, the end vertex will be identical to the start vertex.
//if (fabs(vpos.X())<aSmallNumber&&(fabs(vpos.Y())<aSmallNumber) {
//}
    if ((vpos-pos).Mag() < aSmallNumber && ((epos-pos).Mag() > aSmallNumber || (*_mcps)[i]->getDaughters().size()==0) ) { // require the daughter fly somewhat (aSmallNumber) or no daughters.
      dau_candidates.insert(getParentIfSamePDG((*_mcps)[i]));
//    std::cerr << "i = " << i << " " << (*_mcps)[i]->getPDG() << " (vpos-pos).Mag() = " << (vpos-pos).Mag() << " (epos-pos).Mag() = " << (epos-pos).Mag() << std::endl;
//if (fabs(vpos.X())<aSmallNumber) std::cerr << " added" << std::endl;
    }
  }
bool test = false;
//std::cerr << std::endl;
  std::set<const MCParticle*>::iterator dauitr;
  for (dauitr = dau_candidates.begin(); dauitr != dau_candidates.end(); dauitr++) {
    vtx->addDaughter(*dauitr); // set tracks if the origin is near enough.
    *vtx += *(const_cast<MCParticle*>(*dauitr)); 
//if ((*dauitr)->Pz()>40) {
//std::cerr << "#### pz = " << (*dauitr)->Pz();
TVector3 pos((*dauitr)->getVertex());
TVector3 epos((*dauitr)->getEndVertex());
std::map<const MCParticle*,int>::iterator mcindexitr;
mcindexitr = _mcpIndex.find(*dauitr);
int index = -1;
if (mcindexitr != _mcpIndex.end()) index = mcindexitr->second;
//std::cerr << "  "<< vpos.Z() - pos.Z() << " flight len = " << (pos-epos).Mag() << " pdg = " << (*dauitr)->getPDG() << " index = " << index << std::endl;
//test = true;
//}
  }
//if (test) std::cerr << " pz sum = " << vtx->Pz() << std::endl;
//std::cerr << " pz sum = " << vtx->Pz() << std::endl;
//std::cerr << std::endl;
  return vtx; 
// why sum of momutum of primary vertex tracks is not 0? It seems to be related to ISR emission.
} 

//bool MyAnalysis::MCParticleExt::isAncestorOf(const MCParticle* in) {
//  if (!in) return false;
//  if (_ancestor==in) return true;
//  else {
//     const MCParticle* parent = in->getParent();
//     if(isAncestorOf(parent)) return true; 
//  }
//  return false;
//};
//
//bool MyAnalysis::MCParticleExt::isDescendantOf(const MCParticle* in) {
//  if (!in) return false;
//  if (_descendant==in) return true;
//  else {
//    if (in->getDaughters().size()) {
//      for (int i = 0; i < in->getDaughters().size(); i++) {
//        const MCParticle* daughter = in->getDaughters()[i];
//        if(isDescendantOf(daughter)) return true; 
//      }
//    }
//  }
//  return false;
//};
#endif

}
