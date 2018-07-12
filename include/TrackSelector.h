// TrackSelector.h

#ifndef TrackSelector_h
#define TrackSelector_h 1

#include "lcfiplus.h"
#include "algoSigProb.h"
#include <vector>
#include <cmath>

using namespace std;
using namespace lcfiplus;

namespace lcfiplus {

class TrackSelectorConfig {
 public:
  // cuts which are combined using the AND scheme
  double minD0;
  double maxD0;
  double minD0Err;
  double maxD0Err;
  double minD0Sig;
  double maxD0Sig;
  double minZ0;
  double maxZ0;
  double minZ0Err;
  double maxZ0Err;
  double minZ0Sig;
  double maxZ0Sig;
  double minD0Z0Sig;
  double maxD0Z0Sig;
  double minPt;
  double maxInnermostHitRadius;
  // cuts which are combined using the OR scheme, then AND'd with the AND schemes above
  int minTpcHits;
  double minTpcHitsMinPt;
  int minFtdHits;
  int minVtxHits;
  int minVtxPlusFtdHits;

  TrackSelectorConfig() {
    minD0 = 0.;
    maxD0 = 1e+300;
    minD0Err = 0.;
    maxD0Err = 1e+300;
    minD0Sig = 0.;
    maxD0Sig = 1e+300;
    minZ0 = 0.;
    maxZ0 = 1e+300;
    minZ0Err = 0.;
    maxZ0Err = 1e+300;
    minZ0Sig = 0.;
    maxZ0Sig = 1e+300;
    minD0Z0Sig = 0.;
    maxD0Z0Sig = 1e+300;
    minPt = 0.;
    maxInnermostHitRadius = 1e+300;

    minTpcHits = 999999;
    minTpcHitsMinPt = 999999;
    minFtdHits = 999999;
    minVtxHits = 999999;
    minVtxPlusFtdHits = 0;
  }
};

class TrackSelector {
 public:
#if 0
  vector<const Track*> operator () (const vector<const Track*>& tracks, TrackSelectorConfig& config) {
#else
  vector<const Track*> operator () (const vector<const Track*>& tracks, TrackSelectorConfig& config, const Vertex* ip = 0) {
#endif
    vector<const Track*> ret;

    for (unsigned int i=0; i<tracks.size(); i++) {
#if 0
      if (passesCut(tracks[i], config))
#else
      if (passesCut(tracks[i], config, ip))
#endif
        ret.push_back(tracks[i]);
    }

    return ret;
  }

#if 0
  bool passesCut(const Track* trk, const TrackSelectorConfig& cfg) {
    // AND cuts
    
    if (fabs(trk->getD0()) < cfg.minD0) return false;
    if (fabs(trk->getD0()) > cfg.maxD0) return false;
    double d0sig = fabs(trk->getD0()) / sqrt(trk->getCovMatrix()[tpar::d0d0]);
    if (fabs(trk->getZ0()) < cfg.minZ0) return false;
    if (fabs(trk->getZ0()) > cfg.maxZ0) return false;
    double z0sig = fabs(trk->getZ0()) / sqrt(trk->getCovMatrix()[tpar::z0z0]);
    if (trk->getCovMatrix()[tpar::d0d0] < cfg.minD0Err) return false;
    if (trk->getCovMatrix()[tpar::d0d0] > cfg.maxD0Err) return false;

#else
  bool passesCut(const Track* trk, const TrackSelectorConfig& cfg, const Vertex* ip = 0) {
    if (!ip) ip = new Vertex();
    // AND cuts
//std::cerr << "trk->getD0() = " << trk->getD0() << std::endl;
    //double trkx = - trk->getD0() * sin(trk->getPhi());
    //double trky = trk->getD0() * cos(trk->getPhi());
    //double trkz = trk->getZ0();
    //TVector3 trkv(trkx,trky,trkz);
    //TVector3 ipv = ip->getPos();

    bool updateFlt = true;
    TVector3 pca = algoSigProb::trackPositionFromPrimaryVertex(trk, ip, updateFlt);
    double d0 = pca.Perp();
    double z0 = pca.Z();
    
//std::cerr << "trkx = " << trkx << " ipvx = " << ipv.X() << std::endl;
//std::cerr << "trky = " << trky << " ipvy = " << ipv.Y() << std::endl;
//std::cerr << "trkz = " << trkz << " ipvz = " << ipv.Z() << std::endl;
    //if ((trkv-ipv).Perp() < cfg.minD0) return false;
    //if ((trkv-ipv).Perp() > cfg.maxD0) return false;
    if (d0 < cfg.minD0) return false;
    if (d0 > cfg.maxD0) return false;
    //double d0sig = fabs(d0) / sqrt(trk->getCovMatrix()[tpar::d0d0]);
    // error from vertex 
    double x0 = ip->getX();
    double y0 = ip->getY();
    double x02 = x0*x0;
    double y02 = y0*y0;
    // D0^2 = x^2 + y^2
    double pvcov = (ip->getCov()[Vertex::xx]*x02 + ip->getCov()[Vertex::yy]*y02)/(x02+y02);
    double d0err = sqrt(trk->getCovMatrix()[tpar::d0d0]+pvcov);
    double d0sig = fabs(d0) / d0err;
    //if (abs((trkv-ipv).Z()) < cfg.minZ0) return false;
    //if (abs((trkv-ipv).Z()) > cfg.maxZ0) return false;
    if ( abs(z0) < cfg.minZ0) return false;
    if ( abs(z0) > cfg.maxZ0) return false;
    double z0err = sqrt(trk->getCovMatrix()[tpar::z0z0]);
    double z0sig = fabs(z0) / z0err;
    //if (trk->getCovMatrix()[tpar::d0d0] < cfg.minD0Err) return false;
    //if (trk->getCovMatrix()[tpar::d0d0] > cfg.maxD0Err) return false;
    if ( d0err < cfg.minD0Err) return false;
    if ( d0err > cfg.maxD0Err) return false;

#endif

    if ( d0sig < cfg.minD0Sig) return false;
    if ( d0sig > cfg.maxD0Sig) return false;

#if 0
    if (trk->getCovMatrix()[tpar::z0z0] < cfg.minZ0Err) return false;
    if (trk->getCovMatrix()[tpar::z0z0] > cfg.maxZ0Err) return false;
#else
    if ( z0err < cfg.minZ0Err) return false;
    if ( z0err > cfg.maxZ0Err) return false;
#endif
    if ( z0sig < cfg.minZ0Sig) return false;
    if ( z0sig > cfg.maxZ0Sig) return false;

    if (sqrt(d0sig * d0sig + z0sig * z0sig) < cfg.minD0Z0Sig)return false;
    if (sqrt(d0sig * d0sig + z0sig * z0sig) > cfg.maxD0Z0Sig)return false;

    if (trk->Pt() < cfg.minPt) return false;
    if (trk->getRadiusOfInnermostHit() > cfg.maxInnermostHitRadius) return false;

    // OR cuts
    if (trk->getFtdHits() >= cfg.minFtdHits) return true;
    if (trk->getVtxHits() >= cfg.minVtxHits) return true;
    if (trk->getVtxHits() + trk->getFtdHits() >= cfg.minVtxPlusFtdHits) return true;
    if (trk->getTpcHits() >= cfg.minTpcHits && trk->Pt() > cfg.minTpcHitsMinPt) return true;

    return false;
  }

  //c-tor / d-tor
  TrackSelector() {}
  ~TrackSelector() {}
};
}

#endif //TrackSelector_h
