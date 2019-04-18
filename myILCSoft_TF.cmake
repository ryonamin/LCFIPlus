# Original : /cvmfs/ilc.desy.de/sw/x86_64_gcc49_sl6/v02-00-02/init_ilcsoft.sh
# We need a private Python to install tensorflow.
# MarlinUtil and ROOT are also to be recompiled accordingly.

################################################################################
# Environment script generated by ilcsoft-install on Tue Sep  4 11:27:01 2018
# for the installation located at [ /cvmfs/ilc.desy.de/sw/x86_64_gcc49_sl6/v02-00-02 ]
################################################################################

SET( ILC_HOME "/cvmfs/ilc.desy.de/sw/x86_64_gcc49_sl6/v02-00-02" CACHE PATH "Path to ILC Software" FORCE)
MARK_AS_ADVANCED( ILC_HOME )

SET( CMAKE_PREFIX_PATH 
	${ILC_HOME}/lccd/v01-05;
	${ILC_HOME}/Marlin/v01-16;
	${ILC_HOME}/MarlinDD4hep/v00-06;
	${ILC_HOME}/DDMarlinPandora/v00-10;
	${ILC_HOME}/MarlinReco/v01-25;
        /home/ilc/yonamine/work/LCFIPlusDev/MarlinUtil
	${ILC_HOME}/PandoraAnalysis/v02-00-00;
	${ILC_HOME}/PandoraPFANew/v03-09-00;
	${ILC_HOME}/LCFIVertex/v00-07-04;
	${ILC_HOME}/CEDViewer/v01-16;
	${ILC_HOME}/Overlay/v00-21;
	${ILC_HOME}/MarlinFastJet/v00-05-01;
	${ILC_HOME}/LCTuple/v01-11;
	${ILC_HOME}/MarlinKinfit/v00-06;
	${ILC_HOME}/MarlinTrk/v02-07;
	${ILC_HOME}/KiTrack/v01-09;
	${ILC_HOME}/KiTrackMarlin/v01-12;
	${ILC_HOME}/MarlinTrkProcessors/v02-10;
	${ILC_HOME}/MarlinKinfitProcessors/v00-04;
	${ILC_HOME}/ILDPerformance/v01-06;
	${ILC_HOME}/Clupatra/v01-03;
	${ILC_HOME}/Physsim/v00-04-01;
	${ILC_HOME}/LCFIPlus/v00-06-09;
	${ILC_HOME}/FCalClusterer/v01-00;
	${ILC_HOME}/ForwardTracking/v01-13;
	${ILC_HOME}/ConformalTracking/v01-07;
	${ILC_HOME}/LICH/v00-01;
	${ILC_HOME}/pathfinder/v00-06-01;
	${ILC_HOME}/MarlinTPC/v01-04;
	${ILC_HOME}/bbq/v00-01-03;
	${ILC_HOME}/Garlic/v03-01;
	${ILC_HOME}/RAIDA/v01-09;
	${ILC_HOME}/KalTest/v02-05;
	${ILC_HOME}/KalDet/v01-14-01;
	${ILC_HOME}/GBL/V02-01-01;
	${ILC_HOME}/xercesc/3.1.4;
	${ILC_HOME}/DD4hep/v01-07-02;
	${ILC_HOME}/lcgeo/v00-16-03;
	${ILC_HOME}/aidaTT/v00-09;
	${ILC_HOME}/DDKalTest/v01-05;
	${ILC_HOME}/DD4hepExamples/v01-07-02;
	${ILC_HOME}/CED/v01-09-02;
	${ILC_HOME}/lcio/v02-12-01;
	${ILC_HOME}/gear/v01-08;
	${ILC_HOME}/FastJet/3.2.1;
	/home/ilc/yonamine/work/root-6.08.06;
	${ILC_HOME}/CLHEP/2.3.4.3;
	${ILC_HOME}/gsl/2.1;
	${ILC_HOME}/QT/4.7.4;
	${ILC_HOME}/geant4/10.03.p02;
	${ILC_HOME}/CondDBMySQL/CondDBMySQL_ILC-0-9-6;
	${ILC_HOME}/mysql/5.0.45;
	${ILC_HOME}/ilcutil/v01-05;
	/cvmfs/ilc.desy.de/sw/boost/1.58.0;
	/cvmfs/ilc.desy.de/sw/Eigen/3.2.9;
CACHE PATH "CMAKE_PREFIX_PATH" FORCE )


option(USE_CXX11 "Use cxx11" True)
option(Boost_NO_BOOST_CMAKE "dont use cmake find module for boost" ON)
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g" CACHE STRING "" FORCE )

